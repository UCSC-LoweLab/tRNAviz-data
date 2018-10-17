# This script is best run in a REPL using appropriate steps 
library(taxize)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)

genome_table = read.delim('taxids.tsv', col.names = c('dirname', 'dbname', 'taxid'), stringsAsFactors = FALSE) %>%
  mutate(taxid = as.character(taxid))
taxids = genome_table$taxid %>% unlist %>% unname %>% unique

# grab classifications for each genome assembly
classifs = classification(taxids, db = "ncbi")
save.image('taxonomy.RData')

# parse classifications to usable format
taxonomy = data.frame('superkingdom' = character(0), 
                      'kingdom' = character(0), 
                      'subkingdom' = character(0), 
                      'phylum' = character(0), 
                      'subphylum' = character(0), 
                      'class' = character(0), 
                      'subclass' = character(0), 
                      'order' = character(0), 
                      'family' = character(0), 
                      'genus' = character(0), 
                      'species' = character(0),
                      stringsAsFactors = FALSE) %>% 
  rbind.fill(ldply(1:length(classifs), function(i) {
  taxid = names(classifs)[i]
  classif = classifs[[i]]
  if (taxid %in% classif$id) assembly = classif[classif$id == taxid, ]$name
  else assembly = classif[classif$rank == 'species', ]$name
  classif %>% 
    filter(rank %in% c('superkingdom', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus', 'species')) %>%
    select(name, rank) %>%
    rbind(c(name = taxid, rank = 'taxid')) %>%
    rbind(c(name = assembly, rank = 'assembly')) %>%
    spread(rank, name, fill = '')
}))

taxonomy = taxonomy %>% rename(domain = superkingdom)
taxonomy[is.na(taxonomy)] = ''

# add classifications back to original genome table
genome_table = genome_table %>% left_join(taxonomy, by = 'taxid') %>%
  select(dbname, dirname, taxid, assembly, species, genus, family, order, subclass, class, subphylum, phylum, subkingdom, kingdom, domain)

# filter out entries for which files are missing
# add modifier to filename for higher eukaryotes (some include "-detailed" for out files)
genome_table = genome_table %>%
  filter(file.exists(paste0('data/ss/', dbname, '-tRNAs.ss')),
    file.exists(paste0('data/out/', dbname, '-tRNAs.out')) | file.exists(paste0('data/out/', dbname, '-tRNAs-detailed.out')),
    file.exists(paste0('data/iso/', dbname, '-tRNAs.iso')))


# grab classifications for each genome assembly again, but save all taxid info
taxonomy_long = ldply(1:length(classifs), function(i) {
  taxid = names(classifs)[i]
  if (!(taxid %in% genome_table$taxid)) return(data.frame()) # filter out entries for which files are missing
  classif = classifs[[i]]
  if (taxid %in% classif$id) assembly = classif[classif$id == taxid, ]$name
  else assembly = classif[classif$rank == 'species', ]$name
  classif %>% 
    filter(rank %in% c('superkingdom', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus', 'species')) %>%
    rename(taxid = id) %>%
    mutate(rank = ifelse(rank == 'superkingdom', 'domain', rank)) %>%
    rbind(c(name = assembly, rank = 'assembly', taxid = taxid))
})

# remove duplicates and same-name-assembly situations
taxonomy_long = taxonomy_long %>% unique


# then filter out anything that isn't "representative" or "reference"
# stats:
# 3344 assemblies are represented by 1481 representative species. 
# 890 assemblies belong to 687 unrepresented species distributed as such, where "single" means these species have only 1 corresponding assembly.
#      domain  single  multiple
# *     <chr>   <dbl>     <dbl>
# 1   Archaea     140        15
# 2  Bacteria     468        62
# 3 Eukaryota       2         0

# in bacteria, 1059 species are represented, but we only have 1016 assemblies that actually are the representative one, with 89 assemblies orphaned (these will be treated as if they didn't have representative species)
# this discrepancy is because some assemblies for the same species are also considered "representative" (NCBI does this by strain, not species)


# Unfortuantely, this filters out some species. Here's a list of them:
# taxid  species_taxid  name
# 559292  4932  Saccharomyces cerevisiae S288c
# These will have to be manually edited in taxonomy.tsv and genomes.tsv later.

approved_taxids = c('4932')

ncbi_refs = read.delim('/projects/lowelab/db/ncbi/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt', stringsAsFactors = FALSE, sep = '\t', comment.char = '#', header = FALSE) %>% 
  filter(V5 %in% c('representative genome', 'reference genome'))
genome_table = genome_table[taxonomy_long[match(genome_table$species, taxonomy_long$name), ]$taxid %in% ncbi_refs$V7, ] %>%
  filter(taxid %in% ncbi_refs$V6 | taxid %in% approved_taxids)

write.table(genome_table, file = 'genomes.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

# match filtering in taxonomy_all
# easier to just run the classifs processing loop again with the new genome table
taxonomy_all = ldply(1:length(classifs), function(i) {
  taxid = names(classifs)[i]
  if (!(taxid %in% genome_table$taxid)) return(data.frame()) # filter out entries for which files are missing
  classif_set = classifs[[i]]
  if (taxid %in% classif_set$id) assembly = classif_set[classif_set$id == taxid, ]$name
  else assembly = classif_set[classif_set$rank == 'species', ]$name

  classif_set = classifs[[i]] %>% 
    filter(rank %in% c('superkingdom', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus', 'species')) %>%
    rename(taxid = id) %>%
    mutate(rank = ifelse(rank == 'superkingdom', 'domain', rank)) %>%
    rbind(c(name = assembly, rank = 'assembly', taxid = taxid))

  ldply(1:nrow(classif_set), function(j) {
    row = classif_set[j, ]
    parent_taxids = classif_set %>% mutate(taxid = ifelse(row_number() <= which(classif_set$rank == row$rank), taxid, '')) %>%
      select(rank, taxid) %>%
      spread(rank, taxid, fill = '')
    row %>% cbind(parent_taxids)
  })
})

# remove duplicates and order columns
taxonomy_all = taxonomy_all %>% 
  replace(., is.na(.), '') %>%
  unique %>%
  select(name, rank, taxid, kingdom, subkingdom, phylum, subphylum, class, subclass, order, family, genus, species, assembly))

# add missing genomes
taxonomy_all = taxonomy_all %>% rbind(
  taxonomy_all %>% filter(taxid == '4932') %>%
    filter(row_number() == 1) %>%
    mutate(name = 'Saccharomyces cerevisiae S288c', rank = 'assembly', taxid = '559292')
)

write.table(taxonomy_all, file = 'taxonomy.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

