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
save(classifs, file = 'classifs.RData')

# match filtering in taxonomy
# easier to just run the classifs processing loop again with the new genome table
taxonomy = ldply(1:length(classifs), function(i) {
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

# remove duplicates, order / enforce columns
df = data.frame(name = character(), rank = character(), taxid = character(), domain = character(), kingdom = character(), subkingdom = character(), phylum = character(), subphylum = character(), class = character(), subclass = character(), order = character(), family = character(), genus = character(), species = character(), assembly = character(), stringsAsFactors = FALSE)
df = bind_rows(df, taxonomy %>% replace(., is.na(.), '') %>% unique)
taxonomy = df %>% replace(., is.na(.), '')

write.table(taxonomy, file = 'taxonomy.tsv', quote = FALSE, sep = '\t', row.names = FALSE)


# add classifications back to original genome table
genome_table = genome_table %>% left_join(taxonomy %>% filter(assembly != ''), by = 'taxid') %>%
  select(dbname, taxid, name, domain, kingdom, subkingdom, phylum, subphylum, class, subclass, order, family, genus, species, assembly)

# filter out entries for which files are missing
# add modifier to filename for higher eukaryotes (some include "-detailed" for out files)
genome_table = genome_table %>%
  filter(file.exists(paste0('data/ss/', dbname, '-tRNAs.ss')),
    file.exists(paste0('data/out/', dbname, '-tRNAs.out')) | file.exists(paste0('data/out/', dbname, '-tRNAs-detailed.out')),
    file.exists(paste0('data/iso/', dbname, '-tRNAs.iso')))

# Filter out anything that isn't "representative" or "reference"
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
# assembly  species  name
# 559292  4932  Saccharomyces cerevisiae S288c

approved_taxids = c('559292')

ncbi_refs = read.delim('/projects/lowelab/db/ncbi/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt', stringsAsFactors = FALSE, sep = '\t', comment.char = '#', header = FALSE) %>% 
  filter(V5 %in% c('representative genome', 'reference genome'))
genome_table = genome_table %>% filter(taxid %in% ncbi_refs$V6 | taxid %in% ncbi_refs$V7 | taxid %in% approved_taxids)

write.table(genome_table, file = 'genomes.tsv', quote = FALSE, sep = '\t', row.names = FALSE)