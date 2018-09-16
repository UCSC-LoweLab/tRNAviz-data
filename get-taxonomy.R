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
                      'varietas' = character(0),
                      stringsAsFactors = FALSE) %>% 
  rbind.fill(ldply(1:length(classifs), function(i) {
  taxid = names(classifs)[i]
  classif = classifs[[i]]
  if (taxid %in% classif$id) assembly = classif[classif$id == taxid, ]$name
  else assembly = classif[classif$rank == 'species', ]$name
  classif %>% 
    filter(rank %in% c('superkingdom', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus', 'species', 'varietas')) %>%
    select(name, rank) %>%
    rbind(c(name = taxid, rank = 'taxid')) %>%
    rbind(c(name = assembly, rank = 'assembly')) %>%
    spread(rank, name, fill = '')
}))

taxonomy = taxonomy %>% rename(domain = superkingdom)
taxonomy[is.na(taxonomy)] = ''

# add classifications back to original genome table
genome_table = genome_table %>% left_join(taxonomy, by = 'taxid') %>%
  select(dbname, dirname, taxid, assembly, varietas, species, genus, family, order, subclass, class, subphylum, phylum, subkingdom, kingdom, domain)

# filter out entries for which files are missing
# add modifier to filename for higher eukaryotes (some include "-detailed" for out files)
genome_table = genome_table %>%
  filter(file.exists(paste0('data/ss/', dbname, '-tRNAs.ss')),
    file.exists(paste0('data/out/', dbname, '-tRNAs.out')) | file.exists(paste0('data/out/', dbname, '-tRNAs-detailed.out')),
    file.exists(paste0('data/iso/', dbname, '-tRNAs.iso')))

write.table(genome_table, file = 'genomes.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

# grab classifications for each genome assembly again, but save all taxid info
taxonomy_long = ldply(1:length(classifs), function(i) {
  taxid = names(classifs)[i]
  if (!(taxid %in% genome_table$taxid)) return(data.frame()) # filter out entries for which files are missing
  classif = classifs[[i]]
  if (taxid %in% classif$id) assembly = classif[classif$id == taxid, ]$name
  else assembly = classif[classif$rank == 'species', ]$name
  classif %>% 
    filter(rank %in% c('superkingdom', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus', 'species', 'varietas')) %>%
    rename(taxid = id) %>%
    mutate(rank = ifelse(rank == 'superkingdom', 'domain', rank)) %>%
    rbind(c(name = assembly, rank = 'assembly', taxid = taxid))
})

# remove duplicates and same-name-assembly situations
taxonomy_long = taxonomy_long %>% unique

write.table(taxonomy_long, file = 'taxonomy.tsv', quote = FALSE, sep = '\t', row.names = FALSE)
