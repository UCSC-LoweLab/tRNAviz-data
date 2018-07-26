#!/usr/bin/env Rscript
library(taxize)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)

genome_table = read.delim('taxids.tsv', col.names = c('dirname', 'dbname', 'taxid'), stringsAsFactors = FALSE) %>%
  mutate(taxid = as.character(taxid))
fungi_taxids = genome_table$taxid %>% unlist %>% unname %>% unique

# grab classifications for each genome assembly
classifs = classification(fungi_taxids, db = "ncbi")

# parse classifications to usable format
taxonomy = ldply(1:length(classifs), function(i) {
  taxid = names(classifs)[i]
  classif = classifs[[i]]
	classif %>% 
    filter(rank %in% c('superkingdom', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus', 'species', 'varietas')) %>%
    select(name, rank) %>%
    rbind(c(name = taxid, rank = 'taxid')) %>%
    rbind(c(name = classif[classif$id == taxid, ]$name, rank = 'assembly')) %>%
    spread(rank, name, fill = ' ')
})

taxonomy = taxonomy %>%
  rename(domain = superkingdom) %>% 
  replace_na(list(domain = '', kingdom = '', subkingdom = '', phylum = '', subphylum = '', class = '', subclass = '', order = '', family = '', genus = '', species = '', varietas = '', assembly = ''))

# add classifications back to original genome table
genome_table = genome_table %>% left_join(taxonomy, by = 'taxid') %>%
  select(dbname, dirname, taxid, assembly, varietas, species, genus, family, order, subclass, class, subphylum, phylum, subkingdom, kingdom, domain)

write.table(genome_table, file = 'genomes.tsv', quote = FALSE, sep = '\t', row.names = FALSE)


# grab classifications for each genome assembly again, but save all taxid info
taxonomy_long = ldply(1:length(classifs), function(i) {
  taxid = names(classifs)[i]
  classif = classifs[[i]]
  classif %>% 
    filter(rank %in% c('superkingdom', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus', 'species', 'varietas')) %>%
    rename(taxid = id) %>%
    rbind(c(name = classif[classif$id == taxid, ]$name, rank = 'assembly', taxid = taxid))
})

# remove duplicates and same-name-assembly situations
taxonomy_long = taxonomy_long %>% unique

write.table(taxonomy_long, file = 'taxonomy.tsv', quote = FALSE, sep = '\t', row.names = FALSE)