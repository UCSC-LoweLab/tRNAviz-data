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

# Unfortuantely, this filters out some species. Here we add them back in. It's just yeast plus 20 Thaumarcheaota manually entered from RefSeq.
# species assembly name
# 4932  559292  Saccharomyces cerevisiae S288c
# 1088740  1001994  Nitrosarchaeum koreense MY1 (archaea)
# 1007084   859192  Candidatus Nitrosoarchaeum limnia BG20 (archaea)
# 1170320   859350  Candidatus Nitrosopumilus salaria BD31 (archaea)
# 1027373  1027373  Nitrosopumilus sp. AR (archaea)
# 1027374  1027374  Nitrosopumilus sp. SJ (archaea)
# 1353246  1353246  Candidatus Nitrosotenuis chungbukensis (archaea)
# 1034015   926571  Nitrososphaera viennensis EN76 (archaea)
# 1407055  1407055  Thaumarchaeota archaeon N4 (archaea)
# 1410606  1410606  Candidatus Nitrosopelagicus brevis (archaea)
# 1776294  1776294  Nitrosopumilus sp. Nsub (archaea)
# 1898749  1898749  Candidatus Nitrosomarinus catalina (archaea)
# 1846278  1846278  Candidatus Nitrosotenuis aquarius (archaea)
# 2045011  2045011  Candidatus Nitrosocaldus islandicus (archaea)
# 1410606  1410606  Candidatus Nitrosopelagicus brevis (archaea)
# 718286   718286 Candidatus Nitrosopumilus sp. NM25 (archaea)
# 1465461  1465461  Thaumarchaeota archaeon SCGC AC-337_F14 (archaea)
# 1499975  1499975  Nitrosotalea sp. Nd2 (archaea)
# 1903277  1903277  Nitrosotalea sp. SbT1 (archaea)

# we want to keep 3 additional archaeal genomes already in the set of tRNAscan-SE runs:
# 1229909  1229909  Candidatus Nitrosopumilus sediminis (archaea)
# 1603555  1603555  Candidatus Nitrosotenuis cloacae (archaea)
# 1580092  1580092  Candidatus Nitrosopumilus adriaticus (archaea)

# Archaeal-eukaryote tRNA analysis contains these additional species
# eukaryotes:
# Dictyostelium discoideum AX4 (cellular slime molds) dicDis1 352472
# Heterostelium album PN500 (cellular slime molds)  hetAlb1 670386
# Dictyostelium purpureum (cellular slime molds)  dicPur1 5786
# Cavenderia fasciculata (cellular slime molds) cavFas1 261658
# Entamoeba histolytica HM-1:IMSS (eukaryotes)  entHis1 294381
# Entamoeba dispar SAW760 (eukaryotes)  entDis1 370354
# Entamoeba nuttalli P19 (eukaryotes) entNut1 1076696
# Acanthamoeba castellanii str. Neff (eukaryotes) acaCas1 1257118
# Entamoeba invadens IP1 (eukaryotes) entInv1 370355
# Acytostelium subglobosum LB1 (cellular slime molds) acySub1 1410327
# lokiarchaea:
# Lokiarchaeum sp. GC14_75 (archaea)  lokiSp_GC14_75  1538547
# Candidatus Lokiarchaeota archaeon CR_4 (archaea)  candLoki_CR_4 1849166
# Candidatus Lokiarchaeota archaeon (archaea) candLoki_B53_G9 2053489

approved_taxids = c('559292', '1001994', '859192', '859350', '1229909', '1027373', '1027374', '1353246', '926571', '1407055', '1410606', '1603555', '1580092', '1776294', '1898749', '1846278', '2045011', '1410606', '718286', '1465461', '1499975', '1903277')
approved_taxids = c(approved_taxids, '352472', '670386', '5786', '261658', '294381', '370354', '1076696', '1257118', '370355', '1410327', '1538547', '1849166', '2053489')

ncbi_refs = read.delim('/projects/lowelab/db/ncbi/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt', stringsAsFactors = FALSE, sep = '\t', comment.char = '#', header = FALSE) %>% 
  filter(V5 %in% c('representative genome', 'reference genome'))
genome_table = genome_table %>% filter(taxid %in% ncbi_refs$V6 | taxid %in% ncbi_refs$V7 | taxid %in% approved_taxids)

write.table(genome_table, file = 'genomes.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

# Re-filter full taxonomy table based on presence in genomes table
filtered_taxids = genome_table %>% select(-dbname, -taxid, -name) %>% 
  gather(rank, taxid) %>%
  filter(taxid != '') %>%
  select(taxid) %>% unique %>% unlist %>% unname
taxonomy = taxonomy %>% filter(taxid %in% filtered_taxids)

write.table(taxonomy, file = 'taxonomy.tsv', quote = FALSE, sep = '\t', row.names = FALSE)
