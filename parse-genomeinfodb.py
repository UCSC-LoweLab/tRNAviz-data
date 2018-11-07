import re

taxids = []
taxid = ''
for line in open('Genome-info-arch-new'):
  if line == '\n' and taxid not in taxids and taxid != '':
    print('{}\t{}\t{}'.format(dirname, dbname, taxid))
    taxids.append(taxid)
  elif line[0:7] == 'Org_abr':
    dirname = line.strip().split()[-1]
  elif line[0:7] == 'DB_name':
    dbname = line.strip().split()[-1]
  elif line[0:6] == 'Tax_ID':
    taxid = line.split()[1][:-1]
