#!/usr/bin/env python
import sys, os
message = lambda s: print(s, file = sys.stderr, end = '', flush = True)
message('Loading packages...')
import pandas as pd
import argparse
import subprocess
from Bio import SeqIO, SeqRecord, Seq
message('done\n')

def main():
  args = parse_args()
  message('Loading genome data into data frame...')
  genomes_df = read_genomes_table(args.genome_table_path)
  message('done\n')


def parse_args():
  parser = argparse.ArgumentParser(description = "Generate table of tRNA features using tRNAscan-SE output")
  parser.add_argument('-g', '--genome_table_path', default = 'genomes.tsv', help = '')
  return(parser.parse_args())

def read_genomes_table(path):
  '''Read genomes table, then add and validate paths to tRNAscan-SE output to genomes data frame'''
  genomes_df = pd.read_table(path, header = None, index_col = False,
    names = ['dbname', 'dirname', 'taxid', 'assembly', 'varietas', 'species', 'genus', 'subclass', 'class', 'suborder', 'order', 'subphylum', 'phylum', 'subkingdom', 'kingdom'],
    dtype = {'taxid': str})

  genomes_df['iso'] = genomes_df.dbname.apply(lambda dbname: 'data/iso/{}-tRNAs.iso'.format(dbname))
  genomes_df['ss'] = genomes_df.dbname.apply(lambda dbname: 'data/ss/{}-tRNAs.ss'.format(dbname))
  genomes_df['out'] = genomes_df.dbname.apply(lambda dbname: 'data/out/{}-tRNAs.out'.format(dbname))
  genomes_df['tRNAs'] = genomes_df.dbname.apply(lambda dbname: 'data/tRNAs/{}-tRNAs.fa'.format(dbname))

  # Validate files, make sure they're there
  for row in genomes_df.itertuples():
    for file_type in ['iso', 'ss', 'out']:
      if not os.path.exists(getattr(row, file_type)):
        raise Exception('.{} file {} does not exist'.format(file_type, getattr(row, file_type)))

  return(genomes_df)


def process_tscan_output(genomes_df):
  seqs = read_tRNAs(genomes_df)
  align_seqs(seqs)
  trnas = parse_alignment_into_trnas()
  trnas = annotate_trnas(trnas)
  save(trnas)

def align_seqs(seqs):
  ''' '''
  pass

def parse_alignment_into_trnas():
  pass


if __name__ == '__main__':
  main()