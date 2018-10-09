#!/usr/bin/env python
import sys, os
message = lambda s: print(s, file = sys.stderr, end = '', flush = True)
message('Loading packages...')
import pandas as pd
import numpy
import argparse
import datetime
from collections import OrderedDict
message('done\n')

def main():
  # Read in tRNAs

  message('Loading tRNAs from {}...'.format(input_tRNAs))
  trnas = pd.read_csv(input_tRNAs, sep = '\t', dtype = {
    'seqname': str,
    'isotype': str,
    'anticodon': str,
    'score': float,
    'primary': bool,
    'best_model': str,
    'isoscore': float,
    'isoscore_ac': float,
    'dbname': str,
    'assembly': str,
    'varietas': str,
    'species': str,
    'genus': str,
    'family': str,
    'order': str,
    'subclass': str,
    'class': str,
    'subphylum': str,
    'phylum': str,
    'subkingdom': str,
    'kingdom': str,
    'domain': str,
    'taxid': str,
    'GCcontent': float,
    'insertions': int,
    'deletions': int,
    'intron_length': int,
    'dloop': int,
    'acloop': int,
    'tpcloop': int,
    'varm': int
  })
  trnas = trnas[trnas.primary].fillna('')
  message('done\n')

  message('Parsing clades and taxonomic ranks...')
  taxonomy = pd.read_csv(taxonomy_file, sep = '\t', dtype = {'name': str, 'rank': str, 'taxid': str})
  message('done\n')

  message('Counting freqs for {} clades...\n'.format(taxonomy[taxonomy['rank'] != 'assembly'].shape[0]))
  if saved and os.path.exists('freqs-saved.pkl'):
    message('\tLoading saved data...')
    freqs = pd.read_pickle('freqs-saved.pkl')
    if (freqs.size == 0): 
      message('found empty freqs file, did not load data\n')
      freqs = pd.DataFrame()
    else: message('done\n')
  else:
    freqs = pd.DataFrame()
  for name, rank, taxid in taxonomy.itertuples(index = False):
    if rank == 'assembly': continue
    try:
      message('\tCounting base frequencies for {} {}...'.format(rank, name))
      if saved and 'taxid' in freqs.columns and 'rank' in freqs.columns and not freqs[(freqs['rank'] == rank) & (freqs['taxid'] == taxid)].empty:
        message('using saved data from previous run\n')
        continue
      current_clade_trnas = trnas[trnas[rank] == name]
      current_clade_freqs = count_freqs_isotypes(current_clade_trnas)
      current_clade_freqs['taxid'] = taxid
      current_clade_freqs['rank'] = rank
      freqs = freqs.append(current_clade_freqs, sort = True)
      message('done\n')
    except Exception as e:
      freqs.to_pickle('freqs-saved.pkl')
      raise e
    # Save every 100 records too. 132 positions * 24 isotypes = 3168 rows
    if saved and freqs.shape[0] % 316800 == 0:
      freqs.to_pickle('freqs-saved.pkl')
      message('- Saving checkpoint - \n')
  message('Done\n')
  
  message('Exporting results to {}...'.format(output_file))
  colorder = ['taxid', 'rank', 'isotype', 'position', 'total'] + features
  freqs = freqs.reindex(columns = colorder)
  freqs.to_csv(path_or_buf = output_file, sep = '\t', index = False)
  message('done\n')


def count_freqs_isotypes(trnas):
  '''Count freqs all isotypes given a set of tRNAs'''
  freqs = count_freqs(trnas)
  freqs['isotype'] = 'All'
  for isotype in isotypes:
    current_isotype_freqs = count_freqs(trnas.loc[trnas.isotype == isotype])
    current_isotype_freqs['isotype'] = isotype
    freqs = freqs.append(current_isotype_freqs, sort = True)
  return(freqs)

def count_freqs(trnas):
  '''Count freqas across positions given a set of tRNAs'''
  freqs = [] # to be converted into a pandas dataframe later
  for position in positions.values():
    current_position = trnas.loc[trnas[position].isin(features), position].value_counts()
    current_position['total'] = sum(current_position)
    current_position['position'] = position
    freqs.append(current_position)

  return pd.DataFrame(freqs, columns = features + ['total']).fillna(0).astype(int).reset_index().rename(columns = {'index': 'position'})

def parse_args():
  parser = argparse.ArgumentParser(description = "Generate table of tRNA features using tRNAscan-SE output")
  parser.add_argument('-t', '--input_tRNAs', default = 'tRNAs.tsv', help = '')
  parser.add_argument('-c', '--taxonomy', default = 'taxonomy.tsv', help = '')
  parser.add_argument('-o', '--output_file', default = 'freqs-{}.tsv'.format(timestamp), help = '')
  parser.add_argument('-s', '--saved', default = False, action = 'store_true', help = 'Use and save data before/after a failed run')
  return parser.parse_args()


if __name__ == '__main__':
  timestamp = '{:%m%d%y-%H%M%S}'.format(datetime.datetime.now())
  args = parse_args()
  input_tRNAs = args.input_tRNAs
  output_file = args.output_file
  taxonomy_file = args.taxonomy
  saved = args.saved

  isotypes = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Ile2', 'fMet', 'iMet', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
  positions = [('p1_72', '1:72'), ('p1', '1'), ('p2_71', '2:71'), ('p2', '2'), ('p3_70', '3:70'), ('p3', '3'), ('p4_69', '4:69'), ('p4', '4'), ('p5_68', '5:68'), ('p5', '5'), ('p6_67', '6:67'), ('p6', '6'), ('p7_66', '7:66'), ('p7', '7'), ('p8', '8'), ('p8_14', '8:14'), ('p9', '9'), ('p9_23', '9:23'), ('p10_25', '10:25'), ('p10', '10'), ('p10_45', '10:45'), ('p11_24', '11:24'), ('p11', '11'), ('p12_23', '12:23'), ('p12', '12'), ('p13_22', '13:22'), ('p13', '13'), ('p14', '14'), ('p15', '15'), ('p15_48', '15:48'), ('p16', '16'), ('p17', '17'), ('p17a', '17a'), ('p18', '18'), ('p18_55', '18:55'), ('p19', '19'), ('p19_56', '19:56'), ('p20', '20'), ('p20a', '20a'), ('p20b', '20b'), ('p21', '21'), ('p22', '22'), ('p22_46', '22:46'), ('p23', '23'), ('p24', '24'), ('p25', '25'), ('p26', '26'), ('p26_44', '26:44'), ('p27_43', '27:43'), ('p27', '27'), ('p28_42', '28:42'), ('p28', '28'), ('p29_41', '29:41'), ('p29', '29'), ('p30_40', '30:40'), ('p30', '30'), ('p31_39', '31:39'), ('p31', '31'), ('p32', '32'), ('p33', '33'), ('p34', '34'), ('p35', '35'), ('p36', '36'), ('p37', '37'), ('p38', '38'), ('p39', '39'), ('p40', '40'), ('p41', '41'), ('p42', '42'), ('p43', '43'), ('p44', '44'), ('p45', '45'), ('pV11_V21', 'V11:V21'), ('pV12_V22', 'V12:V22'), ('pV13_V23', 'V13:V23'), ('pV14_V24', 'V14:V24'), ('pV15_V25', 'V15:V25'), ('pV16_V26', 'V16:V26'), ('pV17_V27', 'V17:V27'), ('pV1', 'V1'), ('pV2', 'V2'), ('pV3', 'V3'), ('pV4', 'V4'), ('pV5', 'V5'), ('pV11', 'V11'), ('pV12', 'V12'), ('pV13', 'V13'), ('pV14', 'V14'), ('pV15', 'V15'), ('pV16', 'V16'), ('pV17', 'V17'), ('pV21', 'V21'), ('pV22', 'V22'), ('pV23', 'V23'), ('pV24', 'V24'), ('pV25', 'V25'), ('pV26', 'V26'), ('pV27', 'V27'), ('p46', '46'), ('p47', '47'), ('p48', '48'), ('p49_65', '49:65'), ('p49', '49'), ('p50_64', '50:64'), ('p50', '50'), ('p51_63', '51:63'), ('p51', '51'), ('p52_62', '52:62'), ('p52', '52'), ('p53_61', '53:61'), ('p53', '53'), ('p54', '54'), ('p54_58', '54:58'), ('p55', '55'), ('p56', '56'), ('p57', '57'), ('p58', '58'), ('p59', '59'), ('p60', '60'), ('p61', '61'), ('p62', '62'), ('p63', '63'), ('p64', '64'), ('p65', '65'), ('p66', '66'), ('p67', '67'), ('p68', '68'), ('p69', '69'), ('p70', '70'), ('p71', '71'), ('p72', '72'), ('p73', '73')]
  features = ['A', 'C', 'G', 'U', '-', 'A:U', 'U:A', 'G:C', 'C:G', 'G:U', 'U:G', '-:-', 'A:A', 'A:C', 'A:G', 'C:A', 'C:C', 'C:U', 'G:A', 'G:G', 'U:C', 'U:U', 'A:-', '-:A', 'C:-', '-:C', 'G:-', '-:G', 'U:-', '-:U']
  positions = OrderedDict(positions)
  main()