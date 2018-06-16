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
    'stemGC': float,
    'insertions': int,
    'deletions': int,
    'intron_length': int,
    'dloop': int,
    'acloop': int,
    'tpcloop': int,
    'varm': int
  })
  message('done\n')

  trnas = trnas[trnas.primary]


  message('Parsing taxonomic ranks...')
  classifs = OrderedDict()
  for taxclass in ['genus', 'family', 'order', 'subclass', 'class', 'subphylum', 'phylum', 'subkingdom', 'kingdom', 'domain']:
    for classif in trnas[taxclass].unique():
      classifs[classif] = taxclass
  message('done\n')


  message('Calculating consensus features for {} taxonomic classifications...\n'.format(len(classifs)))
  consensus = pd.DataFrame()
  for classif in classifs:
    message('\tResolving consensus for {} {}...'.format(classifs[classif], classif))
    current_taxclass_trnas = trnas[trnas[classifs[classif]] == classif]
    current_taxclass_consensus = resolve_consensus_isotypes(current_taxclass_trnas)
    current_taxclass_consensus['classif'] = classif
    current_taxclass_consensus['rank'] = classifs[classif]
    consensus = consensus.append(current_taxclass_consensus)
    message('done\n')
  message('Done\n')
  
  message('Exporting results to {}...'.format(output_file))
  consensus.to_csv(path_or_buf = output_file, sep = '\t', index = False)
  message('done\n')


def resolve_consensus_isotypes(trnas):
  '''Resolve consensus features across all isotypes given a set of tRNAs'''
  consensus = resolve_consensus(trnas)
  consensus['isotype'] = 'All'
  for isotype in isotypes:
    current_isotype_consensus = resolve_consensus(trnas.loc[trnas.isotype == isotype])
    current_isotype_consensus['isotype'] = isotype
    consensus = consensus.append(current_isotype_consensus)
  return(consensus)

def resolve_consensus(trnas):
  '''Resolve consensus features across positions given a set of tRNAs'''
  consensus = [] # to be converted into a pandas dataframe later

  for position in positions:
    current_position = {'position': positions[position]}
    freqs = trnas.loc[:, positions[position]].value_counts(normalize = True)
    candidate_features = get_candidate_features(freqs.keys(), combos)
    for candidate in candidate_features:
      freq_check = freqs[freqs.index.isin(combos[candidate]) & (freqs >= 0.05)].sum() > 0.9
      species_check = all(
        trnas.loc[:, [positions[position], 'assembly']].groupby(
          'assembly', group_keys = False
        ).apply(
          lambda subset, position, combo: any(subset.loc[:, position].isin(combo)), 
          position = positions[position],
          combo = combos[candidate]
        )
      )
      if species_check and freq_check:
        current_position['consensus'] = candidate
        break
    if 'consensus' not in current_position:
      if ':' in positions[position]:
        current_position['consensus'] = 'N:N'
      else:
        current_position['consensus'] = 'N'
    consensus.append(current_position)
  return pd.DataFrame(consensus)

def get_candidate_features(features, combos):
  candidates = []
  for combo in combos: #  e.g., ('A', 'G')
    # current combo is a candidate feature if each letter in the combo exists in the feature set
    if numpy.all(numpy.isin(combos[combo], features)):
      candidates.append(combo)
  return candidates


def parse_args():
  parser = argparse.ArgumentParser(description = "Generate table of tRNA features using tRNAscan-SE output")
  parser.add_argument('-t', '--input_tRNAs', default = 'tRNAs.tsv', help = '')
  parser.add_argument('-o', '--output_file', default = 'consensus-{}.tsv'.format(timestamp), help = '')
  return parser.parse_args()


if __name__ == '__main__':
  timestamp = '{:%m%d%y-%H%M%S}'.format(datetime.datetime.now())
  args = parse_args()
  input_tRNAs = args.input_tRNAs
  output_file = args.output_file

  isotypes = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'iMet', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
  positions = [('p1_72', '1:72'), ('p1', '1'), ('p2_71', '2:71'), ('p2', '2'), ('p3_70', '3:70'), ('p3', '3'), ('p4_69', '4:69'), ('p4', '4'), ('p5_68', '5:68'), ('p5', '5'), ('p6_67', '6:67'), ('p6', '6'), ('p7_66', '7:66'), ('p7', '7'), ('p8', '8'), ('p8_14', '8:14'), ('p9', '9'), ('p9_23', '9:23'), ('p10_25', '10:25'), ('p10', '10'), ('p10_45', '10:45'), ('p11_24', '11:24'), ('p11', '11'), ('p12_23', '12:23'), ('p12', '12'), ('p13_22', '13:22'), ('p13', '13'), ('p14', '14'), ('p15', '15'), ('p15_48', '15:48'), ('p16', '16'), ('p17', '17'), ('p17a', '17a'), ('p18', '18'), ('p18_55', '18:55'), ('p19', '19'), ('p19_56', '19:56'), ('p20', '20'), ('p20a', '20a'), ('p20b', '20b'), ('p21', '21'), ('p22', '22'), ('p22_46', '22:46'), ('p23', '23'), ('p24', '24'), ('p25', '25'), ('p26', '26'), ('p26_44', '26:44'), ('p27_43', '27:43'), ('p27', '27'), ('p28_42', '28:42'), ('p28', '28'), ('p29_41', '29:41'), ('p29', '29'), ('p30_40', '30:40'), ('p30', '30'), ('p31_39', '31:39'), ('p31', '31'), ('p32', '32'), ('p33', '33'), ('p34', '34'), ('p35', '35'), ('p36', '36'), ('p37', '37'), ('p38', '38'), ('p39', '39'), ('p40', '40'), ('p41', '41'), ('p42', '42'), ('p43', '43'), ('p44', '44'), ('p45', '45'), ('pV11_V21', 'V11:V21'), ('pV12_V22', 'V12:V22'), ('pV13_V23', 'V13:V23'), ('pV14_V24', 'V14:V24'), ('pV15_V25', 'V15:V25'), ('pV16_V26', 'V16:V26'), ('pV17_V27', 'V17:V27'), ('pV1', 'V1'), ('pV2', 'V2'), ('pV3', 'V3'), ('pV4', 'V4'), ('pV5', 'V5'), ('pV11', 'V11'), ('pV12', 'V12'), ('pV13', 'V13'), ('pV14', 'V14'), ('pV15', 'V15'), ('pV16', 'V16'), ('pV17', 'V17'), ('pV21', 'V21'), ('pV22', 'V22'), ('pV23', 'V23'), ('pV24', 'V24'), ('pV25', 'V25'), ('pV26', 'V26'), ('pV27', 'V27'), ('p46', '46'), ('p47', '47'), ('p48', '48'), ('p49_65', '49:65'), ('p49', '49'), ('p50_64', '50:64'), ('p50', '50'), ('p51_63', '51:63'), ('p51', '51'), ('p52_62', '52:62'), ('p52', '52'), ('p53_61', '53:61'), ('p53', '53'), ('p54', '54'), ('p54_58', '54:58'), ('p55', '55'), ('p56', '56'), ('p57', '57'), ('p58', '58'), ('p59', '59'), ('p60', '60'), ('p61', '61'), ('p62', '62'), ('p63', '63'), ('p64', '64'), ('p65', '65'), ('p66', '66'), ('p67', '67'), ('p68', '68'), ('p69', '69'), ('p70', '70'), ('p71', '71'), ('p72', '72'), ('p73', '73')]
  combos = [('A', ('A',)), ('C', ('C',)), ('G', ('G',)), ('U', ('U',)), ('Gap', ("-", ".", "-:-")), 
    ('GC', ("G:C",)), ('AU', ("A:U",)), ('UA', ("U:A",)), ('CG', ("C:G",)), ('GU', ("G:U",)), ('UG', ("U:G",)), 
    ('Purine', ("A", "G")), ('Pyrimidine', ("C", "U")), ('StrongPair', ("G:C", "C:G")), 
    ('WeakPair', ("A:U", "U:A")), ('WobblePair', ("G:U", "U:G")),
    ('Weak', ("A", "U")), ('Strong', ("G", "C")), 
    ('Amino', ("A", "C")), ('Keto', ("G", "U")), 
    ('PurinePyrimidine', ("A:U", "G:C")), ('PyrimidinePurine', ("U:A", "C:G")),
    ('AminoKeto', ('A:U', 'C:G')), ('KetoAmino', ('U:A', 'G:C')), 
    ('B', ("C", "G", "U")), ('D', ("A", "G", "U")), ('H', ("A", "C", "U")), ('V', ("A", "C", "G")), 
    ('Paired', ("A:U", "U:A", "C:G", "G:C", "G:U", "U:G")), 
    ('Mismatched', ("A:A", "G:G", "C:C", "U:U", "A:G", "A:C", "C:A", "C:U", "G:A", "U:C")), 
    ('Bulge', ("A:-", "U:-", "C:-", "G:-", "-:A", "-:G", "-:C", "-:U"))]
  positions = OrderedDict(positions)
  combos = OrderedDict(combos)
  main()