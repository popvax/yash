"""PopRNA synonymous-encoding heuristic for the fixed vaccine protein.

Adapted from the user-supplied PopRNAOptimizer_final.py. Parameters and original
source SHA-256 are recorded in backbone.json; the builder verifies translation,
rare-codon exclusion and motif constraints independently.
"""
import python_codon_tables as pct
from Bio.Seq import Seq
import argparse
import random
import re
import math
import json
from functools import lru_cache
import RNA as vienna
#from numba import jit

# creat help message
def parse_arguments():
  parser = argparse.ArgumentParser(
    description="This program takes a protein sequence and find synthesizable, efficiend and high yeild RNA for mRNA therapy",
    formatter_class=argparse.RawTextHelpFormatter
  )

# Parse command-line arguments
  parser = argparse.ArgumentParser(description="Find optimal sequence for a given protein sequence")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument("-p", "--protein", action='store_true', help="Input is protein sequence.")
  group.add_argument("-d", "--dna", action='store_true', help="Input is DNA sequence.")
  group.add_argument("-r", "--rna", action='store_true', help="Input is RNA sequence.")
  parser.add_argument("-i", "--in_put", required=True, help="Input sequence")
  parser.add_argument("-o", "--output", required=False, help="Output text file")
  parser.add_argument("-c", "--cutoff", required=False, default=0.1, type=float, help="Cut off value for codons")
  parser.add_argument("-or", "--organism", required=False, default='h_sapiens', help="Available organism name : b_subtilis_1423, c_elegans_6239, d_melanogaster_7227, e_coli_316407, g_gallus_9031, h_sapiens_9606, m_musculus_10090, m_musculus_domesticus_10092, s_cerevisiae_4932")
  parser.add_argument("-s", "--sequences", required=False, default=['UUUU', 'UUUC', 'UCUAGA', 'GGUCUC', 'GAGACC', 'AAAAAAAA', 'GGGGGGGG', 'CCCCCCCC'], nargs='+', help="Sequence to avoid")
  parser.add_argument("--seed", type=int, default=0)
  parser.add_argument("--fold-iterations", type=int, default=100)
  args = parser.parse_args()

# Display help message
  if hasattr(args, "help"):
    parser.print_help()
    exit()

  return args

args = parse_arguments()
random.seed(args.seed)
vienna.cvar.temperature = 37.0
table = pct.get_codons_table(args.organism, replace_U_by_T=False)
AA2codon = {'A': ('GCA', 'GCC', 'GCG', 'GCU'), 'C': ('UGC', 'UGU'), 'D': ('GAC', 'GAU'), 'E': ('GAA', 'GAG'), 'F': ('UUC', 'UUU'), 'G': ('GGA', 'GGC', 'GGG', 'GGU'), 'H': ('CAC', 'CAU'), 'I': ('AUA', 'AUC', 'AUU'), 'K': ('AAA', 'AAG'), 'L': ('CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'), 'M': ['AUG'], 'N': ('AAC', 'AAU'), 'P': ('CCA', 'CCC', 'CCG', 'CCU'), 'Q': ('CAA', 'CAG'), 'R': ('AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'), 'S': ('AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'), 'T': ('ACA', 'ACC', 'ACG', 'ACU'), 'V': ('GUA', 'GUC', 'GUG', 'GUU'), 'W': ['UGG'], 'Y': ('UAC', 'UAU')}
AA_names = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
codon2AA = {'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 'UGC':'C', 'UGU': 'C', 'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 'UUC':'F', 'UUU':'F', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G', 'CAC':'H', 'CAU':'H', 'AUA':'I', 'AUC':'I', 'AUU':'I', 'AAA':'K', 'AAG':'K', 'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 'UUA':'L', 'UUG':'L', 'AUG':'M', 'AAC':'N', 'AAU':'N', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P', 'CAA':'Q', 'CAG':'Q', 'AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R', 'AGC':'S', 'AGU':'S', 'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T', 'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 'UGG':'W', 'UAC':'Y', 'UAU':'Y'}
excluded_codon = [x for x in codon2AA.keys() if table[codon2AA[x]][x]<args.cutoff] # collection of rare codons
AA2_freq_codon = {key: tuple(x for x in values if x not in excluded_codon) if isinstance(values, tuple) else values for key, values in AA2codon.items()}
AA2_freq_codon['M'] = ['AUG']
AA2_freq_codon['W'] = ['UGG']

def modify_table(table, cutoff=0.1): #modify codon table to remove rare codons verified
  for aa in table.keys():
    adjust = 0
    for cod in table[aa].keys():
      if table[aa][cod] < cutoff:
        adjust += table[aa][cod]
        table[aa][cod] = 0
    if adjust:
      adjust_fact = 1-adjust
      for cod1 in table[aa].keys():
        table[aa][cod1] = round(table[aa][cod1]/adjust_fact, 2)

  return table

modi_table = modify_table(table, args.cutoff)

#@jit
def modify_RNA(RNA, start, sequ):
  return RNA[:start] + sequ + RNA[start+len(sequ):]

#@jit
def check_codon_replacement(RNA, protein, aa_num, new_codon, seq2avoid):
  """Reject forbidden motifs across either side of the changed codon."""
  position = 3 * aa_num
  candidate = modify_RNA(RNA, position, new_codon)
  padding = max(map(len, seq2avoid), default=1) - 1
  window = candidate[max(0, position-padding):position+3+padding]
  return not any(motif in window for motif in seq2avoid)

def replace_seq(RNA, n, protein):
  """Bounded initial repair; do not silently leave forbidden motifs behind."""
  def count(sequence):
    return sum(sum(sequence.startswith(motif, i) for i in range(len(sequence)))
               for motif in seq2avoid)
  for iteration in range(n, 101):
    total = count(RNA)
    if total == 0:
      return RNA
    changed = False
    for motif in seq2avoid:
      for pos in range(len(RNA)-len(motif)+1):
        if not RNA.startswith(motif, pos):
          continue
        positions = list(range(pos//3, (pos+len(motif)-1)//3+1))
        random.shuffle(positions)
        for aa_num in positions:
          choices = [codon for codon in AA2_freq_codon[protein[aa_num]]
                     if codon != RNA[3*aa_num:3*aa_num+3]]
          random.shuffle(choices)
          for codon in choices:
            candidate = modify_RNA(RNA, aa_num*3, codon)
            if count(candidate) < total:
              RNA = candidate
              changed = True
              break
          if changed:
            break
        if changed:
          break
      if changed:
        break
    if not changed:
      raise RuntimeError('Cannot repair forbidden motifs with a single allowed synonymous change')
  raise RuntimeError('Initial motif repair exceeded 100 edits')

#@jit
def gc_count(codon):
    return (codon.count('G') + codon.count('C'))

def find_CAI_codon(table): #find most frequent code for each amino acids verified
  aa2CAI={}
  for aa in AA_names:
    aa2CAI[aa] = max(table[aa], key=table[aa].get)

  return aa2CAI

def optimise_CAI(seq, aa2codon): #initialise RNA sequence, not required if the input sequence is DNA or RNA verified
  RNA_seq = ''
  for aa in seq:
    RNA_seq += aa2codon[aa]

  return RNA_seq

#@jit
def opt_first_five_codon(RNA, protein, table):# simultaneously optimise the codon usage and GC content for first 5 aminoacid after ATG
  for i in range(1, 6):
  #codon preference = gc content * use frequency, choose the best codon
    codons = []
    for cod in table[protein[i]]:
      if cod in excluded_codon:
        continue
      codons.append([cod, table[protein[i]][cod]*gc_count(cod)])
    codons = sorted(codons, key=lambda x: x[1], reverse=True)
    for y in codons:
      if check_codon_replacement(RNA, protein, i, y[0], seq2avoid):
        RNA = RNA[:i*3] + y[0] + RNA[i*3+3:]
        break
  return RNA

#@jit
def find_patern_pose(protein, first, second):
  positions = []

  for i in range(len(protein) - 1):
    current_char = protein[i]
    next_char = protein[i + 1]

    if (current_char in first) and (next_char in second):
      positions.append(i)

  return positions

#@jit
def add_stop_one(RNA, protein, fixed=set(), prefered=set()):
  first = ('L', 'V')
  second = ('G', 'D', 'E', 'V', 'A')
  L = ('CUG', 'UUG')
  V = ['GUG']
  patern = find_patern_pose(protein, first, second)
  for entry in patern:
    if all(entry*3 + i not in fixed for i in range(6)):
      match protein[entry]:
        case 'L':
          match protein[entry+1]:
            case 'G':
              for code in L:
                for cod in AA2_freq_codon['G']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break
            case 'D':
              for code in L:
                for cod in AA2_freq_codon['D']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break

            case 'E':
              for code in L:
                for cod in AA2_freq_codon['E']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break

            case 'V':
              for code in L:
                for cod in AA2_freq_codon['V']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break

            case 'A':
              for code in L:
                for cod in AA2_freq_codon['A']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break
        case 'V':
          match protein[entry+1]:
            case 'G':
              for code in V:
                for cod in AA2_freq_codon['G']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break

            case 'D':
              for code in V:
                for cod in AA2_freq_codon['D']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break

            case 'E':
              for code in V:
                for cod in AA2_freq_codon['E']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break

            case 'V':
              for code in V:
                for cod in AA2_freq_codon['V']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break

            case 'A':
              for code in V:
                for cod in AA2_freq_codon['A']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break

  first = ('I', 'V')
  second = ('G', 'D', 'E', 'V', 'A', 'N', 'K', 'I', 'M', 'T')
  patern = find_patern_pose(protein, first, second)
  for entry in patern:
    if all(entry*3 + i not in fixed for i in range(6)):
      match protein[entry]:
        case 'I':
          match protein[entry+1]:
            case 'G':
              for cod in AA2_freq_codon['G']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'V':
              for cod in AA2_freq_codon['V']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'A':
              for cod in AA2_freq_codon['A']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'I':
              for cod in AA2_freq_codon['I']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'M':
              for cod in AA2_freq_codon['M']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'T':
              for cod in AA2_freq_codon['T']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'V':
          match protein[entry+1]:
            case 'G':
              for cod in AA2_freq_codon['G']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'V':
              for cod in AA2_freq_codon['V']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'A':
              for cod in AA2_freq_codon['A']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'I':
              for cod in AA2_freq_codon['I']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'M':
              for cod in AA2_freq_codon['M']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'T':
              for cod in AA2_freq_codon['T']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
  second = ('R', 'S')
  R=['AGA', 'AGG']
  S=['AGC', 'AGU']
  patern = find_patern_pose(protein, first, second)
  for entry in patern:
    if all(entry*3 + i not in fixed for i in range(6)):
      match protein[entry]:
        case 'I':
          match protein[entry+1]:
            case 'R':
              random.shuffle(R)
              for cod in R:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  prefered = prefered.union(set([entry*3, entry*3+1, entry*3+2, entry*3+3, entry*3+4, entry*3+5]))
                  break
            case 'S':
              random.shuffle(S)
              for cod in S:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUA' + cod + RNA[entry*3+6:]
                  prefered = prefered.union(set([entry*3, entry*3+1, entry*3+2, entry*3+3, entry*3+4, entry*3+5]))
                  break
        case 'V':
          match protein[entry+1]:
            case 'R':
              random.shuffle(R)
              for cod in R:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  prefered = prefered.union(set([entry*3, entry*3+1, entry*3+2, entry*3+3, entry*3+4, entry*3+5]))
                  break
            case 'S':
              random.shuffle(S)
              for cod in S:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GUA'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GUA' + cod + RNA[entry*3+6:]
                  prefered = prefered.union(set([entry*3, entry*3+1, entry*3+2, entry*3+3, entry*3+4, entry*3+5]))
                  break

  return RNA, fixed, prefered

def add_stop_two(RNA, protein, fixed=set(), prefered=set()):
  first = ('N', 'S', 'T', 'Y', 'C', 'D', 'G', 'A', 'H', 'P')
  second = ('D', 'E', 'N', 'K')
  S = ['AGU', 'UCU']
  patern = find_patern_pose(protein, first, second)
  for entry in patern:
    if all(entry*3 + i not in fixed for i in range(6)):
      match protein[entry]:
        case 'N':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'AUU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'AUU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'T':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'ACU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'ACU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'ACU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'ACU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'ACU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'ACU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'ACU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'ACU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'Y':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'UAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'UAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'UAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'UAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'UAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'UAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'UAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'UAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'C':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'UGU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'UGU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'UGU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'UGU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'UGU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'UGU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'UGU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'UGU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'S':
          random.shuffle(S)
          match protein[entry+1]:
            case 'D':
              for code in S:
                for cod in AA2_freq_codon['D']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break
            case 'E':
              for code in S:
                for cod in AA2_freq_codon['E']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break
            case 'N':
              for code in S:
                for cod in AA2_freq_codon['N']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break
            case 'K':
              for code in S:
                for cod in AA2_freq_codon['K']:
                  if check_codon_replacement(modify_RNA(RNA, entry*3, code), protein, entry+1, cod, seq2avoid):
                    RNA = RNA[:entry*3] + code + cod + RNA[entry*3+6:]
                    fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                    break
                break
        case 'D':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'G':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GGU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GGU' + cod + RNA[entry*3+6:]
                  prefered = prefered.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GGU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GGU' + cod + RNA[entry*3+6:]
                  prefered = prefered.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GGU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GGU' + cod + RNA[entry*3+6:]
                  prefered = prefered.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GGU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GGU' + cod + RNA[entry*3+6:]
                  prefered = prefered.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'A':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GCU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GCU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GCU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GCU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GCU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GCU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'GCU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'GCU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'H':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'CAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'CAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'CAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'CAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'CAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'CAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'CAU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'CAU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
        case 'P':
          match protein[entry+1]:
            case 'D':
              for cod in AA2_freq_codon['D']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'CCU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'CCU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'E':
              for cod in AA2_freq_codon['E']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'CCU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'CCU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'N':
              for cod in AA2_freq_codon['N']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'CCU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'CCU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
            case 'K':
              for cod in AA2_freq_codon['K']:
                if check_codon_replacement(modify_RNA(RNA, entry*3, 'CCU'), protein, entry+1, cod, seq2avoid):
                  RNA = RNA[:entry*3] + 'CCU' + cod + RNA[entry*3+6:]
                  fixed = fixed.union(set([entry*3, entry*3+1, entry*3+2]))
                  break
  return RNA, fixed, prefered

#@jit
def cal_freq(codon, RNA): #calculate frequency of each codon present in given sequence verified
  if codon in ('AUG', 'UGG'):
    return 1
  else:
    codon_lis = [RNA[i:i+3] for i in range(0, len(RNA), 3)]
    codon_count = codon_lis.count(codon)
    aa = codon2AA[codon]
    syn_codon = AA2codon[aa]
    total_count = 0
    for cod in syn_codon:
      total_count += codon_lis.count(cod)
    if total_count:
      return codon_count/total_count
    else:
      return 0

@lru_cache(maxsize=1024)
def rna_sec_structure(sequence):
  """Same default MFE folding model as RNAfold; use its pinned Python API."""
  structure, energy = vienna.fold(sequence)
  if len(structure) != len(sequence):
    raise RuntimeError('ViennaRNA returned an unexpected structure length')
  return structure

#@jit
def find_continuous_repeats(string, min_length): # find continous repeats of ( or ) verified
  pattern = re.compile(f'[()]{{{min_length},}}')
  matches = pattern.finditer(string)

  repeats = []
  for match in matches:
    repeat_start = match.start()
    repeat_length = match.end() - repeat_start
    if repeat_length >= min_length:
      repeats.append([repeat_start, repeat_length])

  return repeats

#@jit
def replace_string_element(RNA, protein, aa_num, exclude_elements, seq2avoid):
  """Try each allowed synonym once; keep the sequence if none is feasible."""
  excluded = set(exclude_elements or ()) | set(excluded_codon)
  choices = [codon for codon in AA2codon[protein[aa_num]] if codon not in excluded]
  random.shuffle(choices)
  for codon in choices:
    if check_codon_replacement(RNA, protein, aa_num, codon, seq2avoid):
      return modify_RNA(RNA, 3*aa_num, codon)
  return RNA

#@jit
def chose_alt_cod(RNA, protein, pos, seq2avoid, exclude_elements=None):
  excluded = set(excluded_codon if exclude_elements is None else exclude_elements)
  excluded.add(RNA[pos*3:pos*3+3])
  return replace_string_element(RNA, protein, pos, excluded, seq2avoid)

#@jit
def reduce_secondary_structure(RNA, table, fixed=(), prefered=()):
  """Preserve the original paired-bracket-run heuristic, including its 14-nt target."""
  protected = set(fixed) | set(prefered)
  def longest(sequence):
    return max((length for _, length in find_continuous_repeats(rna_sec_structure(sequence), 10)), default=0)
  best = RNA
  best_run = longest(best)
  if best_run <= 14:
    print(f'folding: target already met; longest paired-bracket run={best_run}', flush=True)
    return best
  for n in range(args.fold_iterations):
    for start_nt, length_nt in find_continuous_repeats(rna_sec_structure(RNA), 10):
      start = start_nt//3 if start_nt%3 == 2 else math.ceil(start_nt/3)
      length = math.ceil(length_nt/3)
      for pos in range(start, min(start+length-1, len(protein))):
        if any(3*pos+i in protected for i in range(3)):
          continue
        RNA = chose_alt_cod(RNA, protein, pos, seq2avoid)
    current_run = longest(RNA)
    if current_run < best_run:
      best, best_run = RNA, current_run
    if n == 0 or (n+1)%10 == 0:
      print(f'folding: iteration={n+1}; current_run={current_run}; best_run={best_run}', flush=True)
    if current_run <= 14:
      break
  return best

#@jit
def chose_freq_codon(RNA, aa, table): # provide a list of codon to chose best to worst based on codon table and codon usage in RNA verified
  cod_use = {}
  if aa == 'M':
    return ['AUG']
  elif aa == 'W':
    return ['UGG']
  else:
    for cod in AA2codon[aa]:
      if cod in excluded_codon:
        continue
      cod_use.update({cod:table[aa][cod]-cal_freq(cod, RNA)})
    cod_use = dict(sorted(cod_use.items(), key=lambda item: item[1], reverse=True))

    return list(cod_use.keys())

#@jit
def optimise_codon_usage(RNA, table, fixed=(), prefered=(), threshold=1.1): # improve statistical use of codons verified
  SS = rna_sec_structure(RNA)
  paired = {pos for pos, char in enumerate(SS) if char in ('(', ')')}
  all_ss_pose, semi_ss_pose, non_ss_pose = [], [], []
  protected = set(fixed) | set(prefered)
  for start in range(0, len(RNA), 3):
    count = sum(start+i in paired for i in range(3))
    target = all_ss_pose if count == 3 else non_ss_pose if count == 0 else semi_ss_pose
    target.append(start)

  for pose in all_ss_pose:
    if any(pose+i in protected for i in range(3)):
      continue
    aa_num = math.floor(divmod(pose, 3)[0])
    pad = max(len(ele) for ele in seq2avoid) - 2
    for cod in chose_freq_codon(RNA, protein[aa_num], table):
      if cod == RNA[aa_num*3:aa_num*3+3]:
        break
      if check_codon_replacement(RNA, protein, aa_num, cod, seq2avoid):
        RNA = RNA[:aa_num*3] + cod + RNA[aa_num*3+3:]
        break

  for pose in semi_ss_pose:
    if any(pose+i in protected for i in range(3)):
      continue
    aa_num = math.floor(divmod(pose, 3)[0])
    pad = max(len(ele) for ele in seq2avoid) - 2
    for cod in chose_freq_codon(RNA, protein[aa_num], table):
      if cod == RNA[aa_num*3:aa_num*3+3]:
        break
      if check_codon_replacement(RNA, protein, aa_num, cod, seq2avoid):
        RNA = RNA[:aa_num*3] + cod + RNA[aa_num*3+3:]
        break

  for pose in non_ss_pose:
    if any(pose+i in protected for i in range(3)):
      continue
    aa_num = math.floor(divmod(pose, 3)[0])
    pad = max(len(ele) for ele in seq2avoid) - 2
    for cod in chose_freq_codon(RNA, protein[aa_num], table):
      if cod == 'AUA':
        if cal_freq(cod, RNA) < cal_freq(RNA[aa_num*3:aa_num*3+3], RNA): # if code is 'AUA' threshold ignore
          break
      elif cal_freq(cod, RNA) < cal_freq(RNA[aa_num*3:aa_num*3+3], RNA)*threshold: # if code is 'AUA' threshold ignore
        break
      if check_codon_replacement(RNA, protein, aa_num, cod, seq2avoid):
        RNA = RNA[:aa_num*3] + cod + RNA[aa_num*3+3:]
        break

  return RNA

#@jit
def codn_usage_efficiency(RNA, table): #find the codon usage efficiency verified
  total_squr_error = 0
  aa_error = {}
  for aa in AA_names:
    for cod in AA2codon[aa]:
      error = cal_freq(cod, RNA)-table[aa][cod]
      if aa in aa_error.keys():
        aa_error[aa].update({cod:error})
      else:
        aa_error[aa]= {cod:error}
      total_squr_error += error**2
  mean_error = math.sqrt(total_squr_error/20)
  return mean_error, aa_error

def remove_adjascent(in_lis, min=3, max=24):
  i = 0
  while i < len(in_lis) - 1:
    # Calculate the difference between adjacent in_lis
    difference = abs(in_lis[i + 1] - in_lis[i])
    # If the difference is less than 30 and more than 3, remove the number
    if difference < max and difference > min:
      del in_lis[i + 1]
    else:
      i += 1
  return in_lis

def find_repeats(RNA, windows=15):
  """Enumerate every exact repeated window, including the final window."""
  seen, repeats = {}, []
  for pos in range(len(RNA)-windows+1):
    word = RNA[pos:pos+windows]
    for previous in seen.get(word, ()):
      repeats.append((previous, pos))
    seen.setdefault(word, []).append(pos)
  return sorted(repeats)

def remove_repeats(RNA, table, windows, n, fixed=[]):
  repeats = find_repeats(RNA, windows)
  for rep in repeats:
    start = random.choice(rep)
    start = math.floor(start/3)
    repl_lis = [i for i in range(start, start+int(windows/3))]
    random.shuffle(repl_lis)
    for pos in repl_lis:
      if any(pos*3 + i in fixed for i in range(3)) or protein[pos] in ['M', 'W']:
        continue
      RNAi = chose_alt_cod(RNA, protein, pos, seq2avoid)
      if RNAi != RNA:
        RNA = RNAi
        break
  if find_repeats(RNA, windows) and n<=10:
    n += 1
    RNA = remove_repeats(RNA, table, windows, n, fixed)

  return RNA

def RNA2protein(sequence):
  if len(sequence)%3 or set(sequence)-set('ACGU'):
    raise ValueError('Expected coding RNA of a whole number of codons')
  return str(Seq(sequence).translate())

if __name__ == '__main__':
  seq2avoid = [x.upper().replace('T', 'U') for x in args.sequences]
  if args.protein:
    protein = args.in_put.strip().upper()
    if len(protein) < 6 or set(protein)-set(AA_names):
      raise ValueError("Provide at least six standard amino acids, without a stop")
    RNA = optimise_CAI(protein, find_CAI_codon(table))
  elif args.dna:
    RNA = args.in_put.upper().replace('T', 'U')
    protein = RNA2protein(RNA)
  else:
    RNA = args.in_put.upper()
    protein = RNA2protein(RNA)
  print('----------')
  RNA = replace_seq(RNA, 0, protein)
  RNA = opt_first_five_codon(RNA, protein, modi_table)
  fixed=set(x for x in range(18))
  print('*', end='')
  RNA_f_p = add_stop_one(RNA, protein, fixed)
  RNA = RNA_f_p[0]
  fixed = RNA_f_p[1]
  prefered = RNA_f_p[2]
  print('*', end='')
  RNA_f_p = add_stop_two(RNA, protein, fixed, prefered)
  RNA = RNA_f_p[0]
  fixed = RNA_f_p[1]
  prefered = RNA_f_p[2]
  print('*', end='')
  RNA = reduce_secondary_structure(RNA, modi_table, fixed, prefered)
  print('*', end='')
  RNA = reduce_secondary_structure(RNA, modi_table, fixed)
  print('*', end='')
  RNA = optimise_codon_usage(RNA, modi_table, fixed, prefered)
  print('*', end='')
  RNA = optimise_codon_usage(RNA, modi_table, fixed)
  print('*', end='')
  fixed = sorted(fixed)
  fixed = remove_adjascent(fixed, 3, 24)
  RNA = optimise_codon_usage(RNA, modi_table, fixed)
  print('*', end='')
  RNA = remove_repeats(RNA, modi_table, 12, 5, fixed)
  print('*', end='')
  RNA = remove_repeats(RNA, table, 12, 0)
  print('*')

  translated = str(Seq(RNA).translate())
  if protein != translated:
    raise RuntimeError('Protein translation changed')
  forbidden = [motif for motif in seq2avoid if motif in RNA]
  if forbidden:
    raise RuntimeError('Forbidden coding-sequence motifs remain: '+repr(forbidden))
  rare = [RNA[i:i+3] for i in range(0, len(RNA), 3) if RNA[i:i+3] in excluded_codon]
  if rare:
    raise RuntimeError('Excluded rare codons remain: '+repr(sorted(set(rare))))
  result = RNA.replace('U', 'T')
  print('Coding DNA generated: '+str(len(result))+' nt', flush=True)
  print(json.dumps({'codon_distribution_error':codn_usage_efficiency(RNA, modi_table)[0],
                    'fold_cache':rna_sec_structure.cache_info()._asdict(),
                    'seed':args.seed}), flush=True)
  if args.output:
    with open(args.output, 'w') as file:
      file.write(result+'\n')
  print('END', flush=True)
