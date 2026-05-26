import argparse

parser=argparse.ArgumentParser(description='Program description')
parser.add_argument('--msa', required=True, help='msa to delete positions with more than a percentage of - (gaps)')
parser.add_argument('--out_si', required = True, help='output file to write the new msa that satisfyes the th')
parser.add_argument('--out_no', required=True, help='output file to write the sequences that do not satisfy the th')
parser.add_argument('--len', required=True, help='Length of the sequence that will be the input of AF (e.g. 1RGS)')
args=parser.parse_args()

longitud = int(args.len)

with open(args.msa, 'r') as f:
  lines = f.readlines()


for i, line in enumerate(lines):
  if i % 2 == 1:
    if len(line.strip()) == longitud:
      with open(args.out_si, 'a') as new:
        new.writelines(lines[i-1:i+1])
    else:
      with open(args.out_no, 'a') as new:
        new.writelines(lines[i-1:i+1])