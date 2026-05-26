import argparse

parser=argparse.ArgumentParser(description='Pass the CSV file to a FASTA format so it can be used in the pipeline')
parser.add_argument('--csv', required=True, help='Input CSV file with protein sequence')
parser.add_argument('--fasta', required=True, help='Output file to write the FASTA format')
args=parser.parse_args()

with open(args.csv, 'r') as csv_file:
    lines = csv_file.readlines()

second_line = lines[1].strip().split(',')
with open(args.fasta, 'w') as fasta_file:
    fasta_file.write(f'>{second_line[0]}\n')
    fasta_file.write(f'{second_line[1]}\n')