import argparse
from flair.pycbio.sys import fileOps


def parse_args():
    parser = argparse.ArgumentParser(description='''for annotating FLAIR aaseq predictions with UniProt names''')
    parser.add_argument('-i', '--input_aaseq', required=True, help='protein sequence fasta file - sequence names should be aaseqID_geneID')
    parser.add_argument('-r', '--reference_seq', required=True, help='reference protein sequence fasta file - sequence names should be db|refID|geneName_organism, ex: >sp|Q8NH21|OR4F5_HUMAN. This script will extract the refID. this can be fasta or fasta.gz')
    parser.add_argument('-o', '--output', required=True, help='output name - should be a fasta file')
    args = parser.parse_args()
    return args


def parse_fasta(filename):
    last, seq = None, []
    for line in fileOps.opengz(filename):
        if line[0] == '>':
            if last:
                yield last, ''.join(seq)
                seq = []
            last = line[1:].rstrip()
        else:
            seq.append(line.rstrip())
    if last:
        yield last, ''.join(seq)


def process_reference(ref_file):
    ref_seq_to_name = {}
    for info, seq in parse_fasta(ref_file):
        id = info.split('|')[1]
        ref_seq_to_name[seq] = id
    return ref_seq_to_name


def annotate_input(input_aaseq, ref_seq_to_name, output_name):
    out = open(output_name, 'w')

    for info, seq in parse_fasta(input_aaseq):
        id, otherinfo = info.split('_', 1)
        if seq in ref_seq_to_name:
            id = ref_seq_to_name[seq]
        info = id + '_' + otherinfo
        out.write('>' + info + '\n' + seq + '\n')
    out.close()

def main():
    args = parse_args()
    ref_seq_to_name = process_reference(args.reference_seq)
    annotate_input(args.input_aaseq, ref_seq_to_name, args.output)

if __name__ == '__main__':
    main()