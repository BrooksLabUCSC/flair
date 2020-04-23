#!/usr/bin/env python3
import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	genome = open(sys.argv[3])
except:
	sys.stderr.write('usage: script.py reads.psl annotation.gtf genome.fa > reads.productivity.psl \n')
	sys.stderr.write('assumes the psl has the correct strand information\n')
	sys.exit()

aa_dict = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"Z", "UAG":"Z",  # using Z for stop
    "UGU":"C", "UGC":"C", "UGA":"Z", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",} 

def revcomp(seq):
	seq = seq.replace('A', 'X').replace('T', 'A').replace('X', 'T')
	seq = seq.replace('G', 'X').replace('C', 'G').replace('X', 'C')
	return seq[::-1]

def translate_seq(seq):
	seq = seq.replace('T', 'U')
	aa_string = ''
	for i in range(0,len(seq),3):
		if i + 3 > len(seq):
			break
		if 'N' in seq[i:i+3]:
			aa_string += 'O'
		else:
			aa_string += aa_dict[seq[i:i+3]]
	return aa_string

def get_sequence(entry, seq, strand='+'):
	blocksizes = [int(x) for x in entry[18].split(',')[:-1]]
	blockstarts = [int(x) for x in entry[20].split(',')[:-1]]
	if strand == '+':  # dont include the last exon per the rule
		blocksizes = blocksizes[:-1]
		blockstarts = blockstarts[:-1]
	else:
		blocksizes = blocksizes[1:]
		blockstarts = blockstarts[1:]
	pulledseq = ''
	for block in range(len(blockstarts)):
		pulledseq += seq[blockstarts[block]:blockstarts[block]+blocksizes[block]]
	return pulledseq.upper()

def get_sequence_after_pos(entry, seq, pos, strand='+'):
	blocksizes = [int(x) for x in entry[18].split(',')[:-1]]
	blockstarts = [int(x) for x in entry[20].split(',')[:-1]]

	pulledseq = ''
	if strand == '+':
		if len(blocksizes) > 1:
			blocksizes = blocksizes[:-1]  # because of the last exon exon junction rule
			blockstarts = blockstarts[:-1]  # this also prevents assessment of single exon transcripts
		pull = False
		for start, size in zip(blockstarts, blocksizes):
			if pull:
				pulledseq += seq[start:start+size]
			elif start <= pos and pos <= size + start:
				pull = True
				pulledseq += seq[pos:start+size]
	else:
		blocksizes = blocksizes[1:]
		blockstarts = blockstarts[1:]
		for start, size in zip(blockstarts, blocksizes):
			if start <= pos and pos <= size + start:
				pulledseq += seq[start:pos+2]
				break
			pulledseq += seq[start:start+size]

	return pulledseq.upper()

def longest_orf(protein):
	if 'M' not in protein:
		return 0, protein
	protm = protein[protein.find('M'):]  # find first M
	if 'Z' not in protm:
		return len(protm), protm
	protz = protm[:protm.find('Z')]  # find first stop after M
	if protm[-1] == 'Z':
		protafterz = ''
	else:
		protafterz = protm[protm.find('Z')+1]
	if 'M' in protafterz:
		if 'Z' in protafterz:
			protafterzz = protafterz[:protafterz.find('Z')]
		else:
			protafterzz = protafterz
		if len(protafterzz) > len(protz):
			sys.stderr.write('did it\n')
			return len(protafterzz), protafterzz
	return len(protz), protz

def seq_to_longest_prot(seq):
	proteins = translate_seq(seq), translate_seq(seq[1:]), translate_seq(seq[2:])
	orflens = [longest_orf(proteins[0]), longest_orf(proteins[1]), longest_orf(proteins[2])]
	# bestprot = sorted(zip(proteins, orflens), key=lambda k: k[1][0])[-1]
	bestprot = sorted(orflens, key=lambda k: k[0])[-1]
	return bestprot[1], bestprot[0]  # tuple of protein and length

def find_tss_pos(entry, tss):  # tss is a dictionary
	entrys, entrye, strand, chrom = int(entry[15]), int(entry[16]), entry[8], entry[13]
	pos_candidates = []
	blocksizes = [int(x) for x in entry[18].split(',')[:-1]]
	blockstarts = [int(x) for x in entry[20].split(',')[:-1]]
	if chrom not in tss:
		return -1
	oldstrand = strand
	if oldstrand == '.':
		strand = '-'  # try both strands, start with '-' first
	elif strand not in tss[chrom]:  # there is assigned strand but no start codon on that chromosome strand
		return -1

	for pos in tss[chrom][strand]:
		if entrys <= pos and pos <= entrye:
			for size, start in zip(blocksizes, blockstarts):
				if start <= pos and pos <= start + size: 
					pos_candidates += [pos]
	if oldstrand == '.':
		strand = '+'  # try other strand
		for pos in tss[chrom][strand]:
			if entrys <= pos and pos <= entrye:
				for size, start in zip(blocksizes, blockstarts):
					if start <= pos and pos <= start + size: 
						pos_candidates += [pos]
	if not pos_candidates:
		# sys.stderr.write('no candidates {}:{}-{}\n'.format(chrom, entrys, entrye))
		return -1 
	elif oldstrand == '.':
		return pos_candidates
	elif len(pos_candidates) == 1:
		return pos_candidates[0]
	elif strand == '+':
		return min(pos_candidates)
	return max(pos_candidates)

psldata = {}
for entry in psl:
	entry = entry.rstrip().split('\t')
	if entry[13] in psldata:  # line[13] is chromosome
		psldata[entry[13]] += [entry]
	else:
		psldata[entry[13]] = [entry]

tss = {}  # start codons
for line in gtf:
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'start_codon':
		continue
	if chrom not in tss:
		tss[chrom] = {}
	if strand not in tss[chrom]:
		tss[chrom][strand] = []
	tss[chrom][strand] += [int(start)]

unproductive = 0
valid_transcripts = 0.
seq, chrom = '', ''
noM = 0  # no start codon

for line in genome:
	line = line.rstrip()
	if line.startswith('>'):
		if not chrom:
			chrom = line.split()[0][1:]
			continue
		if chrom in psldata:
			for entry in psldata[chrom]:
				strand = entry[8]
				pos = find_tss_pos(entry, tss)
				if pos == -1:
					print('\t'.join(entry+['2']))
					continue
				if strand == '-':
					revseq = revcomp(get_sequence_after_pos(entry, seq, pos, '-'))
					bestprot = translate_seq(revseq)
				elif strand == '+':
					forseq = get_sequence_after_pos(entry, seq, pos-1)
					bestprot = translate_seq(forseq)
				elif strand == '.':  # infer strand based on ORF length...
					revseq = revcomp(get_sequence_after_pos(entry, seq, max(pos), '-'))
					bestprotrev = translate_seq(revseq)
					forseq = get_sequence_after_pos(entry, seq, min(pos)-1)
					bestprotfor = translate_seq(forseq)
					if len(bestprotrev) > len(bestprotfor):
						bestprot = bestprotrev
						entry[8] = '-'
					else:
						bestprot = bestprotfor
						entry[8] = '+'

				if not bestprot:  # single exon transcript
					print('\t'.join(entry+['2']))  # single exon transcript, can't be assessed using existing method
					continue
				if not bestprot or bestprot[0] != 'M':
					noM += 1
					print('\t'.join(entry+['2']))  # lncRNA of sorts
					continue
				valid_transcripts += 1
				bestprot = bestprot[bestprot.find('M'):]
				protquery = bestprot[:-int(55/3)]  # everything up from N terminal to 55 nt upstream of last exon-exon nucleotide
				if 'Z' in protquery:  # Z meaning a stop codon
					print('\t'.join(entry+['1']))
					unproductive += 1
					continue
				print('\t'.join(entry+['0']))
		chrom = line.split()[0][1:]
		seq = ''
	else:
		seq += line

# sys.stderr.write('# reads with no M: ' + str(noM)+'\n')
# sys.stderr.write('valid transcripts # ' + str(valid_transcripts)+'\n')
sys.stderr.write('Unproductive proportion estimate ' + str(unproductive / (valid_transcripts + unproductive)) + '\n')

