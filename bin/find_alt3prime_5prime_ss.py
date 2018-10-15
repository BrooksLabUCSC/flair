import sys, csv
import scipy.stats as sps

try:
	psl = open(sys.argv[1])
	colnum = int(sys.argv[2])
	outfilename3 = sys.argv[3]
	outfilename5 = sys.argv[4]
except:
	print('usage: script.py isoforms.psl colnum alt_acceptor.txt alt_donor.txt')
	sys.exit()

def get_junctions_psl(starts, sizes):
	junctions = set()
	if len(starts) != 1:
		for b in range(len(starts)-1):
			junctions.add((starts[b]+sizes[b], starts[b+1]))
		return junctions

def pslreader(psl, index=0, fiveprimeon=False):
	junctiondict = {}
	for line in psl:
		line = line.rstrip().split('\t')
		chrom, name, strand, start, end = line[13], line[9], line[8], int(line[15]), int(line[16])
		count0, count1 = int(line[colnum]), int(line[colnum + 1])
		chrom = strand + chrom  # stranded comparisons
		sizes = [int(n) for n in line[18].split(',')[:-1]]
		starts = [int(n) for n in line[20].split(',')[:-1]]
		junctions = get_junctions_psl(starts, sizes)
		if chrom not in junctiondict:
			junctiondict[chrom] = {}
		for j in junctions:
			if fiveprimeon:
				fiveprime, threeprime = j[1], j[0]
			else:
				fiveprime, threeprime = j[0], j[1]
			if fiveprime not in junctiondict[chrom]:
				junctiondict[chrom][fiveprime] = {}  # 5' end anchor
			if threeprime not in junctiondict[chrom][fiveprime]:
				junctiondict[chrom][fiveprime][threeprime] = [0,0, name]
			junctiondict[chrom][fiveprime][threeprime][0] += float(count0)
			junctiondict[chrom][fiveprime][threeprime][1] += float(count1)
	return junctiondict

def find_altss(alljuncs, writer, fiveprimeon=False):  # comments are for alt3' splicing
	far = 0
	maxdist = 100
	nocanonical = 0
	totalnumsites = 0
	for chrom in alljuncs:
		for fiveprime in alljuncs[chrom]:
			if len(alljuncs[chrom][fiveprime]) == 1:  # if there is only one 3' end, not an alt 3' junction
				continue
			entries = []
			for threeprime1 in alljuncs[chrom][fiveprime]:  # for each potential cryptic ss
				thisjunc = alljuncs[chrom][fiveprime][threeprime1][:2]
				if thisjunc[0] < 3 and thisjunc[1] < 3:
					continue
				allothercounts = [0,0]  # for junctions anchored with this 5' site
				oro3p = ('', 0)  # overrepresented other 3 prime site, used to calculate distance from cryptic SS
				farused = False
				for threeprime2 in alljuncs[chrom][fiveprime]:
					if threeprime1 == threeprime2:
						continue
					if abs(int(threeprime1) - int(threeprime2)) > 200:  # likely an exon skipping and not a alt 3/5' site
						if not farused:
							far += 1
							farused = True
						continue
					if int(alljuncs[chrom][fiveprime][threeprime2][0]) > oro3p[1]:
						oro3p = (threeprime2, int(alljuncs[chrom][fiveprime][threeprime2][0]))
					allothercounts[0] += alljuncs[chrom][fiveprime][threeprime2][0]
					allothercounts[1] += alljuncs[chrom][fiveprime][threeprime2][1]				
				samp1coverage = thisjunc[0]+allothercounts[0]
				samp2coverage = thisjunc[1]+allothercounts[1]
				if samp1coverage < 10 or samp2coverage < 10:  # if not enough alternative ss usage
					continue
				ctable = [thisjunc, allothercounts]
				name = alljuncs[chrom][fiveprime][threeprime1][2]
				if oro3p[1] == 0:
					continue
				totalnumsites += 1
				if chrom[0] == '-':
					diff = int(threeprime1)-int(oro3p[0])
					if fiveprimeon:
						entries += [[chrom[1:], fiveprime, threeprime1, sps.fisher_exact(ctable)[1], 
						chrom[0]] + ctable[0] + ctable[1] + [name] + [diff] + [oro3p[0]]]
					else:
						entries += [[chrom[1:], threeprime1, fiveprime, sps.fisher_exact(ctable)[1], 
						chrom[0]] + ctable[0] + ctable[1] + [name] + [diff] + [oro3p[0]]]
				else:
					diff = int(oro3p[0])-int(threeprime1)
					if fiveprimeon:
						entries += [[chrom[1:], threeprime1, fiveprime, sps.fisher_exact(ctable)[1], 
						chrom[0]] + ctable[0] + ctable[1] + [name] + [diff] + [oro3p[0]]]
					else:
						entries += [[chrom[1:], fiveprime, threeprime1, sps.fisher_exact(ctable)[1], 
						chrom[0]] + ctable[0] + ctable[1] + [name] + [diff] + [oro3p[0]]]
			if entries:
				entries = sorted(entries, key=lambda x: x[10])
				entries = sorted(entries, key=lambda x: x[3])
				for e in entries:
					writer.writerow(e)
					break

alljuncs = pslreader(psl)
with open(outfilename3, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	find_altss(alljuncs, writer)


alljuncs = pslreader(open(sys.argv[1]), fiveprimeon=True)
with open(outfilename5, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	find_altss(alljuncs, writer, fiveprimeon=True)






