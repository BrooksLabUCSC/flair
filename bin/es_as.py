#!/usr/bin/env python3

from __future__ import print_function
import os, sys
import re


class Gene(object):
	'''
	'''
	def __init__(self, geneID=None, chromosome=None, strand=None):
		self.name = geneID
		self.isoforms = dict()
		self.strand = strand
		self.chrom = chromosome
		if self.strand == "+":
			self.acceptor, self.donor = 0, 1
		else:
			self.acceptor, self.donor = 1, 0

	def buildGraph(self):
		self.exonGraph = dict()
		self.knownJuncs = dict()
		self.spliceSites = dict()

		isos = sorted(list(self.isoforms.keys()))
		for i in isos:

			for num,j in enumerate(self.isoforms[i]):
				acceptor,donor = j[self.acceptor],j[self.donor]

				if num == 0:
					# this exon cannot be skipped, so just continue
					continue

				if num+2 > len(self.isoforms[i]):
					# this exon also cannot be skipped, so just continue
					continue


				previousDonor = self.isoforms[i][num-1][self.donor]
				nextAcceptor = self.isoforms[i][num+1][self.acceptor]

				if acceptor not in self.spliceSites:
					self.spliceSites[acceptor] = SpliceSite(acceptor,"acceptor")
				if donor not in self.spliceSites:
					self.spliceSites[donor] = SpliceSite(donor,"donor")
				if nextAcceptor not in self.spliceSites:
					self.spliceSites[nextAcceptor] = SpliceSite(nextAcceptor,"acceptor")
				if previousDonor not in self.spliceSites:
					self.spliceSites[previousDonor] = SpliceSite(previousDonor,"donor")

				acceptor,donor = self.spliceSites[acceptor], self.spliceSites[donor] 
				previousDonor = self.spliceSites[previousDonor]
				nextAcceptor = self.spliceSites[nextAcceptor]

				acceptor.up.add(previousDonor)
				acceptor.down.add(donor)

				donor.up.add(acceptor)
				donor.down.add(nextAcceptor)

				previousDonor.down.add(acceptor)
				nextAcceptor.up.add(donor)

				exon = tuple(sorted([acceptor.name,donor.name]))
				if exon not in self.exonGraph:
					self.exonGraph[exon] = Exon(exon,acceptor,donor)

	def buildGraphv2(self):
		self.exonGraph = dict()
		self.knownJuncs = dict()
		self.spliceSites = dict()

		isos = sorted(list(self.isoforms.keys()))
		for i in isos:
			for num,j in enumerate(self.isoforms[i]):
				acceptor,donor = j[self.acceptor],j[self.donor]

				if num == 1 and len(self.isoforms[i]) == 2:
					# two-exon isoform, only add the splice junction and not the exons
					if acceptor not in self.spliceSites:
						self.spliceSites[acceptor] = SpliceSite(acceptor,"acceptor")
					previousDonor = self.isoforms[i][num-1][self.donor]
					if previousDonor not in self.spliceSites:
						self.spliceSites[previousDonor] = SpliceSite(previousDonor,"donor")
					j1 = (self.spliceSites[previousDonor], self.spliceSites[acceptor])
					if j1 not in self.knownJuncs:
						self.knownJuncs[j1] = set()
					self.knownJuncs[j1].add(i)
					continue
				elif num == 0 or num + 1 >= len(self.isoforms[i]):
					# first or last exons cannot be skipped, continue
					continue

				previousDonor = self.isoforms[i][num-1][self.donor]
				nextAcceptor = self.isoforms[i][num+1][self.acceptor]

				if acceptor not in self.spliceSites:
					self.spliceSites[acceptor] = SpliceSite(acceptor,"acceptor")
				if donor not in self.spliceSites:
					self.spliceSites[donor] = SpliceSite(donor,"donor")
				if nextAcceptor not in self.spliceSites:
					self.spliceSites[nextAcceptor] = SpliceSite(nextAcceptor,"acceptor")
				if previousDonor not in self.spliceSites:
					self.spliceSites[previousDonor] = SpliceSite(previousDonor,"donor")

				acceptor,donor = self.spliceSites[acceptor], self.spliceSites[donor] 
				previousDonor = self.spliceSites[previousDonor]
				nextAcceptor = self.spliceSites[nextAcceptor]
			

				exon = tuple(sorted([acceptor.name,donor.name]))
				if exon not in self.exonGraph:
					self.exonGraph[exon] = Exon(exon,acceptor,donor)

				self.exonGraph[exon].inclusionJuncs.add(((previousDonor,acceptor),(donor,nextAcceptor)))

				j1,j2 = (previousDonor,acceptor), (donor,nextAcceptor)
				if j1 not in self.knownJuncs:
					self.knownJuncs[j1] = set()
				if j2 not in self.knownJuncs:
					self.knownJuncs[j2] = set()
				self.knownJuncs[j1].add(i)
				self.knownJuncs[j2].add(i)

	def findSkippedExonsV3(self):
		'''
		'''
		for i,ss in self.spliceSites.items():
			if ss.ssType != "acceptor":
				continue
			exons = [(i,x) for x in ss]

			for exon in exons:
				associatedAcceptors = exon[-1].down
				associatedDonors = exon[0].up

				allAcceptors = associatedAcceptors.union(associatedDonors.down)
				allDonors = associatedDonors.union(associatedAcceptors.up)

	def findSkippedExonsV1(self):
		'''
		'''

		for e,obj in self.exonGraph.items():
			donor,acceptor = obj.donor, obj.acceptor
			if donor == None or acceptor == None:
				# first or last exon
				continue
			inclusionIsos = set()
			exclusionIsos = set()
			for juncs in obj.inclusionJuncs:
				j1,j2 = juncs
				if j1 in self.knownJuncs and j2 in self.knownJuncs:
					inclusionIsos = inclusionIsos.union(self.knownJuncs[j1].intersection(self.knownJuncs[j2]))
				if (j1[0],j2[-1]) in self.knownJuncs:
					exclusionIsos = exclusionIsos.union(self.knownJuncs[j1[0],j2[-1]])
			print("%s:%s-%s" % (self.chrom,acceptor.name,donor.name), self.strand, len(inclusionIsos), \
				len(exclusionIsos), ",".join(inclusionIsos),",".join(exclusionIsos), sep="\t")
			#print(self.chrom, "\t".join(str(x) for x in sorted([acceptor.name,donor.name])), "%s:%s-%s" % (self.chrom,acceptor.name,donor.name), self.name, self.strand, sep="\t")


	def findSkippedExonsV2(self):
		'''
		'''
		for i,ss in self.spliceSites.items():
			if ss.ssType != "donor":
				continue
			nextSS = ss.down
			if len(nextSS) == 1:
				#single path
				#print("chr6:%s-%s" % (ss.name,list(nextSS)[0].name))
				continue
			else:
				paths=list()
				for j in nextSS:

					paths.append(self.getPaths(ss,j,[]))
				for num,j in num,paths:
					print()

			# for j in nextSS:
			# 	print("chr6:%s-%s" % (ss.name,j.name), sep="\t")

	def getPaths(self, sink=None, currentSS=None, pathSoFar=None):
		'''
		'''
		if len(pathSoFar)<1:
			pathSoFar = [currentSS.name]
		else:
			pathSoFar = [currentSS.name] + pathSoFar

		while currentSS != sink:

			possiblePaths = list(currentSS.up)[0]
			currentSS = possiblePaths
			pathSoFar = [currentSS.name] + pathSoFar
			

		return pathSoFar
	def findSkippedExons(self):
		'''
		'''
		for exon,obj in self.exonGraph.items():

			donor,acceptor = obj.donor, obj.acceptor
			if donor == None or acceptor == None:
				# first or last exon
				continue

			downstreamAcceptors = donor.down
			downAcceptorPathsIn = set([y.name for x in downstreamAcceptors for y in x.up])
			

			upstreamDonors = acceptor.up
			upDonorPathsIn = set([y.name for x in upstreamDonors for y in x.down])
			if len(downAcceptorPathsIn) == 1 and len(upDonorPathsIn) == 1:
				# no alternative routes
				print(exon)
				continue

			print(exon,upDonorPathsIn,downAcceptorPathsIn,sep="\t")
			# upstreamDonors = obj.acceptor.up
			# altAcceptors = set()
			# for donor in upstreamDonors:
			# 	for altAcceptor in donor.down:
			# 		if altAcceptor == obj.acceptor:
			# 			continue
			# 		else:
			# 			altAcceptors.add((altAcceptor,altAcceptor.name))
			
			# downstreamAcceptors = obj.donor.down
			# altDonors = set()
			# for acceptor in downstreamAcceptors:
			# 	for altDonor in acceptor.up:
			# 		if altDonor == obj.donor:
			# 			continue
			# 	else:
			# 		altDonors.add((altDonor,altDonor.name))

			# print(exon,altDonors)#.intersection(upstreamDonors))


class SpliceSite(object):
	'''
	'''
	def __init__(self,ssName = None, ssType = None):
		self.name = ssName
		self.ssType = ssType
		self.up = set()
		self.down = set()

class Exon(object):
	'''
	'''

	def __init__(self, exonID=None, acceptor=None, donor=None):
		self.name =  exonID 
		self.donor = donor
		self.acceptor = acceptor
		self.inclusionJuncs = set()

class Junction(object):
	'''
	'''

	def __init__(self, jcnID= None):
		self.name = jcnID
		self.inExon = set()
		self.outExon = set()
		self.weight = int()

#functions
def bed12toExons(start,starts,sizes):
    '''
    Take bed12 entry and convert block/sizes to exon coordinates.
    '''
    start = int(start)
    sizes, starts = list(map(int,sizes)), list(map(int,starts))
    exons = list()
    for num, st in enumerate(starts,0):
        c1 = st + start
        c2 = c1 + sizes[num]
        exons.append((c1,c2))
    return exons


def parse_gene_id(iso_gene):
	if '_chr' in iso_gene:
		gene = iso_gene[iso_gene.rfind('_chr')+1:]
	elif '_XM' in iso_gene:
		gene = iso_gene[iso_gene.rfind('_XM')+1:]
	elif '_XR' in iso_gene:
		gene = iso_gene[iso_gene.rfind('_XR')+1:]
	elif '_NM' in iso_gene:
		gene = iso_gene[iso_gene.rfind('_NM')+1:]
	elif '_NR' in iso_gene:
		gene = iso_gene[iso_gene.rfind('_NR')+1:]
	elif '_R2_' in iso_gene:
		gene = iso_gene[iso_gene.rfind('_R2_')+1:]
	else:
		gene = iso_gene[iso_gene.rfind('_')+1:]
	return gene

#### main ####


def main():

	flairIsoforms = sys.argv[1]
	genes = dict()
	with open(flairIsoforms) as fin:
		for line in fin:
			cols = line.rstrip().split()
			iso, start, starts, sizes = cols[3], cols[1], cols[11], cols[10]
			chrom,strand = cols[0], cols[5]
			starts = starts.rstrip(",").split(",")
			sizes = sizes.rstrip(",").split(",")

			geneID = parse_gene_id(iso)
			exons = bed12toExons(start,starts,sizes)
			
			if chrom not in genes:
				genes[chrom] = Gene(geneID, chrom, strand)
			geneObj = genes[chrom]
			geneObj.isoforms[iso] = exons
	
	for chrom, gobj in genes.items():
		gobj.buildGraphv2()
		gobj.findSkippedExonsV1()

if __name__ == "__main__":
	main()
