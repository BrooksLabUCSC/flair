import sys



alignedbedfile = sys.argv[1]
referencegtffile = sys.argv[2]
outfilename = sys.argv[3]



fusiontoannotsj = {}
for line in open(referencegtffile):
    line = line.rstrip().split('\t', 5)
    if line[2] == 'exon':
        if line[0] not in fusiontoannotsj: fusiontoannotsj[line[0]] = set()
        fusiontoannotsj[line[0]].add(int(line[3]))
        fusiontoannotsj[line[0]].add(int(line[4]))

splicejunctosupport = {}
for line in open(alignedbedfile):
    line = line.rstrip().split('\t')
    thischr, iso, dir, start, esizes, estarts = line[0], line[3], line[5], int(line[1]), \
                            [int(x) for x in line[10].rstrip(',').split(',')], [ int(x) for x in line[11].rstrip(',').split(',')]
    for i in range(len(esizes) - 1):
        thisintron = [thischr, start + estarts[i] + esizes[i] , start + estarts[i + 1] + 1]
        if thischr in fusiontoannotsj:
            for i in range(1,3):
                for sj in fusiontoannotsj[thischr]:
                    if abs(thisintron[i]-sj) <= 3:
                        thisintron[i] = sj
                        break
        thisintron[2] -= 1
        thisintron = tuple(thisintron)
        if thisintron not in splicejunctosupport: splicejunctosupport[thisintron] = 0
        splicejunctosupport[thisintron] += 1

goodsj, lowsupsj = {}, {}

for j in splicejunctosupport:
    if splicejunctosupport[j] >= 2:
        goodsj[j] = splicejunctosupport[j]
    else: lowsupsj[j] = splicejunctosupport[j]


for chr, start, end in lowsupsj:
    for i in range(-2, 3):
        for j in range(-2, 3):
            if (chr, start+i, end+j) in goodsj:
                goodsj[(chr, start+i, end+j)] += lowsupsj[(chr, start, end)]




out = open(outfilename, 'w')
for j in goodsj:
    # if j[0] == 'TP53*ENSG00000141510.18.chr17.7671131.7688547--chr6-30221000.chr6.30223661.30220450':
    #     print(goodsj[j], j)
    if goodsj[j] >= 3:
        out.write('\t'.join([str(x) for x in j]) + '\t+\n')
out.close()
