import sys
from flair.pycbio.hgdata.bed import BedReader

expected, output = sys.argv[1], sys.argv[2]


def parse_bed(file):
    sjc_to_ends = {}
    for bed in BedReader(file, fixScores=True):
        junctions = tuple((bed.blocks[i].end, bed.blocks[i + 1].start) for i in range(len(bed.blocks) - 1))
        if junctions not in sjc_to_ends:
            sjc_to_ends[junctions] = {}
        sjc_to_ends[junctions][(bed.chromStart, bed.chromEnd)] = bed.score
    return sjc_to_ends


def parse_dict_keys_to_shared(a, b):
    ak = set(a.keys())
    bk = set(b.keys())

    shared = ak & bk
    a_only = ak - bk
    b_only = bk - ak
    return shared, a_only, b_only


expected_sjc_to_ends = parse_bed(expected)
output_sjc_to_ends = parse_bed(output)

shared_sjc, expected_sjc, output_sjc = parse_dict_keys_to_shared(expected_sjc_to_ends, output_sjc_to_ends)

s_ends, e_ends, o_ends = 0, 0, 0
e_sjc, o_sjc = 0, 0
for sjc in shared_sjc:
    shared_ends, expected_ends, output_ends = parse_dict_keys_to_shared(expected_sjc_to_ends[sjc], output_sjc_to_ends[sjc])
    for s, e in shared_ends:
        s_ends += 1
        ec, oc = expected_sjc_to_ends[sjc][(s, e)], output_sjc_to_ends[sjc][(s, e)]
        print(f'shared SJC and ends: e_s:{s} o_s:{s} e_e:{e} o_e:{e} exp_c:{ec} out_c:{oc} {sjc}')

    if len(expected_ends) == 1 and len(output_ends) == 1:
        e_s, e_e = list(expected_ends)[0]
        o_s, o_e = list(output_ends)[0]
        if abs(e_s - o_s) < 300 and abs(e_e - o_e) < 300:
            s_ends += 1
            print(f'shared SJC and ends: e_s:{e_s} o_s:{o_s} e_e:{e_e} o_e:{o_e} exp_c:{ec} out_c:{oc} {sjc}')
        else:
            e_ends += 1
            o_ends += 1
            print(f'shared SJC different ends: e_s:{e_s} o_s:{o_s} e_e:{e_e} o_e:{o_e} exp_c:{ec} out_c:{oc} {sjc}')

    else:
        for s, e in expected_ends:
            e_ends += 1
            ec = expected_sjc_to_ends[sjc][(s, e)]
            print(f'shared SJC, ends expected only: s:{s} e:{e} exp_c:{ec} {sjc}')
        for s, e in output_ends:
            o_ends += 1
            oc = output_sjc_to_ends[sjc][(s, e)]
            print(f'shared SJC, ends output only: s:{s} e:{e} out_c:{oc} {sjc}')

for sjc in expected_sjc:
    for s, e in expected_sjc_to_ends[sjc]:
        e_sjc += 1
        ec = expected_sjc_to_ends[sjc][(s, e)]
        print(f'expected only SJC and ends: s:{s} e:{e} exp_c:{ec} {sjc}')

for sjc in output_sjc:
    for s, e in output_sjc_to_ends[sjc]:
        o_sjc += 1
        oc = output_sjc_to_ends[sjc][(s, e)]
        print(f'output only SJC and ends: s:{s} e:{e} out_c:{oc} {sjc}')

print(f'shared SJC and ends: {s_ends} isoforms')
print(f'shared SJC, ends expected only: {e_ends} isoforms')
print(f'shared SJC, ends output only: {o_ends} isoforms')
print(f'expected only SJC and ends: {e_sjc} isoforms')
print(f'expected only SJC and ends: {o_sjc} isoforms')
