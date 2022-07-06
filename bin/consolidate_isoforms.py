#!/usr/bin/env python3
import sys, os, glob, argparse

parser = argparse.ArgumentParser(description='options', \
	usage='consolidate_isoforms.py tempdir run_id outfilename [options]')
required = parser.add_argument_group('required named arguments')
required.add_argument('temp_dir', action='store', \
	type=str, help='temp_dir')
required.add_argument('run_id', action='store', \
	type=str, help='output')
parser.add_argument('outbase', action='store', \
	type=str, help='Number of threads (default=4)')
parser.add_argument('--generate_map', action='store', dest='generate_map', \
	required=False, help='''Specify this argument to force overwriting of files in
	an existing output directory''')
args = parser.parse_args()

temp_dir, run_id, outbase = args.temp_dir, args.run_id, args.outbase

# try:
# 	temp_dir = sys.argv[1]
# 	run_id = sys.argv[2]
# 	outbase = sys.argv[3]
# except:
# 	sys.stderr.write('usage: consolidate_isoforms.py tempdir run_id outfilename\n')
# 	sys.stderr.write('script for collapse-range\n')
	# sys.exit(1)

# subprocess.call(['cat']+glob.glob(args.temp_dir+run_id+'*isoforms.fa'),  stdout=open(args.o+'.isoforms.fa', 'w'))
# if args.f:
# 	subprocess.call(['cat']+glob.glob(args.temp_dir+run_id+'*isoforms.gtf'),  stdout=open(args.o+'.isoforms.gtf', 'w'))
# subprocess.call(['rm']+glob.glob(args.temp_dir+run_id+'*'))


f = open(outbase+'.isoforms.bed', 'wt')
bedfiles = glob.glob(temp_dir+run_id+'*isoforms.bed')
for bed in bedfiles:
	f.write(open(bed,'r').read())
f.close()

f = open(outbase+'.isoforms.fa', 'wt')
fafiles = glob.glob(temp_dir+run_id+'*isoforms.fa')
for fa in fafiles:
	f.write(open(fa,'r').read())
f.close()

f = open(outbase+'.isoforms.gtf', 'wt')
gtffiles = glob.glob(temp_dir+run_id+'*isoforms.gtf')
for gtf in gtffiles:
	f.write(open(gtf,'r').read())
f.close()

print(args.generate_map)
if args.generate_map != 'False':
	f = open(outbase+'.map.txt', 'wt')
	mapfiles = glob.glob(temp_dir+run_id+'*map.txt')
	for mf in mapfiles:
		f.write(open(mf,'r').read())
	f.close()

all_temp = glob.glob(temp_dir+run_id+'*')
for t in all_temp:
	os.remove(t)
