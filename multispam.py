# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details at
# http://www.gnu.org/copyleft/gpl.html

def check_positive(value):
	import argparse
	ivalue = int(value)
	if ivalue <= 0:
		raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
	return ivalue


def build_opts():
	import argparse
	parser = argparse.ArgumentParser(description = '''Build phylogenetic trees from DNA sequence data with Maximum-Likelihood
		 and a quartet based super tree approach''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", nargs="?", type=str, help="Input file (multi-fasta file)", required=True)
	parser.add_argument("-o", nargs="?", type=str, help="Output file (tree file in newick-format)", default="outtree")
	parser.add_argument("-k", metavar="w", nargs="?", type=int, choices=range(6,17), help="weight of the pattern \
		(the number of matching positions)", default=10)
	parser.add_argument("-d", nargs="?", type=check_positive, help="number of don't care positions in the pattern", default=100)
	parser.add_argument("-n", nargs="?", type=check_positive, help="number of sampled quartet blocks", default=1000000)
	parser.add_argument("-t", nargs="?", type=check_positive, help="number of threads", default=1)
	parser.add_argument("-m", metavar="--mem-save", nargs="?", type=bool, \
		const=True, default=False, help="memory saving mode")
	parser.add_argument('--version', action='version', version='Multi-SpaM 1.0')

	return parser.parse_args()

def build_arg_list(opts, tmp_file):
	argl = ["./bin/multi-SpaM"]
	if opts.m == True:
		argl.append("-m")
	del opts.m
	for opt in vars(opts):
		argl.append("-" + opt)
		argl.append(str(getattr(opts,opt)))
	argl[next(idx for idx, arg in enumerate(argl) if arg == "-o") + 1] = str(tmp_file)
	return argl

def build_id_map(filename):
	id_map = {}
	last_id = 0
	with open(opts.i) as fasta:
		for line in fasta:
			if line.startswith(">"):
				id_map[last_id] = line[1:-1]
				last_id += 1
	return id_map


def fix_tree(tree_file, id_map, out_file):
	import re
	print(id_map)
	with open(tree_file) as file, open(out_file, "w+") as out:
		for line in file:
			parts = re.split('(,|\(|\)|;)',line)
			for idx, p in enumerate(parts):
				if p.isdigit():
					value = id_map[int(p)]
					parts[idx] = value
			out.write(''.join(parts))


if __name__ == "__main__":
	import sys
	import subprocess as sp
	import uuid
	import os

	id = str(uuid.uuid4())
	if not (sys.platform == "linux" or sys.platform == "linux2") or not (sys.maxsize > 2 ** 32):
		raise Exception("Multi-SpaM currently works only on 64-bit linux!")
	opts = build_opts()

	try:
		open("./bin/multi-SpaM").close()
	except IOError:
		raise Exception("Multi-SpaM was not installed yet. Run \"make\", then run the script again.")
	
	quartet_file = id + ".quartets"
	argl = build_arg_list(opts, quartet_file)

	sp.check_call(argl)

	id_map = build_id_map(opts.i)
	id_tree_file = id + ".tree"

	sp.check_call(["./bin/max-cut-tree", "qrtt=" + quartet_file, "weights=off", "otre=" + id_tree_file])

	fix_tree(id_tree_file, id_map, opts.o)

	os.remove(quartet_file)
	os.remove(id_tree_file)