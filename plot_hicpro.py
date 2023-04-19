#!/usr/bin/env python3

import sys
import csv
import argparse
import numpy as np
import pandas as pd
import heapq
import re
import math
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

def parse_args(args):
	parser = argparse.ArgumentParser(description="Check the help flag")
	parser.add_argument("-ch",
						"--chr_sizes",
						dest="sizes",
						help="REQUIRED: chromosome sizes (.chrom.sizes)",
						required=True)
	parser.add_argument("-m",
						"--matrix",
						dest="matrix",
						help="REQUIRED: tab-delimited matrix output from HiC-Pro.",
						required=True)
	parser.add_argument("-s",
						"--sample",
						dest="sample",
						help="Sample name",
						required=False)
	parser.add_argument("-o",
						"--output",
						dest="output",
						help="Output directory.",
						required=False)
	parser.add_argument("-r",
						"--resolution",
						dest="res",
						help="Binning resolution.",
						required=False,
						default=10000)
	parser.add_argument("-n",
						"--normalize",
						dest="norm",
						help="Perform counts-per-million normalization.",
						required=False,
						action="store_false",
						default=None)
	parser.add_argument("-g",
						"--gene_list",
						dest="genes",
						help="List of genes to annotate in gff format.",
						required=False)
	parser.add_argument("-c",
						"--centromeres",
						dest="centromeres",
						help="List of centromere locations.",
						required=False)
	return parser.parse_args()

def input_params(args):
	if args.output is None:
		out = os.path.dirname(args.matrix)
	else:
		out = args.output
	if args.sample is None:
		sample = os.path.split(os.path.dirname(os.path.abspath(args.matrix)))[-1]
	else:
		sample = args.sample
	if args.res is not None:
		res = resolution(args.res)
	return out, sample, res

def resolution(res):
	res = res.lower()
	if "kb" in res:
		return int(res.split("kb")[0]) * 1000
	elif "mb" in res:
		return int(res.split("mb")[0]) * 1000000
	else:
		return int(res)

def get_chrom_starts(sizes):
	start_index = 0
	chrom_starts = {}
	for chrom in sizes.itertuples():
		chrom_starts[chrom[1]] = start_index
		start_index += chrom[2]
	return chrom_starts

def matrix2array(num, mtx, norm):
	mat = np.zeros((num, num))
	with open(mtx, 'r') as m:
		for i in csv.reader(m, delimiter = "\t"):
			if int(i[0]) <= num and int(i[1]) <= num:
				mat[int(i[0]) - 1][int(i[1]) - 1] = i[2]
				mat[int(i[1]) - 1][int(i[0]) - 1] = i[2]
		if norm is not None:
			reads = np.sum(mat)
			mat = mat / (reads * 0.000001)
	return mat

def generate_intrachromosomal_array(chrom_starts,arr,sizes,pdf,sample,res,gl,cent):
	maxColorIntra = []
	mC = 0
	for i in range(0,len(arr)):
		arr[i,i] = 0
	for c in chrom_starts:
		if gl is not None:
			gene_bins = gl[gl[0] == c]
			gene_bins = np.unique(gene_bins[[3,4]]).tolist()
		else:
			gene_bins = None
		if cent is not None:
			cent_bins = cent[cent[0] == c]
			cent_bins = list(range(int(cent_bins[1]),int(cent_bins[2]+1)))
		else:
			cent_bins = None

		chrom_sizes = sizes[sizes[0] == c]
		start = chrom_starts[c]
		end = start + int(chrom_sizes[1])
		m = arr[start:end, start:end]

		if len(str(c)) > 2:
			ch = re.split('_', str(c))[1]
		else:
			ch = str(c)

		if not maxColorIntra:
			allContacts = []
			for i, row in enumerate(m):
				for j, elem in enumerate(row):
					allContacts.append(elem)
			top10percent = int(math.ceil(len(allContacts) * 0.10))
			maxColor = min(heapq.nlargest(top10percent, allContacts))
			plot_intra(pdf,m,maxColor,sample,ch,res,gene_bins,cent_bins)
		else:
			plot_intra(pdf,m,maxColorIntra[mC],sample,ch,res,gene_bins,cent_bins)
		mC += 1

def plot_intra(pdf,m,maxColor,sample,ch,res,gene_bins,cent_bins):
	if gene_bins is not None:
		gene_tick_width = 65/np.shape(m)[0]*4
	if cent_bins is not None:
		cent_tick_width = 65/np.shape(m)[0]*4
	fig, ax = plt.subplots()
	plt.imshow(m, cmap='Reds', origin='lower', vmax=maxColor)
	cbar = plt.colorbar(cmap = 'Reds')
	cbar.ax.set_ylabel('Normalized contact counts', rotation=270, labelpad=12, size=10)
	plt.title("%s chr%s %skb \n intrachromosomal interactions" % (sample, ch, res // 1000), pad=12)
	if gene_bins is not None:
		secaxx = ax.secondary_xaxis('top')
		secaxx.set_xticks(gene_bins)
		secaxx.tick_params(width=gene_tick_width, length=8, color='#009e00', labeltop=False)
		secaxy = ax.secondary_yaxis('right')
		secaxy.set_yticks(gene_bins)
		secaxy.tick_params(width=gene_tick_width, length=8, color='#009e00', labelright=False)
	if cent_bins is not None:
		secaxx.set_xticks(cent_bins, minor=True)
		secaxx.tick_params(which='minor', width=cent_tick_width, length=8, color='#0000ff', labeltop=False)
		secaxy.set_yticks(cent_bins, minor=True)
		secaxy.tick_params(which='minor', width=cent_tick_width, length=8, color='#0000ff', labelright=False)
	plt.tight_layout()
	pdf.savefig()
	plt.close()

def generate_interchromosomal_array(chrom_starts,arr,sizes,pdf,sample,res):
	maxColorInter = 10.679982
	chrom_sizes = {}
	for chrom in sizes.itertuples():
		chrom_sizes[chrom[1]] = chrom[2]
	for key, val in chrom_starts.items():
		arr[chrom_starts[key]:(chrom_starts[key]+chrom_sizes[key]),
		chrom_starts[key]:(chrom_starts[key]+chrom_sizes[key])] = 0
	if not maxColorInter:
		allContacts = []
		for i, row in enumerate(arr):
			for j, elem in enumerate(row):
				allContacts.append(elem)
		top10percent = int(math.ceil(len(allContacts) * 0.10))
		maxColor = min(heapq.nlargest(top10percent, allContacts))
		plot_inter(pdf,arr,maxColor,sample,res)
	else:
		plot_inter(pdf,arr,maxColorInter,sample,res)

def plot_inter(pdf,arr,maxColor,sample,res):
	plt.imshow(arr, cmap = 'Reds', origin = 'lower', vmax = maxColor)
	cbar = plt.colorbar(cmap = 'Reds')
	cbar.ax.set_ylabel('Normalized contact counts', rotation = 270, labelpad = 12, size = 10)
	plt.title("%s %skb \n interchromosomal interactions" % (sample, res // 1000))
	plt.tight_layout()
	pdf.savefig()
	plt.close()

def main():
	args = parse_args(sys.argv[1:])
	out, sample, res = input_params(args)

	sizes = pd.read_csv(args.sizes, sep='\t', header=None)
	sizes[1] = sizes[1] // res + 1
	
	num_bins = np.sum(sizes[1])

	chrom_starts = get_chrom_starts(sizes)

	arr = matrix2array(num_bins, args.matrix, args.norm)

	if args.genes:
		gl = pd.read_csv(args.genes, sep="\t", header=None)
		gl = gl[[0,3,4]]
		gl[3] = gl[3]//res
		gl[4] = gl[4]//res
	else:
		gl = None

	if args.centromeres:
		cent = pd.read_csv(args.centromeres, sep="\t", header=None)
		cent[1] = cent[1]//res
		cent[2] = cent[2]//res
	else:
		cent = None

	pdf = PdfPages('%s/%s_%skb.pdf' % (out, sample, res // 1000))
	generate_intrachromosomal_array(chrom_starts,arr,sizes,pdf,sample,res,gl,cent)
	generate_interchromosomal_array(chrom_starts,arr,sizes,pdf,sample,res)
	pdf.close()

if __name__ == "__main__":
	main()