import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors as colors
import os
import math
import heapq

res = 10000

sample = sys.argv[1]

sizes = pd.read_csv(sys.argv[2], sep='\t', header=None)
sizes[2] = sizes[1] // res + 1
inter_ticks = []
prev_size = 0
for i in sizes.itertuples():
	inter_ticks.append(i[3] + prev_size)
	prev_size += i[3]

df = pd.read_csv(sys.argv[3], sep='\t')
df["bin1"] = df["bin1"] // res
df["bin2"] = df["bin2"] // res

if len(sys.argv) > 4:
	cent = pd.read_csv(sys.argv[4], sep="\t", header=None)
	cent[1] = cent[1]//res
	cent[2] = cent[2]//res

	gl = pd.read_csv(sys.argv[5], sep="\t", header=None)
	gl = gl[[0,3,4]]
	gl[3] = gl[3]//res
	gl[4] = gl[4]//res

pdf = PdfPages('%s_%skb_diff.pdf' % (sample, res // 1000))
df_inter = df[df.chr1 != df.chr2]
if len(df_inter) > 1:
	sizes_inter = sizes.copy()
	starts = []
	ends = []
	start_index = 0
	for i in sizes_inter.itertuples():
		starts.append(start_index)
		start_index += i[3]
	sizes_inter[1] = starts
	sizes_inter[2] = sizes_inter[2] + sizes_inter[1] - 1

	sizes_inter.columns = ["chr1", "start", "end"]
	df_inter = df[df.chr1 != df.chr2]
	df_inter = df_inter.merge(sizes_inter, on="chr1")
	df_inter["bin1"] = df_inter["bin1"] + df_inter["start"]
	df_inter.drop(["start", "end"], axis=1, inplace=True)

	sizes_inter.columns = ["chr2", "start", "end"]
	df_inter = df_inter.merge(sizes_inter, on="chr2")
	df_inter["bin2"] = df_inter["bin2"] + df_inter["start"]
	df_inter.drop(["start", "end"], axis=1, inplace=True)
	df_inter = df_inter[df_inter.chr1 != df_inter.chr2]

	maxColor = min(heapq.nlargest(math.ceil(len(df_inter) * 0.10), abs(df_inter.log2FC)))

	num = np.max(sizes_inter.end)
	arr = np.zeros((num,num))
	for row in df_inter.itertuples():
		arr[int(row.bin1) - 1][int(row.bin2) - 1] = row.log2FC
		arr[int(row.bin2) - 1][int(row.bin1) - 1] = row.log2FC
	for j in range(0,len(arr)):
		arr[j,j] = 0
	norm = colors.TwoSlopeNorm(vmin=-maxColor, vcenter=0, vmax=maxColor)
	fig, ax = plt.subplots()
	plt.imshow(arr, cmap = 'RdBu_r', origin = 'lower', norm=norm)
	cbar = plt.colorbar(cmap = 'RdBu_r')
	ax.set_xticks(inter_ticks, minor=True)
	ax.set_yticks(inter_ticks, minor=True)
	ax.tick_params(which="minor", length=0, labeltop=False)
	ax.tick_params(which="minor", length=0, labeltop=False)
	ax.xaxis.grid(True, linestyle="dotted", linewidth=0.5, color="black", which="minor")
	ax.yaxis.grid(True, linestyle="dotted", linewidth=0.5, color="black", which="minor")
	plt.tight_layout()
	pdf.savefig()
	plt.close()

df_intra = df[df.chr1 == df.chr2]
ch = 1
for i in df_intra.chr1.unique():
	sdf = df_intra[df_intra.chr1 == i]

#	maxColor = min(heapq.nlargest(math.ceil(len(sdf) * 0.10), abs(sdf.log2FC)))
	maxColor = max(sdf.log2FC)

	num = int(sizes.loc[sizes[0] == i][2].values)
	arr = np.zeros((num,num))
	for row in sdf.itertuples():
		arr[int(row.bin1) - 1][int(row.bin2) - 1] = row.log2FC
		arr[int(row.bin2) - 1][int(row.bin1) - 1] = row.log2FC
	for j in range(0,len(arr)):
		arr[j,j] = 0

	norm = colors.TwoSlopeNorm(vmin=-maxColor, vcenter=0, vmax=maxColor)

	fig, ax = plt.subplots()
	plt.imshow(arr, cmap = 'RdBu_r', origin = 'lower', norm=norm)
	cbar = plt.colorbar(cmap = 'RdBu_r')
	cbar.ax.set_ylabel('log$_\mathrm{2}$ fold change', rotation=270, labelpad=12, size=10)
	plt.title("%s chr%s %skb\ndifferential intrachromosomal interactions" % (sample, ch, res // 1000), pad=12)
	
	try:
		cent_gene_tick_width = 65/len(arr)*4
		gene_bins = gl[gl[0] == i]
		gene_bins = np.unique(gene_bins[[3,4]]).tolist()
	except:
		pass

	try:
		cent_gene_tick_width = 65/len(arr)*4
		cent_bins = cent[cent[0] == i]
		cent_bins = list(range(int(cent_bins[1]),int(cent_bins[2]+1)))
		secaxx = ax.secondary_xaxis('top')
		secaxx.set_xticks(gene_bins)
		secaxx.tick_params(width=cent_gene_tick_width, length=8, color='#009e00', labeltop=False)
		secaxx.set_xticks(cent_bins, minor=True)
		secaxx.tick_params(which='minor', width=cent_gene_tick_width, length=8, color='#808080', labeltop=False)

		secaxy = ax.secondary_yaxis('right')
		secaxy.set_yticks(gene_bins)
		secaxy.tick_params(width=cent_gene_tick_width, length=8, color='#009e00', labelright=False)
		secaxy.set_yticks(cent_bins, minor=True)
		secaxy.tick_params(which='minor', width=cent_gene_tick_width, length=8, color='#808080', labelright=False)
	except:
		pass
	
	ch += 1
	plt.tight_layout()
	pdf.savefig()
	plt.close()
pdf.close()
