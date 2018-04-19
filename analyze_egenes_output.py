#! /home/yhua295/.linuxbrew/bin/python3
import sys, os
import operator
import pandas as pd 
import argparse
import numpy as np
import gzip

def find_top_n_eqtl(f_In, n=1000):
	"""
	get top n strongest gene-variants cis-eqtl pairs
	use q-value as the criteria
	"""
	genes_eqtl_pairs = pd.read_csv(f_In, index_col=0, sep="\t")
	genes_vars_pair = dict(genes_eqtl_pairs['variant_id'])
	genes_qvals_pair = dict(genes_eqtl_pairs['qval'])
	genes_qvals_pair_list = list(genes_qvals_pair.items())
	genes_qvals_pair_list.sort(key=operator.itemgetter(1))
	top_n_genes = genes_qvals_pair_list[:n]
	top_n_genes_vars_pair = dict((k[0], genes_vars_pair.get(k[0]))for k in top_n_genes)
	return top_n_genes_vars_pair

def find_eqtl_according_to_qval(f_In, neglogq=1):
	"""
	get gene-variants cis-eqtl pairs below threshold
	use q-value as the criteria
	"""
	qval = 10**(-neglogq)
	genes_eqtl_pairs = pd.read_csv(f_In, index_col=0, sep="\t")
	genes_eqtl_pairs = genes_eqtl_pairs[genes_eqtl_pairs['qval'] < qval] 
	significant_genes_vars_pair = dict(genes_eqtl_pairs['variant_id'])
	return significant_genes_vars_pair

def significant_eqtl_writer(genes_vars_pairs, output_dir):
	"""
	f_out = open("significant_%i_eqtl_info.csv"%len(genes_vars_pairs), "w")
	f_out.write("gene" + "," + "variant" + "," + "chrom" + "," + "pos" + "\n")
	for genes_vars_pair in genes_vars_pairs:
		chrom = genes_vars_pairs[genes_vars_pair].split('_')[0]
		pos = genes_vars_pairs[genes_vars_pair].split('_')[1]
		f_out.write(genes_vars_pair + "," + genes_vars_pairs[genes_vars_pair] + "," + 
			str(chrom) + "," + str(pos) + "\n")
	"""
	f_out_eqtl = open(os.path.join(output_dir, "significant_%i_eqtl_info.txt"%len(genes_vars_pairs)), "w")
	f_out_gene = open(os.path.join(output_dir, "gene_list.txt"), "w")
	for genes_vars_pair in genes_vars_pairs:
		chrom = genes_vars_pairs[genes_vars_pair].split('_')[0]
		pos = genes_vars_pairs[genes_vars_pair].split('_')[1]
		f_out_eqtl.write(str(chrom) + "\t" + str(pos) + "\n") #for the input of bcftools
		f_out_gene.write(str(genes_vars_pair) + "\t"+ str(chrom) + "\t" + str(pos) + '\n')
	f_out_eqtl.close()


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Get eligible cis-eqtl from fastQTL report')
	parser.add_argument('egenes_report', help='tissueName.egenes.txt.gz')
	parser.add_argument('-q', '--logq', type=float, action='store', default=None, help='threshold mode: minus log q-value')
	parser.add_argument('-n', type=int, action='store', default=None, help='threshold mode: top n significant eqtls')
	parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	args = parser.parse_args()

	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)
	egenes_report = args.egenes_report
	nModeThresh = args.n
	qModeThresh = args.logq
	output_dir = args.output_dir
	if egenes_report.endswith("txt.gz"):
		f_In = gzip.open(egenes_report, 'r')
	else:
		f_In = open(egenes_report, 'r')
	#find significant gene-eqtl pairs
	if nModeThresh:
		top_n_genes_vars_pair = find_top_n_eqtl(f_In, nModeThresh)
		f_In.close()
		significant_eqtl_writer(top_n_genes_vars_pair, output_dir)
		print("Top %i mode:" % nModeSwitch)
		print("Top %i genes-variance pairs found" % nModeThresh)
	if qModeThresh:
		significant_genes_vars_pair = find_eqtl_according_to_qval(f_In, qModeThresh)
		f_In.close()
		significant_eqtl_writer(significant_genes_vars_pair, output_dir)
		print("Q-value 10^-%i mode:" % qModeThresh)
		print("Top %i genes-variance pairs found" % len(significant_genes_vars_pair))
