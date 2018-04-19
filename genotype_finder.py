#! /home/yhua295/.linuxbrew/bin/python3
import argparse
import numpy as np
import sys,os

parser = argparse.ArgumentParser(description='Find genotypes of significant eqtls from the vcf file')
parser.add_argument('vcf_file', help='gtex_vcf')
parser.add_argument('eqtl_info', help='eqtl_info')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')

def find_genotypes(eqtl_info, train_id, test_id, vcf, output_dir):
	"""
	find genotypes of given eqtls using bcftools
	bcftools manuals can be found at <https://samtools.github.io/bcftools/bcftools.html>
	"""
	train_id = os.path.join(output_dir, train_id)
	test_id = os.path.join(output_dir, test_id)
	train_set_genotype = os.path.join(output_dir, "train_set_genotype")
	test_set_genotype = os.path.join(output_dir, "test_set_genotype")
	os.system("bcftools query -R %s -S %s -f '%%CHROM %%POS GT:[ %%GT]\n' %s>%s" % (eqtl_info, train_id, vcf, train_set_genotype))
	os.system("bcftools query -R %s -S %s -f '%%CHROM %%POS GT:[ %%GT]\n' %s>%s" % (eqtl_info, test_id, vcf, test_set_genotype))

if __name__ == "__main__":
	args = parser.parse_args()
	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)
	output_dir = args.output_dir
	vcf = args.vcf_file
	eqtl_info = args.eqtl_info
	train_id = "train_individual_ids.txt"
	test_id = "test_individual_ids.txt"
	find_genotypes(eqtl_info, train_id, test_id, vcf, output_dir)

