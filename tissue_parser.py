import os
import pandas as pd
import numpy as np
import argparse
import gzip

def train_test_split(samples_list, train_proportion, output_dir):
	f_In = open(samples_list, 'r')
	nTrain = 0
	nTest = 0
	train_id = "train_individual_ids.txt"
	test_id = "test_individual_ids.txt"
	train_sample = "train_sample_ids.txt"
	test_sample = "test_sample_ids.txt"
	f_out_train_sample = open(os.path.join(output_dir, train_sample), "w")
	f_out_test_sample = open(os.path.join(output_dir, test_sample), "w")
	f_out_train_id = open(os.path.join(output_dir, train_id), "w")
	f_out_test_id = open(os.path.join(output_dir, test_id), "w")
	for line in f_In:
		ind_id = '-'.join(line.split('-')[0:2])
		if(np.random.uniform(0, 1) <= train_proportion):
			f_out_train_sample.write(line)
			f_out_train_id.write(ind_id + '\n')
			nTrain += 1
		else:
			f_out_test_sample.write(line)
			f_out_test_id.write(ind_id + '\n')
			nTest += 1
	print("nTrain: %i" % nTrain)
	print("nTest: %i" % nTest)
	f_In.close()
	f_out_train_sample.close()
	f_out_test_sample.close()
	f_out_train_id.close()
	f_out_test_id.close()
	return nTrain, nTest

def num_peer_factors(num_of_samples):
	"""
	Decide the number of peer factors needed according to number of samples used
	"""
	if num_of_samples < 150:
		return 15
	elif num_of_samples < 250:
		return 30
	else:
		return 35

if __name__ == "__main__":
	np.random.seed(10) #ensure reproducibility
	parser = argparse.ArgumentParser(description='Generate covariates for FastQTL cis-eqtl mapping')
	parser.add_argument('sample_attr', help='phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz')
	parser.add_argument('sample_gct', help='sample gct file used for getting the header')
	parser.add_argument('cov_list', help='covariates from the official site')
	parser.add_argument('tissue_name', default=None, help='')
	parser.add_argument('pc_prefix', help='pc_prefix')
	parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
	args = parser.parse_args()

	if not os.path.exists(args.output_dir):
	    os.makedirs(args.output_dir)

	gtex_ids = []
	"""
	'Whole Blood', 'Adipose - Subcutaneous', 'Muscle - Skeletal', 'Artery - Tibial', 
	'Artery - Coronary', 'Heart - Atrial Appendage', 'Adipose - Visceral (Omentum)', 'Ovary', 'Uterus', 
	'Vagina', 'Breast - Mammary Tissue', 'Skin - Not Sun Exposed (Suprapubic)', 'Minor Salivary Gland', 'Brain - Cortex',
	 'Adrenal Gland', 'Thyroid', 'Lung', 'Spleen', 'Pancreas', 'Esophagus - Muscularis', 'Esophagus - Mucosa', 
	 'Esophagus - Gastroesophageal Junction', 'Stomach', 'Colon - Sigmoid', 'Small Intestine - Terminal Ileum', 
	 'Colon - Transverse', 'Prostate', 'Testis', 'Skin - Sun Exposed (Lower leg)', 'Nerve - Tibial', 'Heart - Left Ventricle', 
	 'Pituitary', 'Brain - Cerebellum', 'Cells - Transformed fibroblasts', 'Artery - Aorta', 'Cells - EBV-transformed lymphocytes', 
	 'Liver', 'Kidney - Cortex', 'Brain - Hippocampus', 'Brain - Substantia nigra', 'Brain - Anterior cingulate cortex (BA24)', 
	 'Brain - Frontal Cortex (BA9)', 'Brain - Cerebellar Hemisphere', 'Brain - Caudate (basal ganglia)', 
	 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Hypothalamus', 
	 'Brain - Spinal cord (cervical c-1)', 'Brain - Amygdala', 'Fallopian Tube', 'Bladder', 'Cervix - Ectocervix', 
	 'Cervix - Endocervix', 'Cells - Leukemia cell line (CML)'
	 """
	with gzip.open(args.sample_attr, "rt") as f_In:
		for line in f_In:
			if not (line.startswith("#") or line.startswith("dbGaP_Sample_ID") ):
				line = line.strip("\n").split("\t")
				if (len(line) > 1):
					tissue = line[14]
					SMAFRZE = line[28] #samples included/excluded from the eqtl analysis
					gtex_id = line[1]
					if (tissue == args.tissue_name and SMAFRZE != 'EXCLUDE'):
						gtex_ids.append(gtex_id)
	f_In.close()

	df = pd.read_csv(args.sample_gct, sep='\t', skiprows=2, index_col=0)
	rna_seqed_sample_list = df.columns

	gtex_ids = rna_seqed_sample_list.intersection(gtex_ids)

	#use samples used by official gtex project
	df_getex = pd.read_csv(args.cov_list, index_col=0, sep='\t')
	official_sample_list = df_getex.columns
	offcial_individual_ids = [sample.split('-')[1] for sample in official_sample_list]
	preserved_gtex_idx = []
	original_col_names = []
	for gtex_id in gtex_ids:
		if gtex_id.split('-')[1] in offcial_individual_ids:
			original_col_names.append("GTEX"+"-"+gtex_id.split('-')[1])
			preserved_gtex_idx.append(gtex_id)
	pc_covs_df = df_getex[original_col_names]
	#pc_covs_df.columns = preserved_gtex_idx #if you want sample ids instead of patient ids
	pc_covs_df = pc_covs_df.iloc[[0,1,2,-2,-1],:]
	#save the genotype-based pcs/sex/ sequencing platforms info
	pc_covs_df.to_csv(os.path.join(args.output_dir, 
		args.pc_prefix+"_pcs.txt"), sep='\t')

	with open(os.path.join(args.output_dir, 
		"samples_of_"+"_".join(args.tissue_name.split(" "))+".txt"), "w") as f_out:
		for gtex_id in preserved_gtex_idx:
			f_out.write(gtex_id+"\n")
	f_out.close()
	samples_list = os.path.join(args.output_dir, 
		"samples_of_"+"_".join(args.tissue_name.split(" "))+".txt")
	nTrain, nTest = train_test_split(samples_list, train_proportion=0.8, output_dir=args.output_dir)

	#summary statistic report
	print(args.tissue_name)
	print(str(num_peer_factors(nTrain)) + " peer factors needed")

