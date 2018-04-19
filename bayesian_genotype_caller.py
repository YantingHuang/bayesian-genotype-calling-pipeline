#! /home/yhua295/.linuxbrew/bin/python3

import os
import pandas as pd
from scipy.stats import norm
import argparse
import gzip
import numpy as np
import pymc3 as pm 
import rnaseqnorm
from scipy.stats import gaussian_kde

def read_gt_file(gt_file, ind_id_list, gene_dict):
	"""
	read genotype data into pandas dataframe
	return a dataframe of variants genotypes with index as gene names
	"""
	col_names = ["chrom", "pos"] + ind_id_list[1:]
	rep_gt_df = pd.DataFrame() #for those snps corresponding to different egenes
	gt_df = pd.read_csv(gt_file, sep=" ", header=None)
	gt_df = gt_df.drop(2, axis=1)
	gt_df.columns = col_names
	rows_to_drop = []
	unique_index = []
	# delete unnecessary rows
	for index, row in gt_df.iterrows():
		if (row["chrom"], row["pos"]) not in gene_dict:
			rows_to_drop += [index]
			continue
		if (row["chrom"], row["pos"]) not in unique_index:
			unique_index += [(row["chrom"], row["pos"])]
		else:
			rows_to_drop += [index]
	gt_df = gt_df.drop(rows_to_drop, axis="index")
	# add duplicate snps
	index_list = []
	dup_snp_index_list = []
	for index, row in gt_df.iterrows():
		posInfo = (row["chrom"], row["pos"])
		index_list += [gene_dict[posInfo][0]]
		if len(gene_dict[posInfo]) != 1:
			rep_gt_df = rep_gt_df.append([row]*(len(gene_dict[posInfo])-1), ignore_index=True)
			dup_snp_index_list += gene_dict[posInfo][1:]
	gt_df.index = index_list
	rep_gt_df.index = dup_snp_index_list
	gt_df = pd.concat([gt_df, rep_gt_df], axis=0)
	return gt_df

def fill_missing_val(x, sim_dict):
	"""
	x: rows in gt_df (a pandas series)
	"""
	if (len(sim_dict.get(x.name)) != 0):
		x.loc[x=='./.'] = sim_dict.get(x.name)
		return x
	else:
		return x

def generate_simulation(x):
	sim = []
	if x.loc['./.'] != 0:
		for i in range(int(x.loc['./.'])):
			tmp_0_count = 0
			for j in range(2):
				if np.random.uniform(0, 1) <= x.loc[0]:
					tmp_0_count += 1
			if tmp_0_count == 0:
				sim.append('1/1')
			elif tmp_0_count == 1:
				sim.append('0/1')
			elif tmp_0_count == 2:
				sim.append('0/0')
		return (x.name, sim)
	else:
		return (x.name, sim)

def calcualte_freq(gt_df):
	"""
	calculate the frequency of each genotype (AA, Aa, aa) for each eqtl
	4 possible genotype for each eqtl (misssing value ./. included)

	Missing values are filled according to allele frequency
	"""
	gt_df = gt_df.iloc[:, 2:]
	counts_res = gt_df.apply(lambda x: x.value_counts(), axis=1) #shape: nGenes * 4
	counts_res.fillna(0, inplace=True)
	allele_0_freq = counts_res.loc[:, '0/0']*2 + counts_res.loc[:, '0/1']
	allele_1_freq = counts_res.loc[:, '1/1']*2 + counts_res.loc[:, '0/1']
	allele_0_prob = allele_0_freq/(allele_0_freq+allele_1_freq)
	allele_1_prob = allele_1_freq/(allele_0_freq+allele_1_freq)
	allele_prob = pd.concat([allele_0_prob, allele_1_prob, counts_res], axis=1)
	simulation_gt_for_missingVal = [generate_simulation(row) for index, row in allele_prob.iterrows()]
	simulation_gt_for_missingVal_dict = dict(simulation_gt_for_missingVal)
	filled_gt_df = pd.DataFrame()
	for index, row in gt_df.iterrows():
		filled_gt_df = filled_gt_df.append(fill_missing_val(row, simulation_gt_for_missingVal_dict))
	#filled_gt_df = gt_df.apply(fill_missing_val, sim_dict=simulation_gt_for_missingVal_dict, axis=1)
	updated_counts_res = filled_gt_df.apply(lambda x: x.value_counts(), axis=1)
	updated_counts_res.fillna(0, inplace=True)
	priors = updated_counts_res.div(updated_counts_res.sum(axis=1), axis=0)
	return priors

def calculate_mean(train_gt, train_expression, genotype):
	return train_gt.apply(lambda x: train_expression.loc[x.name][x==genotype].mean(), axis=1).to_frame(name=genotype)

def calculate_std(train_gt, train_expression, genotype):
	return train_gt.apply(lambda x: train_expression.loc[x.name][x==genotype].std(), axis=1).to_frame(name=genotype)

def estimate_statistic_mcmc(data):
	"""
	To be done
	"""
	with pm.Model() as model:
		mu = pm.Normal('mu', mu=0, sd=5)
		std = pm.Normal('std', mu=1, sd=3)
		obs = pm.Normal('obs', mu=mu, sd=std, observed=data)
	with model:
		trace = pm.sample(1000)
	#mu = pm.summary(trace)
	return pm.summary(trace)

def calculate_mean_std_mcmc(train_gt, train_expression, genotype):
	"""
	To be done
	"""
	return train_gt.apply(lambda x: estimate_statistic_mcmc(train_expression.loc[x.name][x==genotype]), axis=1).to_frame(name=genotype)

def get_summary_stats(train_gt, train_expression):
	"""
	calculate the sample mean and std for each genotype of each gene
	"""
	train_gt = train_gt.iloc[:, 2:]
	train_expression.columns = ['-'.join(name.split('-')[:2]) for name in train_expression.columns.tolist()]
	genotypes = ['0/0', '0/1', '1/1']
	sample_mean_df = pd.concat([calculate_mean(train_gt, train_expression, '0/0'), 
		calculate_mean(train_gt, train_expression, '0/1'), calculate_mean(train_gt, train_expression, '1/1')], axis=1)
	sample_std_df = pd.concat([calculate_std(train_gt, train_expression, '0/0'), 
		calculate_std(train_gt, train_expression, '0/1'), calculate_std(train_gt, train_expression, '1/1')], axis=1)
	return sample_mean_df.copy(), sample_std_df.copy()

def prepare_expression(counts_df, tpm_df, mode='tmm'):
    """
    This part and Normalization part is adapted from gtex official pipeline

    Genes are thresholded based on the following expression rules:
      TPM >= tpm_threshold in >= sample_frac_threshold*samples
      read counts >= count_threshold in sample_frac_threshold*samples
    
    vcf_lookup: lookup table mapping sample IDs to VCF IDs
    
    Between-sample normalization modes:
      tmm: TMM from edgeR
      qn:  quantile normalization
    """
    ns = tpm_df.shape[1]
    # apply normalization
    if mode.lower()=='tmm':
        tmm_counts_df = rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True)
        norm_df = rnaseqnorm.inverse_normal_transform(tmm_counts_df)
    elif mode.lower()=='qn':
        qn_df = rnaseqnorm.normalize_quantiles(tpm_df)
        norm_df = rnaseqnorm.inverse_normal_transform(qn_df)
    else:
        raise ValueError('Unsupported mode {}'.format(mode))

    return norm_df

def get_platform_info(cov_list):
	df = pd.read_csv(cov_list, sep='\t', index_col=0)
	platform_info = df.loc["platform", :]
	return dict(platform_info)

def get_specific_platform_only(platform_info, colnames, platform_id, full=True):
	res=[]
	if full:
		for colname in colnames:
			tmpcolname = '-'.join(colname.split('-')[:2])
			if platform_info.get(tmpcolname) == platform_id:
				res += [colname]
	else:
		for colname in colnames:
			if platform_info.get(colname) == platform_id:
				res += [colname]
	return res

def eval(label_series, pred_series):
	eval_tmp_df = pd.DataFrame({"preds": pred_series, "labels": label_series})
	eval_tmp_df = eval_tmp_df[eval_tmp_df.loc[:, 'labels']!=3] #only evaluate the genes with knwon genotypes in the test set
	cmp_res = (eval_tmp_df["labels"] == eval_tmp_df["preds"]).values.astype(int)
	nComparableGenes =  cmp_res.shape[0]
	acc = cmp_res.sum()/nComparableGenes
	return acc

def predict_eval(sample_mean_df, sample_std_df, priors_df, normlized_test_df, test_gt):
	genotypes = ['0/0', '0/1', '1/1']
	"""
	for row, index in normlized_test_df.iterrows():
		geneName = row.name
		likelihood = norm.pdf(row.values, sample_mean_df.loc[geneName, '0/0'], sample_std_df.loc[geneName, '0/0']) 
	"""
	#calculate likelihood based on normal distributed assumption
	#likelihood_0_0 = normlized_test_df.apply(lambda x: norm.pdf(x, sample_mean_df.loc[x.name, '0/0'], sample_std_df.loc[x.name, '0/0']), axis=1)
	#likelihood_0_1 = normlized_test_df.apply(lambda x: norm.pdf(x, sample_mean_df.loc[x.name, '0/1'], sample_std_df.loc[x.name, '0/1']), axis=1)
	#likelihood_1_1 = normlized_test_df.apply(lambda x: norm.pdf(x, sample_mean_df.loc[x.name, '1/1'], sample_std_df.loc[x.name, '1/1']), axis=1)
	likelihood_0_0 = normlized_test_df.apply(lambda x: norm.pdf(x, sample_mean_df.loc[x.name, '0/0'], sample_std_df.loc[x.name]), axis=1)
	likelihood_0_1 = normlized_test_df.apply(lambda x: norm.pdf(x, sample_mean_df.loc[x.name, '0/1'], sample_std_df.loc[x.name]), axis=1)
	likelihood_1_1 = normlized_test_df.apply(lambda x: norm.pdf(x, sample_mean_df.loc[x.name, '1/1'], sample_std_df.loc[x.name]), axis=1)
	#calculate the posterior for 3 genotypes
	posterior_0_0 = likelihood_0_0.mul(priors_df.loc[:,'0/0'], axis=0) #df multiply series
	posterior_0_1 = likelihood_0_1.mul(priors_df.loc[:,'0/1'], axis=0)
	posterior_1_1 = likelihood_1_1.mul(priors_df.loc[:,'1/1'], axis=0)
	posteriors = np.dstack((posterior_0_0.values, posterior_0_1.values, posterior_1_1.values)) #shape: #genes * #samples * 3

	predictions = np.argmax(posteriors, axis=2) #shape: #genes * #samples

	predictions = pd.DataFrame(predictions, index=posterior_0_0.index, columns=posterior_0_0.columns)
	predictions.columns = ['-'.join(label.split('-')[:2]) for label in predictions.columns.tolist()]
	test_gt_labels = test_gt.iloc[:, 2:]
	test_gt_labels[test_gt_labels=="0/0"] = 0
	test_gt_labels[test_gt_labels=="0/1"] = 1
	test_gt_labels[test_gt_labels=="1/1"] = 2
	test_gt_labels[test_gt_labels=="./."] = 3
	import pdb;pdb.set_trace()
	#test result by individuals
	acc_list_ind = []
	individual_list = []
	for individual in predictions.columns.tolist():
		label_series = test_gt_labels.loc[:, individual]
		pred_series = predictions.loc[:, individual]
		acc = eval(label_series, pred_series)
		acc_list_ind += [acc]
		individual_list += [individual]
	acc_dict_ind = dict(zip(individual_list, acc_list_ind))
	with open(os.path.join(args.output_dir, "acc_byIndividuals_testSet.csv"), "w") as f:
		for ind in acc_dict_ind:
			f.write(ind+","+str(acc_dict_ind[ind])+"\n")
	f.close()
	#test result by variants
	acc_list_var = []
	var_list = []
	for var in predictions.index.tolist():
		label_series = test_gt_labels.loc[var, :]
		pred_series = predictions.loc[var, :]
		acc = eval(label_series, pred_series)
		acc_list_var += [acc]
		var_list += [var]
	acc_dict_var = dict(zip(var_list, acc_list_var))	
	with open(os.path.join(args.output_dir, "acc_byVariants_testSet.csv"), "w") as f:
		for var in acc_dict_var:
			f.write(var+","+str(acc_dict_var[var])+"\n")
	f.close()
	return acc_dict_ind, acc_dict_var

if __name__=='__main__':
	np.random.seed(10) #ensure reproducibility

	parser = argparse.ArgumentParser(description='Use Bayesian method to directly call the genotype')
	parser.add_argument('expression_data', help='rpkm/tpm')
	parser.add_argument('read_counts', help='read counts')
	parser.add_argument('train_genotype_data', help='output of genotype_finder')
	parser.add_argument('test_genotype_data', help='output of genotype_finder')
	parser.add_argument('train_id_list', help='train_id_list')
	parser.add_argument('test_id_list', help='test_id_list')
	parser.add_argument('train_sample_list', help='train_sample_list')
	parser.add_argument('test_sample_list', help='test_sample_list')
	parser.add_argument('gene_list', help='gene_list')
	parser.add_argument('cov_list', help='covariate')
	parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
	parser.add_argument('-r', '--rpkm', action='store_true', help='expression data type')
	parser.add_argument('--norm', default='tmm', help='tmm/qn')
	args = parser.parse_args()

	cov_list = args.cov_list
	platform_info = get_platform_info(cov_list)

	#read relevant files into memory
	with open(args.train_id_list) as f:
		train_ids = f.read().strip().split("\n")
		print("%i training samples" % len(train_ids))
	f.close()
	with open(args.test_id_list) as f:
		test_ids = f.read().strip().split("\n")
		print("%i test samples" % len(test_ids))
	f.close()
	with open(args.train_sample_list) as f:
		train_samples = f.read().strip().split("\n")
	f.close()
	with open(args.test_sample_list) as f:
		test_samples = f.read().strip().split("\n")
	f.close()
	gene_dict = {}
	gene_list = []
	with open(args.gene_list) as f:
		for line in f:
			line = line.strip("\n").split("\t")
			posInfo = (int(line[1]), int(line[2]))
			if posInfo in gene_dict:
				gene_dict[posInfo] += [line[0]]
			else:
				gene_dict[posInfo] = [line[0]]
			gene_list.append(line[0])
	f.close()
	#decide what type of expression data to use
	if args.rpkm:
		expression_data_type = "rpkm"
	else:
		expression_data_type = "tpm"
	print("Expression data type: %s"%expression_data_type)
	#Load genotype data
	print("Start reading genotype data...")
	#train_ids = ['tmp'] + get_specific_platform_only(platform_info, train_ids, platform_id=1, full=False)
	#test_ids = ['tmp'] + get_specific_platform_only(platform_info, test_ids, platform_id=1, full=False)
	train_ids = ['tmp'] + train_ids
	test_ids = ['tmp'] + test_ids

	train_gt = read_gt_file(args.train_genotype_data, train_ids, gene_dict)
	test_gt = read_gt_file(args.test_genotype_data, test_ids, gene_dict)
	train_ids = ['chrom', 'pos'] + get_specific_platform_only(platform_info, train_gt.columns, platform_id=1, full=False)
	test_ids = ['chrom', 'pos'] + get_specific_platform_only(platform_info, test_gt.columns, platform_id=1, full=False)
	train_gt = train_gt.loc[:, train_ids]
	test_gt = test_gt.loc[:, test_ids]
	print("Done")
	print("%i genes invloved" % train_gt.shape[0])
	print("Total number of unique variants: %i" % len(gene_dict))

	#calculate the frequency of 0/0 0/1 1/1
	#This step is important since we need to decide the priors of three genotypes for each gene here
	print("Calculate priors...")
	priors_df = calcualte_freq(train_gt)
	print("Done")

	#Load count data
	print("Start reading counts data...")
	train_counts_file = os.path.join(args.output_dir, "train_counts")
	test_counts_file = os.path.join(args.output_dir, "test_counts")
	if (os.path.isfile(train_counts_file) and os.path.isfile(test_counts_file)):
		train_counts = pd.read_csv(train_counts_file, sep='\t', index_col=0)
		test_counts = pd.read_csv(test_counts_file, sep='\t', index_col=0)
	else:
		train_samples = ['Name'] + train_samples
		test_samples = ['Name'] + test_samples
		train_counts = pd.read_csv(args.read_counts, sep='\t', usecols=train_samples, index_col=0, skiprows=2)
		test_counts = pd.read_csv(args.expression_data, sep='\t', usecols=test_samples, index_col=0, skiprows=2)
		#
		train_cols_full = get_specific_platform_only(platform_info, train_counts.columns, platform_id=1, full=True)
		test_cols_full = get_specific_platform_only(platform_info, test_counts.columns, platform_id=1, full=True)

		egenes = train_counts.index.intersection(gene_list)
		train_counts = train_counts.loc[egenes, train_cols_full]
		test_counts = test_counts.loc[egenes, test_cols_full]
		import pdb; pdb.set_trace()
		train_counts.to_csv(train_counts_file, sep='\t')
		test_counts.to_csv(test_counts_file, sep='\t')
	print("Done")
	#read tpm gct file into pandas dataframe
	#use preprocessed expression data if exists to save loading time
	#load original tpm/rpkm data and normalize o/w
	norm_method = args.norm
	pre_train_exp = os.path.join(args.output_dir, "train_expression_%s_%s"%(expression_data_type, norm_method))
	pre_test_exp  = os.path.join(args.output_dir, "test_expression_%s_%s"%(expression_data_type, norm_method))
	if (os.path.isfile(pre_train_exp) and os.path.isfile(pre_test_exp)):
		print("Start reading pre-existing expression data...")
		normlized_train_df = pd.read_csv(pre_train_exp, sep='\t', index_col=0)
		normlized_test_df = pd.read_csv(pre_test_exp, sep='\t', index_col=0)
		print("Done")
	else:
		train_samples = ['Name'] + train_samples
		test_samples = ['Name'] + test_samples
		print("Start reading expression data...")
		#train_expression_df = pd.read_csv(args.expression_data, sep='\t', usecols=train_samples, index_col=0, skiprows=2)
		train_expression_df = pd.read_csv(args.expression_data, sep='\t', usecols=train_cols_full, index_col=0, skiprows=2)
		egenes = train_expression_df.index.intersection(gene_list)
		train_expression_df = train_expression_df.loc[egenes, :]
		#test_expression_df = pd.read_csv(args.expression_data, sep='\t', usecols=test_samples, index_col=0, skiprows=2)
		test_expression_df = pd.read_csv(args.expression_data, sep='\t', usecols=test_cols_full, index_col=0, skiprows=2)
		test_expression_df = test_expression_df.loc[egenes, :]
		print("Done")
		print("Start normalizing data for each gene using %s..." % norm_method)
		train_mean = train_expression_df.mean(axis=1)
		train_std = train_expression_df.std(axis=1)
		normlized_train_df = prepare_expression(train_counts, train_expression_df, mode=norm_method)
		normlized_test_df = prepare_expression(test_counts, test_expression_df, mode=norm_method)
		normlized_train_df.to_csv(pre_train_exp, sep='\t')
		normlized_test_df.to_csv(pre_test_exp, sep='\t')
		"""
		if normalize:
			normlized_train_df = train_df.sub(train_mean, axis=0).div(train_std, axis=0)
			normlized_test_df = test_df.sub(train_mean, axis=0).div(train_std, axis=0)
			print("Start writing normalized expression data for future use")
			normlized_train_df.to_csv(pre_train_exp, sep='\t')
			normlized_test_df.to_csv(pre_test_exp, sep='\t')
		else:
			normlized_train_df = train_df
			normlized_test_df = test_df
			print("Start writing normalized expression data for future use")
			normlized_train_df.to_csv(pre_train_exp, sep='\t')
			normlized_test_df.to_csv(pre_test_exp, sep='\t')
		"""
		#normlized_test_df = test_df.sub(test_df.mean(axis=1), axis=0).div(test_df.mean(axis=1), axis=0)
		print("Done")
	
	#log-transform
	#normlized_train_df = np.sqrt(normlized_train_df.sub(1.01*normlized_train_df.min(axis=1), axis="index"))
	#normlized_test_df = np.sqrt(normlized_test_df.sub(1.01*normlized_test_df.min(axis=1), axis="index"))
	"""
	# mcmc to estimate the mean and variance for each gene
	gName = 'ENSG00000148484.13'
	train_gt = train_gt.iloc[:, 2:]
	normlized_train_df.columns = ['-'.join(name.split('-')[:2]) for name in normlized_train_df.columns.tolist()]
	summary = estimate_statistic_mcmc(normlized_train_df.loc[gName][train_gt.loc[gName]=='0/1'])
	import pdb; pdb.set_trace()
	"""
	

	print("Start calculating summary statistic for training data...")
	sample_mean_df, sample_std_df = get_summary_stats(train_gt, normlized_train_df)
	sample_mean_df.fillna(0, inplace=True)
	#sample_std_df.fillna(1000, inplace=True)

	sample_std_df.fillna(0, inplace=True)
	sample_std_df = sample_std_df.mul(priors_df).sum(axis=1)
	print("Done")
	#prediction and evaluation
	print("Start predicting and evaluating on the test set...")
	acc_dict_ind, acc_dict_var = predict_eval(sample_mean_df, sample_std_df, priors_df, normlized_test_df, test_gt)
	print("Done")
