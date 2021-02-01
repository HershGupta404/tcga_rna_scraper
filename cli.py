#!/usr/bin/env python3


import argparse 
import os 
import sys
import scraper
from shutil import copy2
import subprocess
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

log = logging.getLogger('CLI')

def parse_cli():
	parser = argparse.ArgumentParser(description = "Parse TCGA metadata files and either download files or get summary statistics")
	parser.add_argument("histology", help="Name of histology")
	parser.add_argument("root_folder", help="Parent folder where to build the backend file system.")
	parser.add_argument("json_file", help="Location of json_file. This should be an absolute path. File will be duplicated to histology folder.")
	parser.add_argument("config", help="Location of config file to support code.")
	parser.add_argument("--optional", action = "store_true")
	return parser.parse_args()	


def read_config(config_file):
	pass

def full_load(histology,root_location,json_file,include_others, config_options=None):
	original_location = os.getcwd()
	os.chdir(root_location)
	# Build the file structure as shown in the diagram
	os.mkdir(histology)
	os.chdir(histology)
	os.mkdir("scraped_data")
	os.mkdir("deseq_input")
	os.mkdir("deseq_output")

	# Move to scraping part and begin to scrape the data
	json_file_name = json_file.split("/")[-1]
	copy2(json_file,"scraped_data/"+json_file_name)
	os.chdir("scraped_data")
	file_table = scraper.parse_json(json_file_name, include_all=include_others)

	#Split into normal, tumor, and possibly other files and download the files
	normal_files = file_table[file_table["type"]=="Normal"].index
	tumor_files = file_table[file_table["type"]=="Tumor"].index
	scraper.download_files(normal_files, "normal_tissue")
	normal_matrix=scraper.combine_raw_rna("normal_tissue")
	scraper.download_files(tumor_files, "tumor_tissue")
	tumor_matrix=scraper.combine_raw_rna("tumor_tissue")
	if include_others:
		other_files = file_table[file_table["type"]=="Other"]
		scraper.download_files(other_files, "other_tissue")
		other_matrix = scraper.combine_raw_rna("other_tissue")
		# Create the deseq matrices
		deseq_matrix, sample_info = scraper.create_deseq_files([tumor_matrix,normal_matrix,other_matrix],["Tumor","Normal","Other"])
	else:
		deseq_matrix, sample_info = scraper.create_deseq_files([tumor_matrix,normal_matrix])

	# Deseq input save off 

	log.info("Running DESeq")
	os.chdir("../deseq_input")
	deseq_matrix.to_csv("raw_counts.csv")
	sample_info.to_csv("sample_table.csv")

	# Run DESeq script
	deseq_directory = os.getcwd()
	os.chdir(original_location)
	subprocess.run("Rscript deseq_script.R "+deseq_directory, shell=True)
	log.info("Files saved off successfully")



if __name__ == '__main__':
	arg_dicts = parse_cli()
	full_load(arg_dicts.histology,arg_dicts.root_folder, arg_dicts.json_file, arg_dicts.optional)


