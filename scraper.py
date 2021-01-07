#!/usr/bin/env python3

import pandas as pd
import numpy as np
from glob import glob
import subprocess
import json
import requests
import os 





def parse_json(case_list, include_all =  False, sample_file_save = True, location = None):
	"""
	Pareses a TCGA RNA-seq case list containg possible HTSEQ counts files for download. 
	include_all -- Include all cases, not just solid tissue tumor primaries and normals. Any other samples will be grouped in one category (default False)
	sample_file_save -- Saves a list of the final files and case/sample information in location (default True)
	location -- location for saved sample file (default Containing folder of case list)
	"""
	with open(case_list) as f:
		file_list_json = json.load(f)

	data_list={}
	for entity in file_list_json:
		case_code = entity["associated_entities"][0]["entity_submitter_id"].split("-")[3][:2]
		if case_code == "01":
			data_list[entity["file_id"]]= ["Tumor",entity["associated_entities"][0]["case_id"]]
		elif case_code == "11":
			data_list[entity["file_id"]]= ["Normal",entity["associated_entities"][0]["case_id"]]
		elif include_all:
			data_list[entity["file_id"]]= ["Other",entity["associated_entities"][0]["case_id"]]
	file_entity_df = pd.DataFrame.from_dict(data_list,orient = "index",columns=["type","sample_name"])
	file_entity_df.index.name="file_name"
	if location == None:
		location = "/".join(case_list.split("/")[:-1])
		location += "/"
	file_entity_df.to_csv(location+"file_sample_list.csv")
	return file_entity_df

#TODO: Find a better way than subprocess and curl to deal with this piece of work
def download_files(tumor_file_list, normal_file_list, location, other_file_list = None):
	os.chdir(location)
	os.mkdir("normal_files")
	os.chdir("normal_files")
	for i,file in enumerate(normal_file_list):
		print("Total files:", len(normal_file_list))
		if i%50 ==0:
			print("On file:",i)
    	subprocess.run("curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/"+file+"'",shell=True)
	os.chdir(location)
	os.mkdir("tumor_files")
	os.chdir("tumor_files")
	for i,file in enumerate(tumor_file_list):
		print("Total files:", len(tumor_file_list))
		if i%50 ==0:
			print("On file:",i)
		subprocess.run("curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/"+file+"'",shell=True)
	if other_file_list != None:
		for i,file in enumerate(other_file_list):
			print("Total files:", len(other_file_list))
			if i%50 ==0:
				print("On file:",i)
			subprocess.run("curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/"+file+"'",shell=True)



