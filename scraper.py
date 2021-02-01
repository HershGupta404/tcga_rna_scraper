#!/usr/bin/env python3

import pandas as pd
import numpy as np
from glob import glob
import subprocess
import json
import requests
import os 
from collections import Counter
import logging



logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

log = logging.getLogger('Scraper')

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
    if sample_file_save:
        if location == None:
            location = "/".join(case_list.split("/")[:-1])
            location += "/"
            if location == "/":
                location = "./"
        file_entity_df.to_csv(location+"file_sample_list.csv")
    log.info("Total counts of generated file types")
    print(Counter(file_entity_df['type']))
    return file_entity_df

#TODO: Find a better way than subprocess and curl to deal with this piece of work
def download_files(file_list, folder_title, location = "./"):
    """
    Downloads passed list of files in a certain location. Default structure is a normal file folder, a tumor folder and a possible other folder
    tumor_file_list -- All file names for tumor files
    normal_file_list -- All file names for normal files 
    location -- location to download files #TODO: Assign a default location? 
    other_files
    """
    og_location=os.getcwd()
    os.chdir(location)
    location = os.getcwd()
    try:
        os.mkdir(folder_title)
        os.chdir(folder_title)
        log.info(msg=("Total " + folder_title + ":" +str(len(file_list))))
        for i,file in enumerate(file_list):
            if i%50 ==0:
                log.info(msg=("On file: "+str(i)))
            subprocess.run("curl --silent --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/"+file+"'",shell=True)
        subprocess.run("gunzip *.counts.gz", shell=True)
    except FileExistsError: 
        log.warning("Folder already exists- using subroutine to find new files")
    os.chdir(og_location)


def combine_raw_rna(loc, ensmbl_hugo_map ="ensmbl_hugo_map.json"):
    """
    Combines raw RNA files at a a loc into one large file. Requires a mapping file to convert from ensmbl to hugo symbols
    """

    #TODO: change this to a config files option   
    with open("/Users/hershgupta/work/ensmbl_hugo_map.json") as f:
        gene_mapping = json.load(f)
    og_location = os.getcwd()
    os.chdir(loc)
    all_samples = pd.DataFrame()
    file_order = glob("*.counts")
    log.info("Combining RNA files at "+loc)
    for i,file in enumerate(file_order):
        if i%50==0:
            log.info("On file:" +str(i))
        curr_sample=pd.read_csv(file,sep=None,index_col=0,header=None,engine="python")
        all_samples=pd.concat([all_samples,curr_sample],axis=1)
    # Strip away the file thing and continue to use file identifer
    all_samples.columns=[i.split(".")[0] for i in file_order]
    #Last 5 rows are not useful
    all_samples=all_samples.iloc[:-5,:]
    #Rename index to hugo symbols unless some not present
    all_samples.index = [i.split(".")[0] for i in all_samples.index]
    new_index = []
    for gene in all_samples.index:
        if gene in gene_mapping:
            new_index.append(gene_mapping[gene])
        else:
            new_index.append(gene)
    assert len(new_index) == len(all_samples.index)
    all_samples.index = new_index
    # Remove all the non-zero samples
    all_samples_nz = all_samples.loc[~(all_samples==0).all(axis=1)]
    os.chdir(og_location)
    return all_samples_nz


def create_deseq_files(count_list, sample_class = ["Tumor","Normal"]):
    """
    Prepare DESEQ files
    """
    assert len(count_list) == len(sample_class)
    full_raw_counts = pd.DataFrame()
    for counts in count_list:
        full_raw_counts = counts.join(full_raw_counts)
    full_raw_counts=full_raw_counts.fillna(0)
    # Remove duplicated genes
    log.warning("The following genes are duplicated, if any are important for downstream analysis manually remove duplicates")
    duplicated_genes = full_raw_counts[full_raw_counts.index.duplicated(keep="first")].index
    log.warning(set(duplicated_genes))
    full_raw_counts = full_raw_counts[~full_raw_counts.index.duplicated(keep="first")]
    #Create sample list
    sample_info =pd.DataFrame(index=full_raw_counts.columns)
    sample_list = []
    for i in sample_info.index:
        for j,data in enumerate(count_list):
            if i in data.columns:
                sample_list.append(sample_class[j])
                continue
    log.info("Total samples:"+str(len(sample_list)))
    sample_info['Type']=sample_list
    return full_raw_counts, sample_info










