#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 2024
"""
import pandas as pd
import scanpy as sc
import anndata
import os
filtered_dir = "/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/GSE201446_RAW/filtered_files"
combined_data = None
# Loop through all files in the filtered directory
for filtered_file in os.listdir(filtered_dir):
    # Load each filtered matrix using Scanpy
    adata= sc.read(filtered_dir+"/"+filtered_file)

    if combined_data is None:
        combined_data = adata
    else:
        combined_data = anndata.concat([combined_data,adata],join="outer")

# Save the combined data
combined_data.write("/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/Anno/Axolotl_combined_data.h5ad")