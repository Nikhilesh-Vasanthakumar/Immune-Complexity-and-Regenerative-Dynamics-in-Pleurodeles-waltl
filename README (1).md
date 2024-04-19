
# Unveiling Immune Complexity and Regenerative Dynamics in Pleurodeles waltl: A Cross-Species Single-Cell Analysis

This repository contains the code and analysis for studying the immune compartments of *Pleurodeles waltl* by performing cross-species comparisions with humans and *Xenopus laevis* .Through this comparative approach, our aim is to uncover the cellular and immune mechanisms that facilitate regeneration and to understand how these processes differ across species with varying regenerative capacities. The analysis pipeline includes differential gene expression analysis, gene ontology enrichment analysis (GOEA) to identify significant expression changes and pathways activated with a specific focus on the adaptive immunity including B and T cells.





## Data Utilized

The raw gene expression count matrices of the spleen samples from Pleurodeles waltl were obtained from a previous research study conducted evaluating the genetic demultiplexing of the scRNA sequencing data for this dataset. The dataset consisted of a total of five biological replicates, which were divided into two count matrices. The first one containing three biological replicates [dataset1](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12186) and the second dataset containing two biological replicates [dataset2](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12182)
## Folder Structure

The repository contains all the necessary python scripts and files to reproduce the analysis. The folder structure is organized as follows:

```python
- data/              (contains the dataset files after ambient RNA filtering)
- scripts/           (contains python scripts for data processing, analysis, and visualization)
- results/           (stores the generated plots and analysis results)
- README.md          (the main readme file)
- Requirements.md     (script to install required python packages)
```
## Installation

Clone the repository to your local machine using the following command:

```python
git clone https://github.com/Nikhilesh-Vasanthakumar/BINP51
```

The following packages were utilised for the analysis and can be found in the requirements file.
The comparative analysis was performed using Python v3.8 and the following primary packages were used Scanpy (v1.9.6) , SAMap (v1.0.15), Cellbender (v0.3.0), scvi tools (v1.0.4) , NumPy (v1.23.5), Anndata (v0.9.2), Harmony (v0.0.9), scikit-learn (v1.3,1) , ESM-2 (v1.0.3), AGAT (v1.3.0), GOATOOLS (v1.1.5).
## Workflow
The analysis workflow involves several steps.Here is an overview of the workflow:

1.  **Ambient RNA filtering using Cellbender**
To filter out the datasets for the high level of ambient RNA we use the Cellbender tool. The script for running the tool is in the scripts folder.

2.  **Data Loading and Preprocessing**: 
We used the scanpy package to load the scRNA-seq datasets. Eggnog gene annotations were incorporated from the Trinotate tool. During preprocessing, we filtered the datasets to remove ambient RNA and doublets, and we excluded cells of low quality. Subsequently, the gene expression data was normalized and scaled to facilitate accurate downstream analyses and comparisons. 

3.  **Differnatial Expression Analysis and Gene Ontology Enrichment Analysis**:
Using the acquired gene annotations we Identified the immune clusters using the expression of marker genes. Further we subsetted the potential B and T subclusters to identify the differentially expressed genes and the pathways that these genes are involved in. We performed DEA and GOEA on the identified subclusters of B and T cells.

4.  **SAMap Cross-species Analysis**:
To utilise SAMap we first have to perform pairwise blast (The script is in the scripts folder ) of the protein sequences of the genes in the gtf files. After which we perform the SAMap analysis.After which we calculate the cell type alignments, commonly enriched genes and gene triangles.


## Data Preprocessing

After running the cellBender and pairwise blast you can start the preprocessing step.
```python
import scanpy as sc
import pandas as pd
import anndata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import requests
import sklearn.decomposition #import TruncatedSVD
import scipy 
import sparse, io
import scvi


#Load in the raw first spleen dataset
adata_pl_1_raw=sc.read_10x_h5("/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/Anno/matrix_filtered_clear/run2/output_filtered.h5") #Load in your cellbender output for dataset1

#Load in the raw second spleen dataset
adata_pl_2_raw= sc.read_10x_h5('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/Anno/matrix_filtered_clear/second_spleen/output_filtered.h5') #Load in your cellbender output for dataset2

#Load in the souporcell annotated spleen dataset
souporcell_1=pd.read_table('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nick/Pleuro/outs/20230704_spleen_R1_soup/20230603_24h_soup_res/clusters.tsv')
souporcell_2=pd.read_table('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nick/Pleuro/outs/20230704_spleen_R2_soup/20230603_24h_soup_res/clusters.tsv')
#Replace the variables load in the triannotate file
mapper=pd.read_table("/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/aPlwal.pri.V2.genome.annots.tsv")
mapper_dict = mapper.set_index('#gene_id')['EggNM.Preferred_name'].to_dict()
status_1_mapper=dict(zip(souporcell_1['barcode'],souporcell_1['status']))
assignment_1_mapper=dict(zip(souporcell_1['barcode'],souporcell_1['assignment']))
status_2_mapper=dict(zip(souporcell_2['barcode'],souporcell_2['status']))
assignment_2_mapper=dict(zip(souporcell_2['barcode'],souporcell_2['assignment']))

#Add in the souporcell annotation to the first and second spleen dataset
adata_pl_1_raw.obs['status']='NA'
adata_pl_1_raw.obs['assignment']='NA'
adata_pl_1_raw.obs['batch']='1'
adata_pl_2_raw.obs['status']='NA'
adata_pl_2_raw.obs['assignment']='NA'
adata_pl_2_raw.obs['batch']='2'
for i in adata_pl_1_raw.obs.index:
    if i in status_1_mapper.keys():
        adata_pl_1_raw.obs.loc[i,'status']=status_1_mapper[i]
        adata_pl_1_raw.obs.loc[i,'assignment']=assignment_1_mapper[i]
for i in adata_pl_2_raw.obs.index:
    if i in status_2_mapper.keys():
        adata_pl_2_raw.obs.loc[i,'status']=status_2_mapper[i]
        adata_pl_2_raw.obs.loc[i,'assignment']=assignment_2_mapper[i]

#Filter out the cells
adata_pl_1_raw.var_names = [mapper_dict.get(x, x) if mapper_dict.get(x, x) != '.' else x for x in adata_pl_1_raw.var_names]
adata_pl_2_raw.var_names = [mapper_dict.get(x, x) if mapper_dict.get(x, x) != '.' else x for x in adata_pl_2_raw.var_names]
#Preprocess the data
mt_gene_patterns = ['COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'NU1M', 'NU2M', 'NU3M', 'NU4M', 'NU4LM', 'NU5M', 'NU6M', 'CYB']
mt_gene_pattern = '|'.join(mt_gene_patterns)
for adata in [adata_pl_1_raw, adata_pl_2_raw]:
    sc.pp.filter_cells(adata, min_genes=400)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.match(mt_gene_pattern)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    adata.var_names_make_unique()

scvi.model.SCVI.setup_anndata(adata_pl_1_raw)
vae = scvi.model.SCVI(adata_pl_1_raw)
vae.train()
solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()
df_1 = solo.predict()
df_1['prediction'] = solo.predict(soft=False)

scvi.model.SCVI.setup_anndata(adata_pl_2_raw)
vae = scvi.model.SCVI(adata_pl_2_raw)
vae.train()
solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()
df_2 = solo.predict()
df_2['prediction'] = solo.predict(soft=False)

#Remove Doublets from the second spleen dataset
adata_pl_2_raw_dob = adata_pl_2_raw[(adata_pl_2_raw.obs['status'] == 'doublet') ]
cellid_soup=adata_pl_2_raw_dob.obs.index
cellid_scvi=df_2[df_2['prediction']=='doublet'].index
common_elements = set(cellid_scvi).intersection(cellid_soup)
cellid_scvi = set(cellid_scvi)
cellid_soup = set(cellid_soup)
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
venn2([cellid_scvi, cellid_soup], ('cellid_scvi', 'cellid_soup'))
doublets= cellid_scvi.union(cellid_soup)
doublets=list(doublets)
doublets_in_adata_pl_2_raw = list(set(doublets).intersection(adata_pl_2_raw.obs_names))
adata_pl_2_raw.obs.loc[doublets_in_adata_pl_2_raw, 'status'] = 'doublet'
adata_pl_2_raw=adata_pl_2_raw[(adata_pl_2_raw.obs['status'] == 'singlet') ]

adata_pl_1_raw_dob = adata_pl_1_raw[(adata_pl_1_raw.obs['status'] == 'doublet') ]
cellid_soup=adata_pl_1_raw_dob.obs.index
cellid_scvi=df_1[df_1['prediction']=='doublet'].index
common_elements = set(cellid_scvi).intersection(cellid_soup)
cellid_scvi = set(cellid_scvi)
cellid_soup = set(cellid_soup)
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
venn2([cellid_scvi, cellid_soup], ('cellid_scvi', 'cellid_soup'))
doublets= cellid_scvi.union(cellid_soup)
doublets=list(doublets)

doublets_in_adata_pl_1_raw= list(set(doublets).intersection(adata_pl_1_raw.obs_names))
adata_pl_1_raw.obs.loc[doublets_in_adata_pl_1_raw, 'status'] = 'doublet'
adata_pl_1_raw=adata_pl_1_raw[adata_pl_1_raw.obs['status'] == 'singlet']

#Rename the animals in the second spleen dataset to be 3 and 4
adata_pl_2_raw.obs['assignment']=adata_pl_2_raw.obs['assignment'].replace({'1':'3','2':'4'})


#Integrate the two spleen datasets
sc.pp.normalize_total(adata_pl_1_raw, target_sum=1e4)
sc.pp.normalize_total(adata_pl_2_raw, target_sum=1e4)
sc.pp.log1p(adata_pl_1_raw)
sc.pp.log1p(adata_pl_2_raw)
sc.pp.highly_variable_genes(adata_pl_1_raw, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(adata_pl_2_raw, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata_pl_1_raw, max_value=10)
sc.pp.scale(adata_pl_2_raw, max_value=10)
sc.tl.pca(adata_pl_1_raw, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_pl_1_raw, log=True)
sc.tl.pca(adata_pl_2_raw, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_pl_2_raw, log=True)

#Clustering the data
sc.pp.neighbors(adata_pl_1_raw, n_neighbors=30, n_pcs=40)
sc.tl.umap(adata_pl_1_raw)
sc.tl.leiden(adata_pl_1_raw,resolution=0.5)
sc.pl.umap(adata_pl_1_raw, color=['leiden'])
sc.pp.neighbors(adata_pl_2_raw, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata_pl_2_raw)
sc.tl.leiden(adata_pl_2_raw,resolution=0.5)
sc.pl.umap(adata_pl_2_raw, color=['leiden'])

#Integrate the two spleen datasets
adata_pl_1_raw.var_names_make_unique()
adata_pl_1_raw.obs["dataset"]="1"
adata_pl_2_raw.var_names_make_unique()
adata_pl_2_raw.obs["dataset"]="2"
var_names= adata_pl_1_raw.var_names.intersection(adata_pl_2_raw.var_names)
adata_pl_1_raw=adata_pl_1_raw[:,var_names]
adata_pl_2_raw=adata_pl_2_raw[:,var_names]
spleen_merged=adata_pl_1_raw.concatenate(adata_pl_2_raw)

#Batch correct using Harmony
sc.external.pp.harmony_integrate(spleen_merged, ['assignment','dataset'])

#Re cluster the integrated spleen dataset
sc.pp.neighbors(spleen_merged, n_neighbors=30, n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(spleen_merged)
sc.tl.leiden(spleen_merged,resolution=0.5)
sc.pl.umap(spleen_merged, color=['leiden'])
```
## Batch Correction and Adaptive Subclustering

After preprocessing we visualise the cluster expression and identify differentially expressed genes and subcluster them based on the visulaised expression

```python
import anndata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.rcParams.update({'font.size': 12})
%config InlineBackend.figure_format = 'retina'

#Load in the integrated spleen dataset
adata_pl_raw=anndata.read_h5ad('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Inter_data/spleen_merged_raw.h5ad')

#Subset for B cells

adata_b_cells=adata_pl_raw[adata_pl_raw.obs['leiden_spleen'].isin(['0', '5', '7','1','8'])]

#Preprocess the data
mt_gene_patterns = ['COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'NU1M', 'NU2M', 'NU3M', 'NU4M', 'NU4LM', 'NU5M', 'NU6M', 'CYB']
mt_gene_pattern = '|'.join(mt_gene_patterns)
sc.pp.filter_cells(adata_b_cells, min_genes=400)
sc.pp.filter_genes(adata_b_cells, min_cells=3)
adata_b_cells.var['mt'] = adata_b_cells.var_names.str.match(mt_gene_pattern)
sc.pp.calculate_qc_metrics(adata_b_cells, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_b_cells = adata_b_cells[adata_b_cells.obs.pct_counts_mt < 10, :]
adata_b_cells.var_names_make_unique()

#Normalize the data
sc.pp.normalize_total(adata_b_cells, target_sum=1e4)
sc.pp.log1p(adata_b_cells)
sc.pp.scale(adata_b_cells, max_value=10)
sc.pp.pca(adata_b_cells, n_comps=50,svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_b_cells, n_pcs=50)

#Clustering

#UMAP
sc.pp.neighbors(adata_b_cells, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_b_cells)
sc.pl.umap(adata_b_cells,color=['leiden_spleen'])

#DE analysis
sc.tl.rank_genes_groups(adata_b_cells, 'leiden_spleen', method='wilcoxon')
sc.pl.rank_genes_groups(adata_b_cells, n_genes=25, sharey=False)
sc.tl.leiden(adata_b_cells, resolution=0.5)
sc.pl.umap(adata_b_cells,color=['leiden','assignment'])

#Batch correction using Harmony
sc.external.pp.harmony_integrate(adata_b_cells, key='assignment', max_iter_harmony=100)
#Re Preprocess the data
sc.pp.neighbors(adata_b_cells, n_neighbors=10, n_pcs=20,use_rep='X_pca_harmony')
sc.tl.umap(adata_b_cells)
sc.tl.leiden(adata_b_cells, resolution=0.5)
sc.pl.umap(adata_b_cells,color=['leiden','leiden_spleen','assignment'])

#Remove cells in cluster 6,11 #This can change based on the version of scanpy
adata_b_cells=adata_b_cells[adata_b_cells.obs['leiden'].isin(['0','1','2','3','5','7','4','8','9','10','12','13'])]
#Recluster
sc.pp.neighbors(adata_b_cells, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_b_cells)
sc.tl.leiden(adata_b_cells, resolution=0.5)
#Batch correction using Harmony
sc.external.pp.harmony_integrate(adata_b_cells, key='assignment', max_iter_harmony=100)
#Re Preprocess the data
sc.pp.neighbors(adata_b_cells, n_neighbors=10, n_pcs=20,use_rep='X_pca_harmony')
sc.tl.umap(adata_b_cells)
sc.tl.leiden(adata_b_cells, resolution=0.5)
sc.pl.umap(adata_b_cells,color=['leiden','leiden_spleen','assignment'])
#Plot the markers
sc.pl.umap(adata_b_cells,color=['PTPRC','CD79A','CD79B','CD3E','leiden'],legend_loc='on data')
#DE analysis
sc.tl.rank_genes_groups(adata_b_cells, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_b_cells, n_genes=25, sharey=False)

#Subset for T cells
adata_t_cells=adata_pl_raw[adata_pl_raw.obs['leiden_spleen'].isin(['2','14','12'])]
#Preprocess the data
mt_gene_patterns = ['COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'NU1M', 'NU2M', 'NU3M', 'NU4M', 'NU4LM', 'NU5M', 'NU6M', 'CYB']
mt_gene_pattern = '|'.join(mt_gene_patterns)
sc.pp.filter_cells(adata_t_cells, min_genes=400)
sc.pp.filter_genes(adata_t_cells, min_cells=3)
adata_t_cells.var['mt'] = adata_t_cells.var_names.str.match(mt_gene_pattern)
sc.pp.calculate_qc_metrics(adata_t_cells, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_t_cells = adata_t_cells[adata_t_cells.obs.pct_counts_mt < 10, :]
adata_t_cells.var_names_make_unique()
#Normalize the data
sc.pp.normalize_total(adata_t_cells, target_sum=1e4)
sc.pp.log1p(adata_t_cells)
sc.pp.scale(adata_t_cells, max_value=10)
sc.pp.pca(adata_t_cells, n_comps=50,svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_t_cells, n_pcs=50)
#UMAP
sc.pp.neighbors(adata_t_cells, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_t_cells)
sc.pl.umap(adata_t_cells,color=['leiden_spleen'])
#DE analysis
sc.tl.rank_genes_groups(adata_t_cells, 'leiden_spleen', method='wilcoxon')
sc.pl.rank_genes_groups(adata_t_cells, n_genes=25, sharey=False)
#Clustering again
sc.tl.leiden(adata_t_cells, resolution=0.5)
sc.pl.umap(adata_t_cells,color=['leiden','leiden_spleen','assignment'])

#Batch correction using Harmony
sc.external.pp.harmony_integrate(adata_t_cells, key='assignment', max_iter_harmony=100)
#Re Preprocess the data
sc.pp.neighbors(adata_t_cells, n_neighbors=10, n_pcs=20,use_rep='X_pca_harmony')
sc.tl.umap(adata_t_cells)
sc.tl.leiden(adata_t_cells, resolution=0.5)
sc.pl.umap(adata_t_cells,color=['leiden','leiden_spleen','assignment'],legend_loc='on data')

#Further subset for T cells based on T cell marker expression
adata_t_cells=adata_t_cells[adata_t_cells.obs['leiden'].isin(['0','3','6','8','9','10','7'])]

#Recluster after batch correction
sc.external.pp.harmony_integrate(adata_t_cells, key='assignment', max_iter_harmony=100)
sc.pp.neighbors(adata_t_cells, n_neighbors=10, n_pcs=20,use_rep='X_pca_harmony')
sc.tl.umap(adata_t_cells)
sc.tl.leiden(adata_t_cells, resolution=0.5)

#DEG
sc.tl.rank_genes_groups(adata_t_cells, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_t_cells, n_genes=25, sharey=False)

adata_b_cells.write_h5ad('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Final_Data/spleen_b_cells.h5ad')
adata_t_cells.write_h5ad('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Final_Data/spleen_t_cells.h5ad')

```
## Gene Ontology Enrichment Analysis

Using the subclustered datasets we perform DE analysis and GO enrichment analysis.
```python
import anndata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


#Load in the integrated spleen dataset
adata_pl_raw=anndata.read_h5ad('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Inter_data/spleen_merged_raw.h5ad')
adata_b_cells=anndata.read_h5ad("/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Final_Data/spleen_b_cells.h5ad")
adata_t_cells=anndata.read_h5ad("/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Final_Data/spleen_t_cells.h5ad")

#GO Analysis
import wget
import os
from goatools import obo_parser
go_db_url=' http://purl.obolibrary.org/obo/go/go-basic.obo'
data_folder="/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/data"# Replace this with your data folder

if(not os.path.isfile(data_folder+'/go-basic.obo')):
    go_obo = wget.download(go_db_url, data_folder+'/go-basic.obo')
else:
    go_obo = data_folder+'/go-basic.obo'

#Create a id2gos format file using Trianotate
triannotate=pd.read_csv('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/aPlwal.pri.V2.genome.annots.tsv',sep='\t')

#Creating a GODAG for Pleurodeles

from collections import defaultdict

eggnog_mapper = defaultdict(list)
preferred_names = defaultdict(int)

for i, row in triannotate.iterrows():
    gene_id = row['#gene_id']
    preferred_name = row['EggNM.Preferred_name']
    go_terms = row['EggNM.GOs']

    if go_terms != '.':
        # Determine the key based on the preferred name and gene ID
        if preferred_name != '.':
            key = preferred_name
            # Update the preferred_names dictionary
            preferred_names[preferred_name] += 1
        else:
            key = gene_id
            # Check if the gene_id has appeared before
            if gene_id in preferred_names:
                # Use the preferred name from the first appearance
                key = f"{gene_id}-{preferred_names[gene_id]}"
                preferred_names[gene_id] += 1

        # Check if the key already exists in the dictionary
        if key in eggnog_mapper:
            # Append a suffix starting from -1
            suffix = 1
            new_key = f"{key}-{suffix}"
            while new_key in eggnog_mapper:
                suffix += 1
                new_key = f"{key}-{suffix}"
            key = new_key

        eggnog_mapper[key].append(go_terms)

# Convert defaultdict to a regular dictionary
eggnog_mapper = dict(eggnog_mapper)

#Remove values with only '.' in the values
eggnog_mapper_filtered = {}
gene_list = []
for key in eggnog_mapper.keys():
    if eggnog_mapper[key] == '.':
        continue
    else:
        eggnog_mapper_filtered[key] = eggnog_mapper[key]
        gene_list.append(key)

#Write the mapper to create a id2gos file
with open('/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Inter_data/id2gos.txt', 'w') as f: #Change this to where you want the dataset
    for key, go_terms in eggnog_mapper.items():
        # Join the list of GO IDs with semi-colons
        go_ids_str = go_terms[0].replace(',', ';') if len(go_terms) == 1 else ';'.join(go_terms)
        # Write the gene name, a tab, and the list of GO IDs
        f.write(f"{key}\t{go_ids_str}\n")
#Collecting the top 200 differentially expressed gene

#For T cells
sc.tl.rank_genes_groups(adata_t_cells,groupby='leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_t_cells, n_genes=25, sharey=False)
cluster0=sc.get.rank_genes_groups_df(adata_t_cells,group='0',key='rank_genes_groups')
cluster1=sc.get.rank_genes_groups_df(adata_t_cells,group='1',key='rank_genes_groups')
cluster2=sc.get.rank_genes_groups_df(adata_t_cells,group='2',key='rank_genes_groups')
cluster3=sc.get.rank_genes_groups_df(adata_t_cells,group='3',key='rank_genes_groups')
cluster4=sc.get.rank_genes_groups_df(adata_t_cells,group='4',key='rank_genes_groups')
cluster5=sc.get.rank_genes_groups_df(adata_t_cells,group='5',key='rank_genes_groups')
cluster6=sc.get.rank_genes_groups_df(adata_t_cells,group='6',key='rank_genes_groups')
cluster7=sc.get.rank_genes_groups_df(adata_t_cells,group='7',key='rank_genes_groups')
cluster8=sc.get.rank_genes_groups_df(adata_t_cells,group='8',key='rank_genes_groups')
cluster0=cluster0['names'][:200]
cluster1=cluster1['names'][:200]
cluster2=cluster2['names'][:200]
cluster3=cluster3['names'][:200]
cluster4=cluster4['names'][:200]
cluster5=cluster5['names'][:200]
cluster6=cluster6['names'][:200]
cluster7=cluster7['names'][:200]
cluster8=cluster8['names'][:200]

#For B cells
sc.tl.rank_genes_groups(adata_b_cells,groupby='leiden',method='wilcoxon')
sc.pl.rank_genes_groups(adata_b_cells,n_genes=25,sharey=False)
cluster0_b=sc.get.rank_genes_groups_df(adata_b_cells,group='0',key='rank_genes_groups')
cluster1_b=sc.get.rank_genes_groups_df(adata_b_cells,group='1',key='rank_genes_groups')
cluster2_b=sc.get.rank_genes_groups_df(adata_b_cells,group='2',key='rank_genes_groups')
cluster3_b=sc.get.rank_genes_groups_df(adata_b_cells,group='3',key='rank_genes_groups')
cluster4_b=sc.get.rank_genes_groups_df(adata_b_cells,group='4',key='rank_genes_groups')
cluster5_b=sc.get.rank_genes_groups_df(adata_b_cells,group='5',key='rank_genes_groups')
cluster6_b=sc.get.rank_genes_groups_df(adata_b_cells,group='6',key='rank_genes_groups')
cluster7_b=sc.get.rank_genes_groups_df(adata_b_cells,group='7',key='rank_genes_groups')
cluster8_b=sc.get.rank_genes_groups_df(adata_b_cells,group='8',key='rank_genes_groups')
cluster9_b=sc.get.rank_genes_groups_df(adata_b_cells,group='9',key='rank_genes_groups')
cluster10_b=sc.get.rank_genes_groups_df(adata_b_cells,group='10',key='rank_genes_groups')
cluster11_b=sc.get.rank_genes_groups_df(adata_b_cells,group='11',key='rank_genes_groups')
cluster12_b=sc.get.rank_genes_groups_df(adata_b_cells,group='12',key='rank_genes_groups')
cluster0_b=cluster0_b['names'][:200]
cluster1_b=cluster1_b['names'][:200]
cluster2_b=cluster2_b['names'][:200]
cluster3_b=cluster3_b['names'][:200]
cluster4_b=cluster4_b['names'][:200]
cluster5_b=cluster5_b['names'][:200]
cluster6_b=cluster6_b['names'][:200]
cluster7_b=cluster7_b['names'][:200]
cluster8_b=cluster8_b['names'][:200]
cluster9_b=cluster9_b['names'][:200]
cluster10_b=cluster10_b['names'][:200]
cluster11_b=cluster11_b['names'][:200]
cluster12_b=cluster12_b['names'][:200]

#Gene Ontology Enrichment Analysis

from goatools.go_enrichment import GOEnrichmentStudy
from goatools.anno.idtogos_reader import IdToGosReader
import scipy
import fisher
from goatools.base import download_go_basic_obo
def go_enrichment(cluster):
    
    id2gos_file_path = '/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Inter_data/id2gos.txt'

    # Path to the Gene Ontology OBO file
    obo_file='/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/data/go-basic.obo'
    obo_file =GODag(obo_file)
# Create an instance of IdToGosReader and read the gene to GO annotations
    id2gos_reader = IdToGosReader(id2gos_file_path)
    id2gos_dict = id2gos_reader.get_id2gos()

    # Background set of genes (all genes in your study)
    background_genes = set(adata_pl_raw.var_names)

    # Your set of genes of interest
    genes_of_interest = list(cluster)  # Replace with your actual gene IDs

    # Create an instance of GOEnrichmentStudy
    go_enrichment = GOEnrichmentStudy(
        background_genes, # List of all genes in my study
        id2gos_dict,# List of Pleuro genes
        obo_file,
        propagate_counts=False,# Set to True if you want to propagate counts up the GO hierarchy
        alpha=0.01,  # Set your desired significance level
        methods=['fdr_bh']  # Set the multiple testing correction method
    )

    # Run the GO enrichment analysis
    go_results = go_enrichment.run_study(genes_of_interest)
    return go_results

#Creating a filtered table for the GO analysis
def create_go_table(go_results):
    goea_results_sig = [r for r in go_results if r.p_fdr_bh < 0.05 and len(r.study_items) > 1]

    # Create DataFrame from the significant results
    GO = pd.DataFrame([{
        'GO': res.GO,
        'term': res.goterm.name,
        'class': res.goterm.namespace,
        'p': res.p_uncorrected,
        'p_corr': res.p_fdr_bh,
        'n_genes': res.ratio_in_study[0],
        'n_study': res.ratio_in_study[1],
        'n_go': len(res.goterm.get_all_parents()) + 1, # Assuming you want to count the term itself plus all parents
        'study_genes': list(res.study_items)
    } for res in goea_results_sig])
    
    return GO


#Plot the results of GOEA
def go_plot(dataframe):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    
    df = dataframe.copy()
    # Convert p-values to a negative log scale to amplify differences
    df['-log10(p_corr)'] = -np.log10(df['p_corr'])

    # Sort the DataFrame based on the significance level of GO terms
    df_sorted = df.sort_values('-log10(p_corr)', ascending=False)

    # Select the top N significant terms for plotting
    N = 20
    df_top = df_sorted.head(N)
    
    # Normalize the number of genes for color mapping across the full range (blue to red)
    gene_counts_norm = (df_top['n_genes'] - df_top['n_genes'].min()) / (df_top['n_genes'].max() - df_top['n_genes'].min())
    colors = plt.cm.coolwarm(gene_counts_norm)

    # Create the bar plot
    plt.figure(figsize=(10, 8))
    ax = sns.barplot(x='-log10(p_corr)', y='term', data=df_top,
                     palette=colors,
                     edgecolor='black')

    # Create color bar for the number of genes
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(vmin=df_top['n_genes'].min(), vmax=df_top['n_genes'].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Number of Genes', rotation=270, labelpad=15)

    # Set the plot labels and title
    plt.xlabel('-log10(Adjusted P-Value)')
    plt.ylabel('GO Term')
    plt.title('Cluster12 Significant GO Terms by Adjusted P-Value')

    plt.show()

#Use these functions for the entire process.
go_results = go_enrichment(cluster12_b)
GO_table = create_go_table(go_results)
go_plot(GO_table)
```
## SAMap Cross-Species Analysis
SAMap requires the raw dataset to stitch in the counts matrix based on the blast result

```python
import anndata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import scanpy as sc
import sklearn.decomposition #import TruncatedSVD for imputing values
import scipy 


#Samap analysis Packages
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles,
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)

matplotlib.rcParams.update({'font.size': 12})
%config InlineBackend.figure_format = 'retina'

org1='/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Inter_data/spleen_merged_raw.h5ad'
org2='/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/Anno/frog_laevis_raw_with_leiden_without_L&S.h5ad'
#org3 ='/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/Anno/Tabula_immune.h5ad'

filenames = {'pw':org1,'xl':org2} #Add org3 if you want to compare the three species
sm = SAMAP(
        filenames,
        f_maps='/home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Inter_data/blast_mappings/maps/',
        save_processed=True, #if False, do not save the processed results to `*_pr.h5ad`
        keys = {'pw':'leiden','xl':'leiden','hs':'cell_ontology_class'},
   )
sm.run()
samap = sm.samap

#Plotting the stiched data
adata_combined=samap.adata
sc.pl.umap(adata_combined,color='species')

#Ct and Mapping Table contains cluster alignment and score
keys_ct = {'pw':'leiden','xl':'leiden','hs':'cell_ontology_class'} #Find mapping scores between all the three species cluster add or remove as needed
Ct,MappingTable_ct = get_mapping_scores(sm,keys_ct,n_top = 100)


#Find Gene pairs in all the clusters in the species selected above
gpf=GenePairFinder(sm,keys_ct)
gene_pairs = gpf.find_all(align_thr=0.2)

#CAN ONLY BE RUN FOR THREE SPECIES
keys={'pw':'leiden','xl':'leiden','hs':'cell_type'}
result = CellTypeTriangles(sm,keys,align_thr=0.15)
result
```
## Conclusion

In summary, our study sheds light on the diverse immune cell populations present in *Pleurodeles waltl* and their significance in the context of its remarkable regenerative abilities. Through advanced techniques such as SAMap analysis, we have identified novel immune cell subsets and elucidated their roles in the immune response and possible tissue regeneration processes.Overall, our study paves the way for future investigations into the therapeutic applications of regenerative biology, rooted in a comprehensive understanding of immune cell dynamics and evolutionary adaptations.