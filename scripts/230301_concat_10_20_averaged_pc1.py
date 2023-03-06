#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import numpy as np
import os
import glob


# In[2]:


cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']


# In[5]:


cpg_type = 'opensea'
data_dir = '/data/project/3dith/data'
matrix_type = 'bdm'


# ---

# # avg10, all 22 autosomes

# In[13]:


flag = 'avg10'
for cohort in cohorts:
    print(f"===\n{cohort}")
    
    cohort_pc1_dir = os.path.join(data_dir, cpg_type, cohort, 'pc1')
    all_pc1_files = glob.glob(os.path.join(cohort_pc1_dir, f'*pc1.standardized-concat-{flag}*'))
    samplenames = []
    for i in range(len(all_pc1_files)):
        f = all_pc1_files[i]
        samplenames.append(os.path.basename(f).split('.pc1')[0])
        samplename = os.path.basename(f).split('.pc1')[0]
        if i == 0:
            all_data = np.load(f, allow_pickle = True)[matrix_type]
        else:
            all_data = np.vstack((all_data, np.load(f, allow_pickle = True)[matrix_type]))
    all_df = pd.DataFrame(all_data, index = samplenames)
    all_df_fname = os.path.join(data_dir, cpg_type, cohort, f'all_samples_pc1_concat_{matrix_type}.s_per_chrom_{flag}.csv')
    all_df.to_csv(all_df_fname)
    print(all_df_fname)


# ---
# # avg20, all 22 autosomes

# In[14]:


flag = 'avg20'
for cohort in cohorts:
    print(f"===\n{cohort}")
    
    cohort_pc1_dir = os.path.join(data_dir, cpg_type, cohort, 'pc1')
    all_pc1_files = glob.glob(os.path.join(cohort_pc1_dir, f'*pc1.standardized-concat-{flag}*'))
    samplenames = []
    for i in range(len(all_pc1_files)):
        f = all_pc1_files[i]
        samplenames.append(os.path.basename(f).split('.pc1')[0])
        samplename = os.path.basename(f).split('.pc1')[0]
        if i == 0:
            all_data = np.load(f, allow_pickle = True)[matrix_type]
        else:
            all_data = np.vstack((all_data, np.load(f, allow_pickle = True)[matrix_type]))
    all_df = pd.DataFrame(all_data, index = samplenames)
    all_df_fname = os.path.join(data_dir, cpg_type, cohort, f'all_samples_pc1_concat_{matrix_type}.s_per_chrom_{flag}.csv')
    all_df.to_csv(all_df_fname)
    print(all_df_fname)


# ---
# # First 8 autosomes

# In[15]:


flag = '8'
for cohort in cohorts:
    print(f"===\n{cohort}")
    
    cohort_pc1_dir = os.path.join(data_dir, cpg_type, cohort, 'pc1')
    all_pc1_files = glob.glob(os.path.join(cohort_pc1_dir, f'*pc1.standardized-concat-{flag}*'))
    samplenames = []
    for i in range(len(all_pc1_files)):
        f = all_pc1_files[i]
        samplenames.append(os.path.basename(f).split('.pc1')[0])
        samplename = os.path.basename(f).split('.pc1')[0]
        if i == 0:
            all_data = np.load(f, allow_pickle = True)[matrix_type]
        else:
            all_data = np.vstack((all_data, np.load(f, allow_pickle = True)[matrix_type]))
    all_df = pd.DataFrame(all_data, index = samplenames)
    all_df_fname = os.path.join(data_dir, cpg_type, cohort, f'all_samples_pc1_concat_{matrix_type}.s_per_chrom_{flag}.csv')
    all_df.to_csv(all_df_fname)
    print(all_df_fname)


# In[ ]:




