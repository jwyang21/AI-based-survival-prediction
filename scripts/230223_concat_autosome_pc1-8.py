#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import glob


# ## global variable

# In[2]:


cohort_fname = '/data/project/3dith/data/cohorts.manifest'
cohorts = pd.read_csv(cohort_fname).cohort.values
tcga_cohorts = []
for c in cohorts:
    if 'TCGA' in c:
        tcga_cohorts.append(c)


# In[3]:


data_dir = '/data/project/3dith/data/'
#cpg_type = 'opensea'


# In[4]:


all_cpg_types = ['opensea','island', 'shelf_shore']


# In[5]:


samplenames_fname = '/data/project/3dith/data/samplenames.npz'


# In[6]:


AUT_LIST = [f'chr{i}' for i in np.arange(1, 23)]
SMALL_AUT_LIST = [f'chr{i}' for i in np.arange(1, 9)]


# In[7]:


def standardize(v):
    return (v-np.mean(v))/np.std(v)


# ---

# ## concat v1: Concat all pc1 vectors of 22 autosomes, per sample, per cohort (raw)

# In[8]:


for cpg_type in all_cpg_types:
    for cohort in tcga_cohorts:
        print(f'===\n{cohort}')
        cohort_data_dir = os.path.join(data_dir, cpg_type, cohort, 'pc1')
        cohort_samples = np.load(samplenames_fname,allow_pickle=True)[cohort]

        all_pc1_files = glob.glob(os.path.join(cohort_data_dir, f'*.npz'))
        bdm_pc1_files = []
        iebdm_pc1_files = []
        for f in all_pc1_files:
            if 'concat' not in f:
                if f.endswith('_inv_exp.npz'):
                    iebdm_pc1_files.append(f)
                else:
                    bdm_pc1_files.append(f)
        assert len(bdm_pc1_files) == len(iebdm_pc1_files) and len(bdm_pc1_files) == len(cohort_samples)

        for sample in cohort_samples:
            bdm_pc1_fname = os.path.join(cohort_data_dir, f'{sample}.npz')
            iebdm_pc1_fname = os.path.join(cohort_data_dir, f'{sample}_inv_exp.npz')
            bdm_pc1 = np.load(bdm_pc1_fname, allow_pickle = True)
            iebdm_pc1 = np.load(iebdm_pc1_fname, allow_pickle = True)

            globals()[f'{sample}'] = {}
            total_n_bins = 0
            
            #for chrom in AUT_LIST:
            for chrom in SMALL_AUT_LIST:
                k = f'{chrom}_pc1'
                
                #if chrom == AUT_LIST[0]:
                if chrom == SMALL_AUT_LIST[0]:
                    all_bdm_pc1 = bdm_pc1[k]
                    all_iebdm_pc1 = iebdm_pc1[k]
                else:
                    all_bdm_pc1 = np.concatenate((all_bdm_pc1, bdm_pc1[k]), axis = None)
                    all_iebdm_pc1 = np.concatenate((all_iebdm_pc1, iebdm_pc1[k]), axis = None)
                assert len(bdm_pc1[k]) == len(iebdm_pc1[k])
                total_n_bins += len(bdm_pc1[k])
            assert len(all_bdm_pc1) == total_n_bins and len(all_iebdm_pc1) == total_n_bins
            globals()[f'{sample}']['bdm'] = all_bdm_pc1
            globals()[f'{sample}']['iebdm'] = all_iebdm_pc1

            result_fname = os.path.join(cohort_data_dir, f'{sample}.pc1.raw-concat-8')
            np.savez(result_fname, **globals()[f'{sample}'])
            #print(f'{sample}: {result_fname}.npz')

        # print number of total autosome bins in current cohort just once (at the end of the loop)

        print(f'number of bins (autosome): {total_n_bins}')
        print(f'number of samples: {len(cohort_samples)}')
        #print(f'{result_name}.npz')


# In[9]:


# 결과를 표로 정리 및 저장. 
for cpg_type in all_cpg_types:
    print(f'===\n{cpg_type}')
    globals()[f'{cpg_type}_concat_pc1_info'] = pd.DataFrame(np.zeros((len(tcga_cohorts), 2), dtype = 'int32'), index = tcga_cohorts, columns = ['n_bins', 'n_samples'])
    for cohort in tcga_cohorts:
        cohort_data_dir = os.path.join(data_dir, cpg_type, cohort, 'pc1')
        all_concat_files = glob.glob(os.path.join(cohort_data_dir, f'*raw-concat-8*'))
        globals()[f'{cpg_type}_concat_pc1_info'].loc[cohort]['n_samples'] = len(all_concat_files)
        tmp = np.load(all_concat_files[0], allow_pickle = True)
        globals()[f'{cpg_type}_concat_pc1_info'].loc[cohort]['n_bins'] = len(tmp['bdm'])
    display(globals()[f'{cpg_type}_concat_pc1_info'])
    fname = f'concat_pc1_info_{cpg_type}-8.csv'
    globals()[f'{cpg_type}_concat_pc1_info'].to_csv(os.path.join(data_dir, fname))
    print(os.path.join(data_dir, fname))


# ## concat v2: for each chromosome, standardize bdm or iebdm pc1 -> and then concatenate 22 standardized pc1 vectors
# - use this !!!

# In[10]:


for cpg_type in all_cpg_types:
    for cohort in tcga_cohorts:
        print(f'===\n{cohort}')
        cohort_data_dir = os.path.join(data_dir, cpg_type, cohort, 'pc1')
        cohort_samples = np.load(samplenames_fname,allow_pickle=True)[cohort]

        all_pc1_files = glob.glob(os.path.join(cohort_data_dir, f'*.npz'))
        bdm_pc1_files = []
        iebdm_pc1_files = []
        for f in all_pc1_files:
            if 'concat' not in f:
                if f.endswith('_inv_exp.npz'):
                    iebdm_pc1_files.append(f)
                else:
                    bdm_pc1_files.append(f)
        assert len(bdm_pc1_files) == len(iebdm_pc1_files) and len(bdm_pc1_files) == len(cohort_samples)

        for sample in cohort_samples:
            bdm_pc1_fname = os.path.join(cohort_data_dir, f'{sample}.npz')
            iebdm_pc1_fname = os.path.join(cohort_data_dir, f'{sample}_inv_exp.npz')
            bdm_pc1 = np.load(bdm_pc1_fname, allow_pickle = True)
            iebdm_pc1 = np.load(iebdm_pc1_fname, allow_pickle = True)

            globals()[f'{sample}'] = {}
            total_n_bins = 0
            
            #for chrom in AUT_LIST:
            for chrom in SMALL_AUT_LIST:
                k = f'{chrom}_pc1'
                # before adding into all_bdm_pc1, standardize pc1 vector (which is computed from a single chromsome)
                #if chrom == AUT_LIST[0]:
                if chrom == SMALL_AUT_LIST[0]:
                    all_bdm_pc1 = standardize(bdm_pc1[k])
                    all_iebdm_pc1 = standardize(iebdm_pc1[k])
                else:
                    all_bdm_pc1 = np.concatenate((all_bdm_pc1, standardize(bdm_pc1[k])), axis = None)
                    all_iebdm_pc1 = np.concatenate((all_iebdm_pc1, standardize(iebdm_pc1[k])), axis = None)
                assert len(bdm_pc1[k]) == len(iebdm_pc1[k])
                total_n_bins += len(bdm_pc1[k])
            assert len(all_bdm_pc1) == total_n_bins and len(all_iebdm_pc1) == total_n_bins
            globals()[f'{sample}']['bdm'] = all_bdm_pc1
            globals()[f'{sample}']['iebdm'] = all_iebdm_pc1

            result_fname = os.path.join(cohort_data_dir, f'{sample}.pc1.standardized-concat-8')
            np.savez(result_fname, **globals()[f'{sample}'])
            #print(f'{sample}: {result_fname}.npz')

        # print number of total autosome bins in current cohort just once (at the end of the loop)
        print(f'number of bins (autosome): {total_n_bins}')
        print(f'number of samples: {len(cohort_samples)}')


# In[ ]:




