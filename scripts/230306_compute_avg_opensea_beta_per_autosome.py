#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import os
#os.system('python utils.py')
from utils import *
import glob


# # pipeline
# - 각 cohort에 대해:
#     - 각 autosome 별 average opensea beta value 계산

# In[3]:


#def standardize(v):
#    return (v-np.mean(v))/np.std(v)


# In[38]:


#CPG_ANNOT = pd.read_csv('/data/project/3dith/data/humanmethylation450_15017482_v1-2.csv', skiprows = [0,1,2,3,4,5,6], index_col=0)
CPG_METADATA = pd.read_csv('/data/project/3dith/data/450k_metadata.open_sea.sorted.bed', header=None, sep = '\t') # columns = ['chrom', 'start', 'end', 'cpg'] # use this!!
# CPG_METADAT가 CPG_ANNOT으로부터 (1) 이름이 'cg'로 시작하고 (2) opensea CpG probe인 것 만 골라내서 bed file로 정리한 것.
CPG_METADATA.columns = ['chrom', 'start', 'end', 'cpg']


# In[103]:


SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}


# In[75]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]


# In[39]:


print(f'CPG_ANNOT: {CPG_ANNOT.shape}')
print(f'CPG_METADATA: {CPG_METADATA.shape}')


# In[105]:


ALL_COHORTS = 'TCGA-BLCA TCGA-BRCA TCGA-CHOL TCGA-COAD TCGA-KIRC TCGA-KIRP TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PAAD TCGA-PRAD TCGA-THCA TCGA-UCEC'.split(' ')


# In[137]:


cpg_type = 'opensea'
result_dir = f'/data/project/3dith/data/{cpg_type}/'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


# In[104]:


def get_sample_list(cohort):
    # sample list of input TCGA cohort
    samples = np.load(SAMPLE_NAME_FILE)[cohort]
    S = samples.tolist()
    if cohort=='PCBC':
        T = []
        N = []
    else: #TCGA cohort
        T = []
        N = []
        for s in samples:
            if int(s[13:15]) >= 1 and int(s[13:15]) <= 9: #tumor barcode: '01' ~ '09'
                T.append(s)
            elif int(s[13:15]) >=10 and int(s[13:15]) <= 19:
                N.append(s)
            else:
                pass
    return T, N, S
    # usage: T, N, S = get_sample_list(args.cohort)


# In[140]:


def get_cpg_list_per_chrom(chr_list, cpg_type):
    if cpg_type == 'opensea':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.open_sea.sorted.bed', sep = '\t', header = None) #chrom // start // end // cpg_id
    elif cpg_type == 'island':
        df = pd.read_csv('/data/project/3dith/data/450k_metadata.island.sorted.bed', sep = '\t', header = None) #chrom // start // end // cpg_id
    else: #일단 shelf, shore, shelr_shore는 보류.
        raise Exception("Wrong cpg_type!")

    df.columns = ['chrom', 'start', 'end', 'cpg']
    total_list = {}
    for chrom in chr_list:
        current_chrom_cpgs = df[df['chrom'] == chrom].copy().cpg.values.flatten()
        total_list[chrom] = current_chrom_cpgs

    return total_list


# In[141]:


def get_avg_beta_per_chrom(cohort, cpg_list, S, chr_list):
    # 현재 cohort의 전체 beta 데이터 중에서, 입력받은 cpg_list들의 beta value들만 반환.
    beta_fname = f'/data/project/3dith/data/450k_xena/{cohort}.HumanMethylation450.tsv' #이 파일로부터 beta value를 읽어옴.
    beta = pd.read_csv(beta_fname, sep = '\t', index_col = 0) # row: CpG probe, column: sample
    
    all_avg_opensea_beta = pd.DataFrame(index = beta.columns.values) # average opensea beta value per chrom, per sample. 
    
    for chrom in chr_list:
        beta_target_cpg_df = beta.loc[cpg_list[chrom]].dropna() # opensea beta values of current chromosomes for all samples
        #print('beta_target_cpg_df')
        #display(beta_target_cpg_df.head(3))
        #print(beta_target_cpg_df.shape)
        #print("\n")
        avg_target_cpg_beta = beta_target_cpg_df.mean().values # average opensea beta value of current chromosome for all samples 
        #print('avg_target_cpg_beta')
        #print(avg_target_cpg_beta.shape)
        for i in range(beta_target_cpg_df.shape[1]):
            assert beta_target_cpg_df.columns.values[i] == beta.columns.values[i]
        all_avg_opensea_beta[chrom] = avg_target_cpg_beta
    
    return all_avg_opensea_beta.loc[S].copy()


# # 0. Check CpG probe list
# - CPG_METADATA의 cpg probe들이 CPG_ANNOT에서 추출된 게 맞는지 probe 총 개수 비교.

# ## 0.1. CPG_ANNOT의 opensea cpg probe 개수

# In[40]:


CPG_ANNOT.head(3)


# In[41]:


CPG_ANNOT.CHR.unique()


# In[42]:


tmp = CPG_ANNOT[CPG_ANNOT['CHR'].notnull()].copy()


# In[43]:


chr_cleaned = [str(x).split('.')[0] for x in tmp.CHR.values]


# In[44]:


tmp['CHR'] = chr_cleaned


# In[47]:


tmp2 = tmp[tmp['CHR']!='X'].copy()
tmp3 = tmp2[tmp2['CHR']!='Y'].copy()


# In[48]:


tmp3.shape


# In[53]:


tmp4 = tmp3[tmp3['Relation_to_UCSC_CpG_Island'].isnull()].copy()


# In[64]:


tmp4.shape


# In[67]:


cg_mask = ['cg' in tmp4.Name.values[i] for i in range(tmp4.shape[0])]


# In[68]:


tmp5 = tmp4.iloc[cg_mask,:].copy()


# In[70]:


tmp5.shape


# In[61]:


mask = [CPG_METADATA.chrom.values[i] in CHR_LIST for i in range(CPG_METADATA.shape[0])]


# ## 0.2. CPG_METADATA의 opensea cpg probe 개수

# In[62]:


CPG_METADATA_autosome = CPG_METADATA.iloc[mask,:]


# In[63]:


CPG_METADATA_autosome.shape


# # 1. Import opensea CpG probe list

# In[89]:


cpg_list = get_cpg_list_per_chrom(CHR_LIST, cpg_type) #key: chromosome name #value: opensea CpG probes in that chromosome


# # 2. Compute average opensea beta value per autosome
# - should iterate for each cohort

# In[143]:


for cohort in ALL_COHORTS:
    print(f"===\n{cohort}")
    T, N, S = get_sample_list(cohort)
    cohort_type = cohort.split('-')[1]
    globals()[cohort_type] = get_avg_beta_per_chrom(cohort, cpg_list, S, CHR_LIST)
    display(globals()[cohort_type].head(3))
    print(globals()[cohort_type].shape)
    cohort_result_dir = os.path.join(result_dir, cohort)
    if not os.path.exists(cohort_result_dir):
        os.makedirs(cohort_result_dir)
    fname = os.path.join(cohort_result_dir, f'avg_opensea_beta_autosome.csv')
    print(f'result fname: {fname}')
    globals()[cohort_type].to_csv(fname)


# ---
