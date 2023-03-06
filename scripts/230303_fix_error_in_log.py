#!/usr/bin/env python
# coding: utf-8

# - v11.2의 로그를 기록할 때, tarin 및 validation c-index를 best performance가 아닌 맨 마지막 epoch 때의 성능으로 기록함. 
# - 각 case마다 저장된 metric을 다시 load해서, 기록 정정
# - bhi4의 v11.2, bhi-gpu3의 v9.1

# In[1]:


import pandas as pd
import numpy as np
import os
import glob
import pickle


# In[2]:


all_lrs = [0.0005, 0.0003, 0.0001, 0.00007, 0.00005, 0.00003, 0.00001]
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']


# In[3]:


server = 'bhi-gpu3'
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result'


# In[4]:


all_version_fnames = ['merge_PFI_w_age_gender_all_samples_pc1_concat_iebdm.s_all_chrom.csv', 'merge_DSS_w_age_gender_all_samples_pc1_concat_iebdm.raw.csv', 'merge_PFI_w_age_gender_all_samples_pc1_concat_bdm.raw.csv', 'merge_PFI_w_age_gender_all_samples_pc1_concat_iebdm.s_per_chrom.csv', 'merge_OS_w_age_gender_all_samples_pc1_concat_bdm.s_per_chrom.csv', 'merge_DFI_w_age_gender_all_samples_pc1_concat_bdm.s_all_chrom.csv', 'merge_DFI_w_age_gender_all_samples_pc1_concat_iebdm.raw.csv', 'merge_PFI_w_age_gender_all_samples_pc1_concat_bdm.s_per_chrom.csv', 'merge_OS_w_age_gender_all_samples_pc1_concat_iebdm.raw.csv', 'merge_PFI_w_age_gender_all_samples_pc1_concat_iebdm.raw.csv', 'merge_DSS_w_age_gender_all_samples_pc1_concat_iebdm.s_per_chrom.csv', 'merge_DFI_w_age_gender_all_samples_pc1_concat_bdm.raw.csv', 'merge_DFI_w_age_gender_all_samples_pc1_concat_iebdm.s_all_chrom.csv', 'merge_DSS_w_age_gender_all_samples_pc1_concat_bdm.s_per_chrom.csv', 'merge_DSS_w_age_gender_all_samples_pc1_concat_bdm.s_all_chrom.csv', 'merge_DSS_w_age_gender_all_samples_pc1_concat_iebdm.s_all_chrom.csv', 'merge_OS_w_age_gender_all_samples_pc1_concat_iebdm.s_all_chrom.csv', 'merge_DFI_w_age_gender_all_samples_pc1_concat_bdm.s_per_chrom.csv', 'merge_DFI_w_age_gender_all_samples_pc1_concat_iebdm.s_per_chrom.csv', 'merge_OS_w_age_gender_all_samples_pc1_concat_iebdm.s_per_chrom.csv', 'merge_OS_w_age_gender_all_samples_pc1_concat_bdm.raw.csv', 'merge_OS_w_age_gender_all_samples_pc1_concat_bdm.s_all_chrom.csv', 'merge_DSS_w_age_gender_all_samples_pc1_concat_bdm.raw.csv', 'merge_PFI_w_age_gender_all_samples_pc1_concat_bdm.s_all_chrom.csv']


# # v1

# In[5]:


version = 1

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v1.2

# In[6]:


version = 1.2

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v2

# In[10]:


version = 2

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v2.2

# In[13]:


version = 2.2

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v3
# - re-running in bhi-gpu1 (23/03/03)

# In[14]:


version = 3

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v3.2

# In[25]:


version = 3.2

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v4
# - running in bhi-gpu1 (23/03/03)

# In[17]:


version = 4

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v4.2

# In[18]:


version = 4.2

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v5

# In[19]:


version = 5

for cohort in all_cohorts:
    print(f'===\n{cohort}')
    #print(f'---\nlr: {lr}')
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
    log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
    #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
    if len(log_list) ==1:
        cohort_log = log_list[0]
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
        print(f"--\nOriginal df:")
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3)) 
        original_columns = globals()[f'{cohort_type}'].columns.values

        new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
        new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
        
        all_items = [x.split('.csv')[0] for x in globals()[f'{cohort_type}'].index.values]
        #for current_item in globals()[f'{cohort_type}'].index.values:
        for current_item in all_items:
            #print(f"-\nCurrent item: {current_item}")
            #metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
            metrics_fname = os.path.join(cohort_dir,  f'{current_item}_metrics')
            assert os.path.exists(metrics_fname)
            with open(metrics_fname, 'rb') as f:
                metrics = pickle.load(f)
            f.close()

            train_tmp = []
            valid_tmp = []
            for l in range(len(metrics['c-index'])):
                train_tmp.append(metrics['c-index'][l][1])
            for l in range(len(metrics['valid_c-index'])):
                valid_tmp.append(metrics['valid_c-index'][l][1])
            #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
            #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')

            new_train_c_df.loc[f'{current_item}.csv'] = np.max(train_tmp)
            new_valid_c_df.loc[f'{current_item}.csv'] = np.max(valid_tmp)

            # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['train_c'] - train_tmp[-1]) <= 1e-10
            #assert (globals()[f'{cohort_type}'].loc[f'{current_item}.csv']['valid_c'] - valid_tmp[-1]) <= 1e-10

        globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)

        print(f"--\nChanged df:")
        globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
        display(globals()[f'{cohort_type}'].head(3))
        display(globals()[f'{cohort_type}'].tail(3))
        globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# ---
# # v6 ~ : not needed (fixed code)

# ---
# # v9.1

# In[20]:


version = 9.1


# In[24]:


for cohort in all_cohorts:
    print(f'===\n{cohort}')
    for lr in all_lrs:
        print(f'---\nlr: {lr}')
        cohort_type = cohort.split('-')[1]
        cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
        log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
        if len(log_list) ==1:
            cohort_log = log_list[0]
            globals()[f'{cohort_type}'] = pd.read_csv(cohort_log, index_col = 0)
            print(f"--\nOriginal df:")
            display(globals()[f'{cohort_type}'].head(3))
            display(globals()[f'{cohort_type}'].tail(3)) 
            original_columns = globals()[f'{cohort_type}'].columns.values
            
            new_train_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['train_c'])
            new_valid_c_df = pd.DataFrame(index = globals()[f'{cohort_type}'].index.values, columns = ['valid_c'])
            
            for current_item in globals()[f'{cohort_type}'].index.values:
                #print(f"-\nCurrent item: {current_item}")
                metrics_fname = os.path.join(cohort_dir,  f'{current_item}_lr_{lr}_metrics')
                with open(metrics_fname, 'rb') as f:
                    metrics = pickle.load(f)
                f.close()
                
                train_tmp = []
                valid_tmp = []
                for l in range(len(metrics['c-index'])):
                    train_tmp.append(metrics['c-index'][l][1])
                for l in range(len(metrics['valid_c-index'])):
                    valid_tmp.append(metrics['valid_c-index'][l][1])
                #print(f'train: max {np.max(train_tmp)}, last {train_tmp[-1]}')
                #print(f'valid: max {np.max(valid_tmp)}, last{valid_tmp[-1]}')
                
                new_train_c_df.loc[current_item] = np.max(train_tmp)
                new_valid_c_df.loc[current_item] = np.max(valid_tmp)
                
                # 기존에 저장되어 있던 metric들이 last epoch 때의 성능이 맞는지 confirm
                #assert (globals()[f'{cohort_type}'].loc[current_item]['train_c'] - train_tmp[-1]) <= 1e-10
                #assert (globals()[f'{cohort_type}'].loc[current_item]['valid_c'] - valid_tmp[-1]) <= 1e-10

            globals()[f'{cohort_type}'].drop(['train_c'], axis = 1, inplace = True)
            globals()[f'{cohort_type}'].drop(['valid_c'], axis = 1, inplace = True)
            globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_train_c_df, left_index = True, right_index = True)
            globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'].merge(new_valid_c_df, left_index = True, right_index = True)
                                
            print(f"--\nChanged df:")
            globals()[f'{cohort_type}'] = globals()[f'{cohort_type}'][original_columns]
            display(globals()[f'{cohort_type}'].head(3))
            display(globals()[f'{cohort_type}'].tail(3))
            globals()[f'{cohort_type}'].to_csv(cohort_log, index = True)


# In[ ]:




