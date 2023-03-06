#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
#os.system('python utils.py')
from utils import *
import glob


# # pipeline
# - 각 cohort에 대해:
#     - 각 sample의 각 autosome의 average opensea beta value를 survival target (OS / DSS / DFI / PFI)의 event 및 time 데이터와 merge 
#         - age + gender 버전일 경우 age, gender 데이터도 merge.

# In[2]:


SAMPLE_NAME_FILE = '/data/project/3dith/data/samplenames.npz'#item: {cohort}


# In[3]:


CHR_LIST = [f'chr{i}' for i in np.arange(1, 23)]


# In[4]:


ALL_COHORTS = 'TCGA-BLCA TCGA-BRCA TCGA-CHOL TCGA-COAD TCGA-KIRC TCGA-KIRP TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PAAD TCGA-PRAD TCGA-THCA TCGA-UCEC'.split(' ')


# In[5]:


cpg_type = 'opensea'
all_cpg_types = ['opensea', 'island', 'shelf_shore']


# In[6]:


result_dir = f'/data/project/3dith/data/{cpg_type}/'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


# In[7]:


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


# # 1. merge average opensea beta value and TCGA clinical data
# - TCGA clinical data: OS, OS.time, DSS, DSS.time, PFI, PFI.time, DFI, DFI.time

# ## (1-1) remove samples which duplicate by the first 12 characters of the TCGA barcode 

# - version 12.1: average opensea beta value of each autosome + age + gender

# In[3]:


# re-load clinical data
clinical_fname = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.xlsx'
clinical = pd.read_excel(clinical_fname, sheet_name = 'TCGA-CDR', index_col = 0)
clinical_target_columns = ['OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'age', 'gender=FEMALE']#age_at_initial_pathologic_diagnosis
survival_targets = ['OS', 'DSS', 'DFI', 'PFI']


# In[21]:


version = '12.1'
import_filename = 'avg_opensea_beta_autosome.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_filename), index_col = 0).sample(frac=1)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)
    
    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    clinical_cohort = clinical_cohort.loc[intersecting_samples].copy()
    
    clinical_cohort.rename(columns = {'age_at_initial_pathologic_diagnosis':'age'}, inplace = True)
    clinical_cohort = clinical_cohort[clinical_cohort.age.notnull() & clinical_cohort.gender.notnull()].copy()
    clinical_cohort['age'] /= max_age
    clinical_cohort['gender'] = clinical_cohort['gender'].map(gender_dict) 
    clinical_cohort.rename(columns = {'gender':'gender=FEMALE'},inplace=True) 
    
    for target in survival_targets:
        print(f'---\nTarget: {target}')
        current_target_col = [target, f'{target}.time', 'age', 'gender=FEMALE']
        clinical_tmp = clinical_cohort[current_target_col].copy()
        feature_clinical = all_df.merge(clinical_tmp, left_index=True, right_index=True)
        feature_clinical = feature_clinical[feature_clinical[target].notnull() & feature_clinical[f'{target}.time'].notnull()].copy()
        feature_clinical = feature_clinical.sample(frac = 1)
        feature_clinical.dropna(inplace = True)
        
        X = feature_clinical.drop([target], axis=1).copy().values
        y = feature_clinical[target].values.flatten()
        skf = StratifiedKFold(n_splits = num_folds)
        
        new_colnames = np.append(feature_clinical.drop([target], axis = 1).columns.values, target)#so that the target can be at the rightmost column
        
        for j, (train_valid_index, test_index) in enumerate(skf.split(X, y)):
            print(f'\nFold{j+1}')

            train_valid_data = np.hstack((X[train_valid_index], y[train_valid_index].reshape(-1,1)))
            train_valid_df = pd.DataFrame(train_valid_data, columns = new_colnames, index = feature_clinical.index.values[train_valid_index])
            train_valid_fname = os.path.join(cohort_save_dir, f'{target}_train_valid_fold{j+1}')
            train_valid_df2 = train_valid_df.rename({f'{target}':'Event', f'{target}.time':'Time'}, axis = 'columns')
            train_valid_df2.to_csv(train_valid_fname)
            del(train_valid_df2)
            
            X2 = train_valid_df.drop([target], axis = 1).copy().values
            y2 = train_valid_df[target].values.flatten()
            skf2 = StratifiedKFold(n_splits = num_folds_train_valid)
            
            test_data = np.hstack((X[test_index], y[test_index].reshape(-1,1)))
            test_df = pd.DataFrame(test_data, columns = new_colnames, index = feature_clinical.index.values[test_index])
            
            test_fname = os.path.join(cohort_save_dir, f'{target}_test_fold{j+1}.csv')
            
            test_df.rename({f'{target}':'Event', f'{target}.time':'Time'}, axis='columns', inplace = True)
            test_df.to_csv(test_fname)
            #print(test_df['Event'].sum()/test_df.shape[0])
            print(f'test fold {j+1}: {test_fname}')
            
            for k, (train_index, valid_index) in enumerate(skf2.split(X2, y2)):
                if k == 0:
                    train_data = np.hstack((X2[train_index], y2[train_index].reshape(-1,1)))
                    valid_data = np.hstack((X2[valid_index], y2[valid_index].reshape(-1,1)))
                    train_df = pd.DataFrame(train_data, columns = new_colnames, index = train_valid_df.index.values[train_index])
                    valid_df = pd.DataFrame(valid_data, columns = new_colnames, index = train_valid_df.index.values[valid_index])
                    
                    train_df.rename({f'{target}':'Event', f'{target}.time':'Time'}, axis = 'columns', inplace = True)
                    valid_df.rename({f'{target}':'Event', f'{target}.time':'Time'}, axis = 'columns', inplace = True)
                    
                    train_fname = os.path.join(cohort_save_dir, f'{target}_train_fold{j+1}.csv')
                    valid_fname = os.path.join(cohort_save_dir, f'{target}_valid_fold{j+1}.csv')
                    
                    train_df.to_csv(train_fname)
                    valid_df.to_csv(valid_fname)
                    print(f'train fold{j+1}: {train_fname}')
                    print(f'valid fold{j+1}: {valid_fname}')
                    #print(train_df['Event'].sum()/train_df.shape[0])
                    #print(valid_df['Event'].sum()/valid_df.shape[0])


#     - version 1-2: standardized concat PC1 + normalized age + gender

# In[3]:


# re-load clinical data
clinical_fname = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.xlsx'
clinical = pd.read_excel(clinical_fname, sheet_name = 'TCGA-CDR', index_col = 0)


# In[4]:


clinical_target_columns = ['OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'age', 'gender=FEMALE']#age_at_initial_pathologic_diagnosis
survival_targets = ['OS', 'DSS', 'DFI', 'PFI']


# In[5]:


clinical_all_columns = ['bcr_patient_barcode', 'type', 'age_at_initial_pathologic_diagnosis', 'gender', 'race', 'ajcc_pathologic_tumor_stage', 'clinical_stage',                         'histological_type', 'histological_grade', 'initial_pathologic_dx_year', 'menopause_status', 'birth_days_to', 'vital_status', 'tumor_status',                         'last_contact_days_to', 'death_days_to', 'cause_of_death', 'new_tumor_event_type', 'new_tumor_event_site', 'new_tumor_event_site_other', 'new_tumor_event_dx_days_to',                         'treatment_outcome_first_course', 'margin_status', 'residual_tumor', 'OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'Redaction']


# In[6]:


all_version_fnames = ['all_samples_pc1_concat_bdm.s_all_chrom.csv', 'all_samples_pc1_concat_iebdm.s_all_chrom.csv', 'all_samples_pc1_concat_bdm.s_per_chrom.csv', 'all_samples_pc1_concat_iebdm.raw.csv', 'all_samples_pc1_concat_iebdm.s_per_chrom.csv', 'all_samples_pc1_concat_bdm.raw.csv']
# all versions of concatenating pc1 vectors of each sample


# In[7]:


#clinical_target_columns = ['OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'age_at_initial_pathologic_diagnosis', 'gender']
clinical.rename(columns = {'age_at_initial_pathologic_diagnosis':'age'}, inplace = True)


# In[8]:


clinical = clinical[clinical.age.notnull() & clinical.gender.notnull()].copy()


# In[9]:


gender_dict = {'MALE': 0, 'FEMALE': 1}
clinical['gender'] = clinical['gender'].map(gender_dict) 
clinical.rename(columns = {'gender':'gender=FEMALE'},inplace=True) 


# In[41]:


clinical['age'] /= 100


# In[42]:


clinical.head(3)


# In[43]:


for cpg_type in all_cpg_types:
    print(f'===\n{cpg_type}')
    for cohort in tcga_cohorts:
        print(f'---\n{cohort}')
        savedir = f'/data/project/3dith/data/{cpg_type}/{cohort}'
        # for version_fname in all_version_fnames:
        for version_fname in all_s_version_fnames:
            print(f'--\nversion: {version_fname}')
            tmp = pd.read_csv(os.path.join(data_dir, cpg_type, cohort, version_fname), index_col = 0)#현재 cohort의 현재 버전 데이터
            
            # concatenated pc1 데이터에서는 TCGA barcode가 full barcode (길이 15)
            # TCGA clinical data에서는 TCGA barcode가 short barcode (길이 12)
            # --> 이 둘을 맞춰줘야 함. 
            short_barcode = [x[:12] for x in tmp.index.values]
            tmp['barcode'] = short_barcode
            tmp.drop_duplicates(subset = ['barcode'], inplace=True)
            tmp['full_barcode'] = tmp.index.values
            tmp.index = tmp['barcode'].values
            tmp.drop(['barcode'], axis=1, inplace = True)
            
            cohort_type = cohort.split('-')[1]
            clinical_cohort = clinical[clinical['type']==cohort_type].copy()
            clinical_cohort.index = clinical_cohort.bcr_patient_barcode.values #bcr_patient_barcode는 short barcode format
            clinical_cohort_target = clinical_cohort[clinical_target_columns].copy()
            
            intersection = np.intersect1d(tmp.index.values, clinical_cohort.index.values)
            
            print(f'{len(intersection)} samples have both pc1 and clinical data')
            
            # merge TCGA clinical data and concatenated pc1 data
            merged_all = pd.merge(tmp.loc[intersection].copy(), clinical_cohort.loc[intersection].copy(), left_index=True, right_index=True)
            merged_target = pd.merge(tmp.loc[intersection].copy(), clinical_cohort_target.loc[intersection].copy(), left_index=True, right_index=True)
            merged_target.index = merged_target['full_barcode'].values.flatten()
            merged_target.drop(['full_barcode'], axis=1, inplace=True)
            
            # merged_all 저장.  #version2와 동일하므로 여기선 저장안하고 v2에서 저장하자
            #merged_all_fname = os.path.join(savedir, f'merge_all_clinical_{version_fname}')# include all columns from clinical data
            #merged_all.to_csv(merged_all_fname)
            #print(f'merged all clinical data and concatenated pc1: {merged_all_fname}')
            
            # merged_target 전체 저장
            merged_target_fname = os.path.join(savedir, f'merge_target_clinical_w_age_gender_{version_fname}')# include target columns from clinical data
            merged_target.to_csv(merged_target_fname)
            print(f'merged target clinical data and concatenate pc1: {merged_target_fname}')
            
            # merged_target을 각 survival event마다 별개의 dataframe으로 나누기. 
            #merged_target_dictionary = {}
            for target in survival_targets:
                current_df = merged_target.copy()
                drop_col = survival_targets.copy()
                drop_col.remove(target)
                for col in drop_col:
                    current_df.drop(col, axis=1, inplace=True)
                    current_df.drop(f'{col}.time', axis=1, inplace=True)
                #merged_target_dictionary[target] = current_df.dropna()
                #del(current_df)
                current_df.dropna(inplace = True)
                current_df.rename({f'{target}':'Event', f'{target}.time':'Time'}, axis='columns', inplace=True)
                current_df_fname = os.path.join(savedir, f'merge_{target}_w_age_gender_{version_fname}')
                current_df.to_csv(current_df_fname)
                print(f'{target}: {current_df_fname}')            


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ---
# - version 2: use concat PC1 data only (not use age, gender)

# In[ ]:


# re-load clinical data
clinical_fname = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.xlsx'
clinical = pd.read_excel(clinical_fname, sheet_name = 'TCGA-CDR', index_col = 0)


# In[3]:


clinical_target_columns = ['OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time']# 실제 생존분석에 필요한 column들. 

survival_targets = ['OS', 'DSS', 'DFI', 'PFI']


# In[4]:


clinical_all_columns = ['bcr_patient_barcode', 'type', 'age_at_initial_pathologic_diagnosis', 'gender', 'race', 'ajcc_pathologic_tumor_stage', 'clinical_stage',                         'histological_type', 'histological_grade', 'initial_pathologic_dx_year', 'menopause_status', 'birth_days_to', 'vital_status', 'tumor_status',                         'last_contact_days_to', 'death_days_to', 'cause_of_death', 'new_tumor_event_type', 'new_tumor_event_site', 'new_tumor_event_site_other', 'new_tumor_event_dx_days_to',                         'treatment_outcome_first_course', 'margin_status', 'residual_tumor', 'OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'Redaction']


# In[5]:


all_version_fnames = ['all_samples_pc1_concat_bdm.s_all_chrom.csv', 'all_samples_pc1_concat_iebdm.s_all_chrom.csv', 'all_samples_pc1_concat_bdm.s_per_chrom.csv', 'all_samples_pc1_concat_iebdm.raw.csv', 'all_samples_pc1_concat_iebdm.s_per_chrom.csv', 'all_samples_pc1_concat_bdm.raw.csv']
# all versions of concatenating pc1 vectors of each sample


# In[ ]:





# In[5]:





# In[6]:


for cpg_type in all_cpg_types:
    print(f'===\n{cpg_type}')
    for cohort in tcga_cohorts:
        print(f'---\n{cohort}')
        savedir = f'/data/project/3dith/data/{cpg_type}/{cohort}'
        for version_fname in all_version_fnames:
            print(f'--\nversion: {version_fname}')
            tmp = pd.read_csv(os.path.join(data_dir, cpg_type, cohort, version_fname), index_col = 0)#현재 cohort의 현재 버전 데이터
            
            # concatenated pc1 데이터에서는 TCGA barcode가 full barcode (길이 15)
            # TCGA clinical data에서는 TCGA barcode가 short barcode (길이 12)
            # --> 이 둘을 맞춰줘야 함. 
            short_barcode = [x[:12] for x in tmp.index.values]
            tmp['barcode'] = short_barcode
            tmp.drop_duplicates(subset = ['barcode'], inplace=True)
            tmp['full_barcode'] = tmp.index.values
            tmp.index = tmp['barcode'].values
            tmp.drop(['barcode'], axis=1, inplace = True)
            
            cohort_type = cohort.split('-')[1]
            clinical_cohort = clinical[clinical['type']==cohort_type].copy()
            clinical_cohort.index = clinical_cohort.bcr_patient_barcode.values #bcr_patient_barcode는 short barcode format
            clinical_cohort_target = clinical_cohort[clinical_target_columns].copy()
            
            intersection = np.intersect1d(tmp.index.values, clinical_cohort.index.values)
            
            print(f'{len(intersection)} samples have both pc1 and clinical data')
            
            # merge TCGA clinical data and concatenated pc1 data
            merged_all = pd.merge(tmp.loc[intersection].copy(), clinical_cohort.loc[intersection].copy(), left_index=True, right_index=True)
            merged_target = pd.merge(tmp.loc[intersection].copy(), clinical_cohort_target.loc[intersection].copy(), left_index=True, right_index=True)
            merged_target.index = merged_target['full_barcode'].values.flatten()
            merged_target.drop(['full_barcode'], axis=1, inplace=True)
            
            # merged_all 저장.  
            merged_all_fname = os.path.join(savedir, f'merge_all_clinical_{version_fname}')# include all columns from clinical data
            merged_all.to_csv(merged_all_fname)
            print(f'merged all clinical data and concatenated pc1: {merged_all_fname}')
            
            # merged_target 전체 저장
            merged_target_fname = os.path.join(savedir, f'merge_target_clinical_{version_fname}')# include target columns from clinical data
            merged_target.to_csv(merged_target_fname)
            print(f'merged target clinical data and concatenate pc1: {merged_target_fname}')
            
            # merged_target을 각 survival event마다 별개의 dataframe으로 나누기. 
            #merged_target_dictionary = {}
            for target in survival_targets:
                current_df = merged_target.copy()
                drop_col = survival_targets.copy()
                drop_col.remove(target)
                for col in drop_col:
                    current_df.drop(col, axis=1, inplace=True)
                    current_df.drop(f'{col}.time', axis=1, inplace=True)
                #merged_target_dictionary[target] = current_df.dropna()
                #del(current_df)
                current_df.dropna(inplace = True)
                current_df.rename({f'{target}':'Event', f'{target}.time':'Time'}, axis='columns', inplace=True)
                current_df_fname = os.path.join(savedir, f'merge_{target}_{version_fname}')
                current_df.to_csv(current_df_fname)
                print(f'{target}: {current_df_fname}')            


# -------

# # End

# In[82]:


# concatenated pc1 데이터에서는 TCGA barcode가 full barcode (길이 15)
# TCGA clinical data에서는 TCGA barcode가 short barcode (길이 12)
# --> 이 둘을 맞춰줘야 함. 
short_barcode = [x[:12] for x in tmp.index.values]
tmp['barcode'] = short_barcode
tmp.drop_duplicates(subset = ['barcode'], inplace=True)
tmp['full_barcode'] = tmp.index.values
tmp.index = tmp['barcode'].values
tmp.drop(['barcode'], axis=1, inplace = True)


# In[83]:


cohort_type = cohort.split('-')[1]
clinical_cohort = clinical[clinical['type']==cohort_type].copy()
clinical_cohort.index = clinical_cohort.bcr_patient_barcode.values #bcr_patient_barcode는 short barcode format
clinical_cohort_target = clinical_cohort[clinical_target_columns].copy()


# In[84]:


intersection = np.intersect1d(tmp.index.values, clinical_cohort.index.values)


# In[85]:


print(f'{len(intersection)} samples have both pc1 and clinical data')


# In[86]:


# merge TCGA clinical data and concatenated pc1 data
merged_all = pd.merge(tmp.loc[intersection].copy(), clinical_cohort.loc[intersection].copy(), left_index=True, right_index=True)
merged_target = pd.merge(tmp.loc[intersection].copy(), clinical_cohort_target.loc[intersection].copy(), left_index=True, right_index=True)
merged_target.index = merged_target['full_barcode'].values.flatten()
merged_target.drop(['full_barcode'], axis=1, inplace=True)


# In[90]:


# 일단 merged_all 먼저 저장.  
merged_all_fname = os.path.join(savedir, f'merge_all_clinical_{version_fname}')# include all columns from clinical data
merged_all.to_csv(merged_all_fname)
print(merged_all_fname)


# In[92]:


# merged_target 전체 저장
merged_target_fname = os.path.join(savedir, f'merge_target_clinical_{version_fname}')# include target columns from clinical data
merged_target.to_csv(merged_target_fname)
print(merged_target_fname)


# In[94]:


merged_target.head(3)


# In[122]:


#merged_target_dictionary = {}
for target in survival_targets:
    current_df = merged_target.copy()
    drop_col = survival_targets.copy()
    drop_col.remove(target)
    for col in drop_col:
        current_df.drop(col, axis=1, inplace=True)
        current_df.drop(f'{col}.time', axis=1, inplace=True)
    #merged_target_dictionary[target] = current_df.dropna()
    #del(current_df)
    current_df.rename({f'{target}':'Event', f'{target}.time':'Time'}, axis='columns', inplace=True)
    current_df_fname = os.path.join(savedir, f'merge_{target}_{version_fname}')
    current_df.to_csv(current_df_fname)
    print(f'{target}: {current_df_fname}')


# In[ ]:




