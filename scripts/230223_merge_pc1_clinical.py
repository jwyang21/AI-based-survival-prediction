#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
#os.system('python utils.py')
from utils import *
import glob


# In[2]:


def standardize(v):
    return (v-np.mean(v))/np.std(v)


# # pipeline
# - 각 cohort에 대해:
#     - clinical data가 존재하는 샘플에 대해:
#         - bdm concatenated pc1과 clinical data를 merge
#         - iebdm concatenated pc1과 clinical data를 merge

# # 1. make a dataframe per cohort (num_samples, concatenated pc1)
# - concat all autosomes

# ## (1-1) concat raw PC1

# In[3]:


for cpg_type in all_cpg_types:
    for cohort in tcga_cohorts:
        print(f'===\n{cohort}')
        cohort_samples = np.load(samplenames_fname, allow_pickle = True)[cohort]
        cohort_data_dir = os.path.join(data_dir, cpg_type, cohort)
        cohort_pc1_data_dir = os.path.join(cohort_data_dir, 'pc1')
        cohort_bdm_pc1 = []
        cohort_iebdm_pc1 = []
        for sample in cohort_samples:
            #print(sample)
            concat_pc1_fname = glob.glob(os.path.join(cohort_pc1_data_dir, f'{sample}*raw-concat*'))
            assert len(concat_pc1_fname)==1 #for one sample, one concat file should exist. 

            bdm_pc1 = np.load(concat_pc1_fname[0], allow_pickle = True)['bdm']
            iebdm_pc1 = np.load(concat_pc1_fname[0], allow_pickle = True)['iebdm']
            #print(bdm_pc1.shape, iebdm_pc1.shape)
            cohort_bdm_pc1.append(bdm_pc1.tolist())
            cohort_iebdm_pc1.append(iebdm_pc1.tolist())
            #print('bdm', np.array(cohort_bdm_pc1).shape)
            #print('iebdm', np.array(cohort_iebdm_pc1).shape)
        cohort_bdm_pc1 = np.vstack(cohort_bdm_pc1)
        cohort_iebdm_pc1 = np.vstack(cohort_iebdm_pc1)
        #print(cohort_bdm_pc1.shape)
        #print(cohort_iebdm_pc1.shape)
        raw_cohort_bdm_pc1_df = pd.DataFrame(cohort_bdm_pc1, index = cohort_samples)
        raw_cohort_iebdm_pc1_df = pd.DataFrame(cohort_iebdm_pc1, index = cohort_samples)

        raw_cohort_bdm_pc1_df_fname = 'all_samples_pc1_concat_bdm.raw.csv'
        raw_cohort_iebdm_pc1_df_fname = 'all_samples_pc1_concat_iebdm.raw.csv'

        raw_cohort_bdm_pc1_df.to_csv(os.path.join(cohort_data_dir, raw_cohort_bdm_pc1_df_fname))
        raw_cohort_iebdm_pc1_df.to_csv(os.path.join(cohort_data_dir, raw_cohort_iebdm_pc1_df_fname))
        print(os.path.join(cohort_data_dir, raw_cohort_bdm_pc1_df_fname))
        print(os.path.join(cohort_data_dir, raw_cohort_iebdm_pc1_df_fname))


# ## (1-2) concat raw PC1s and then standardize across all chrom

# In[4]:


for cpg_type in all_cpg_types:
    for cohort in tcga_cohorts:
        print(f'===\n{cohort}')
        cohort_samples = np.load(samplenames_fname, allow_pickle = True)[cohort]
        cohort_data_dir = os.path.join(data_dir, cpg_type, cohort)
        cohort_pc1_data_dir = os.path.join(cohort_data_dir, 'pc1')
        cohort_bdm_pc1 = []
        cohort_iebdm_pc1 = []
        for sample in cohort_samples:
            #print(sample)
            concat_pc1_fname = glob.glob(os.path.join(cohort_pc1_data_dir, f'{sample}*raw-concat*'))
            assert len(concat_pc1_fname)==1 #for one sample, one concat file should exist. 

            bdm_pc1 = np.load(concat_pc1_fname[0], allow_pickle = True)['bdm']
            iebdm_pc1 = np.load(concat_pc1_fname[0], allow_pickle = True)['iebdm']
            #print(bdm_pc1.shape, iebdm_pc1.shape)
            cohort_bdm_pc1.append(standardize(bdm_pc1.tolist()))
            cohort_iebdm_pc1.append(standardize(iebdm_pc1.tolist()))
            #print('bdm', np.array(cohort_bdm_pc1).shape)
            #print('iebdm', np.array(cohort_iebdm_pc1).shape)
        cohort_bdm_pc1 = np.vstack(cohort_bdm_pc1)
        cohort_iebdm_pc1 = np.vstack(cohort_iebdm_pc1)
        #print(cohort_bdm_pc1.shape)
        #print(cohort_iebdm_pc1.shape)
        s1_cohort_bdm_pc1_df = pd.DataFrame(cohort_bdm_pc1, index = cohort_samples)
        s1_cohort_iebdm_pc1_df = pd.DataFrame(cohort_iebdm_pc1, index = cohort_samples)

        s1_cohort_bdm_pc1_df_fname = 'all_samples_pc1_concat_bdm.s_all_chrom.csv'
        s1_cohort_iebdm_pc1_df_fname = 'all_samples_pc1_concat_iebdm.s_all_chrom.csv'

        s1_cohort_bdm_pc1_df.to_csv(os.path.join(cohort_data_dir, s1_cohort_bdm_pc1_df_fname))
        s1_cohort_iebdm_pc1_df.to_csv(os.path.join(cohort_data_dir, s1_cohort_iebdm_pc1_df_fname))
        print(os.path.join(cohort_data_dir, s1_cohort_bdm_pc1_df_fname))
        print(os.path.join(cohort_data_dir, s1_cohort_iebdm_pc1_df_fname))


# ## (1-3) standardize pc1 per chrom and then concatenate

# In[6]:


for cpg_type in all_cpg_types:
    for cohort in tcga_cohorts:
        print(f'===\n{cohort}')
        cohort_samples = np.load(samplenames_fname, allow_pickle = True)[cohort]
        cohort_data_dir = os.path.join(data_dir, cpg_type, cohort)
        cohort_pc1_data_dir = os.path.join(cohort_data_dir, 'pc1')
        cohort_bdm_pc1 = []
        cohort_iebdm_pc1 = []
        for sample in cohort_samples:
            #print(sample)
            concat_pc1_fname = glob.glob(os.path.join(cohort_pc1_data_dir, f'{sample}*standardized-concat*'))
            assert len(concat_pc1_fname)==1 #for one sample, one concat file should exist. 

            bdm_pc1 = np.load(concat_pc1_fname[0], allow_pickle = True)['bdm']
            iebdm_pc1 = np.load(concat_pc1_fname[0], allow_pickle = True)['iebdm']
            #print(bdm_pc1.shape, iebdm_pc1.shape)
            cohort_bdm_pc1.append(bdm_pc1.tolist())
            cohort_iebdm_pc1.append(iebdm_pc1.tolist())
            #print('bdm', np.array(cohort_bdm_pc1).shape)
            #print('iebdm', np.array(cohort_iebdm_pc1).shape)
        cohort_bdm_pc1 = np.vstack(cohort_bdm_pc1)
        cohort_iebdm_pc1 = np.vstack(cohort_iebdm_pc1)
        #print(cohort_bdm_pc1.shape)
        #print(cohort_iebdm_pc1.shape)
        s2_cohort_bdm_pc1_df = pd.DataFrame(cohort_bdm_pc1, index = cohort_samples)
        s2_cohort_iebdm_pc1_df = pd.DataFrame(cohort_iebdm_pc1, index = cohort_samples)

        s2_cohort_bdm_pc1_df_fname = 'all_samples_pc1_concat_bdm.s_per_chrom.csv'
        s2_cohort_iebdm_pc1_df_fname = 'all_samples_pc1_concat_iebdm.s_per_chrom.csv'

        s2_cohort_bdm_pc1_df.to_csv(os.path.join(cohort_data_dir, s2_cohort_bdm_pc1_df_fname))
        s2_cohort_iebdm_pc1_df.to_csv(os.path.join(cohort_data_dir, s2_cohort_iebdm_pc1_df_fname))
        print(os.path.join(cohort_data_dir, s2_cohort_bdm_pc1_df_fname))
        print(os.path.join(cohort_data_dir, s2_cohort_iebdm_pc1_df_fname))


# ---

# # 2. merge concat_pc1 data and TCGA clinical data
# - TCGA clinical data: OS, OS.time, DSS, DSS.time, PFI, PFI.time, DFI, DFI.time
# 

# ## (2-1) remove samples which duplicate by the first 12 characters of the TCGA barcode 

# In[ ]:


all_version_fnames = ['all_samples_pc1_concat_bdm.s_all_chrom.csv', 'all_samples_pc1_concat_iebdm.s_all_chrom.csv', 'all_samples_pc1_concat_bdm.s_per_chrom.csv', 'all_samples_pc1_concat_iebdm.raw.csv', 'all_samples_pc1_concat_iebdm.s_per_chrom.csv', 'all_samples_pc1_concat_bdm.raw.csv']
# all versions of concatenating pc1 vectors of each sample


# In[20]:


all_raw_version_fnames = []
all_s_version_fnames = []
for f in all_version_fnames:
    if 'raw' in f:
        all_raw_version_fnames.append(f)
    else:
        all_s_version_fnames.append(f)


# - version 1: input covariates are age + gender + concat PC1 vectors

#     - version 1-1: raw concat PC1 + raw age + gender

# In[11]:





# In[12]:


# re-load clinical data
clinical_fname = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.xlsx'
clinical = pd.read_excel(clinical_fname, sheet_name = 'TCGA-CDR', index_col = 0)
clinical_target_columns = ['OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'age', 'gender=FEMALE']#age_at_initial_pathologic_diagnosis
survival_targets = ['OS', 'DSS', 'DFI', 'PFI']


# In[13]:


clinical_all_columns = ['bcr_patient_barcode', 'type', 'age_at_initial_pathologic_diagnosis', 'gender', 'race', 'ajcc_pathologic_tumor_stage', 'clinical_stage',                         'histological_type', 'histological_grade', 'initial_pathologic_dx_year', 'menopause_status', 'birth_days_to', 'vital_status', 'tumor_status',                         'last_contact_days_to', 'death_days_to', 'cause_of_death', 'new_tumor_event_type', 'new_tumor_event_site', 'new_tumor_event_site_other', 'new_tumor_event_dx_days_to',                         'treatment_outcome_first_course', 'margin_status', 'residual_tumor', 'OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'Redaction']


# In[14]:





# In[15]:


#clinical_target_columns = ['OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'age_at_initial_pathologic_diagnosis', 'gender']
clinical.rename(columns = {'age_at_initial_pathologic_diagnosis':'age'}, inplace = True)
clinical = clinical[clinical.age.notnull() & clinical.gender.notnull()].copy()
gender_dict = {'MALE': 0, 'FEMALE': 1}
clinical['gender'] = clinical['gender'].map(gender_dict) 
clinical.rename(columns = {'gender':'gender=FEMALE'},inplace=True) 


# In[16]:





# In[17]:





# In[18]:


clinical.head(3)


# In[21]:


for cpg_type in all_cpg_types:
    print(f'===\n{cpg_type}')
    for cohort in tcga_cohorts:
        print(f'---\n{cohort}')
        savedir = f'/data/project/3dith/data/{cpg_type}/{cohort}'
        #for version_fname in all_version_fnames:
        for version_fname in all_raw_version_fnames:
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




