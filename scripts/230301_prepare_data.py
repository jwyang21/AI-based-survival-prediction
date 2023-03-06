#!/usr/bin/env python
# coding: utf-8

# ## merge {OS/DSS/DFI/PFI} and {OS/DSS/DFI/PFI}.time with input features
# ## split data of each cohort into 5 folds (random split, non-overlapping)
# - stratify by the event

# In[1]:


import pandas as pd
import numpy as np
import os
import random
from sklearn.model_selection import StratifiedKFold


# In[2]:


#https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.StratifiedKFold.html


# In[3]:


num_folds = 5
num_folds_train_valid = 2 #4 #if 2: split train_valid set into half and use the resulting subsets as train set and validation set, respectively. 
seed = 42
cpg_type = 'opensea'
target_cohorts_fname = '/data/project/3dith/data/target_cohorts.manifest'
data_dir = os.path.join('/data/project/3dith/data', cpg_type)
gender_dict = {'MALE': 0, 'FEMALE': 1}
max_age = 100 # max age in TCGA clinical data


# In[4]:


np.random.seed(seed)
random.seed(seed) 
target_cohorts = pd.read_csv(target_cohorts_fname).cohort.values


# In[5]:


# re-load clinical data
clinical_fname = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.xlsx'
clinical = pd.read_excel(clinical_fname, sheet_name = 'TCGA-CDR', index_col = 0)
clinical_target_columns = ['OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'age', 'gender=FEMALE']#age_at_initial_pathologic_diagnosis
survival_targets = ['OS', 'DSS', 'DFI', 'PFI']
clinical.index = clinical.bcr_patient_barcode.values


# ---
# # v6.1
# - standardize(normal distance 22개) + standardize(stem distance 22개) + stem closeness + age + gender

# In[6]:


version = '6.1'
for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, 'optimal_dist_score_merged.csv'), index_col = 0).sample(frac=1)
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


# # v6.2
# - standardize(normal distance 22개) + standardize(stem distance 22개) + stem closeness

# In[7]:


version = '6.2'
for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, 'optimal_dist_score_merged.csv'), index_col = 0).sample(frac=1)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)
    
    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    clinical_cohort = clinical_cohort.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\nTarget: {target}')
        current_target_col = [target, f'{target}.time']
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
            test_df.dropna().to_csv(test_fname)
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
                    
                    train_df.dropna().to_csv(train_fname)
                    valid_df.dropna().to_csv(valid_fname)
                    print(f'train fold{j+1}: {train_fname}')
                    print(f'valid fold{j+1}: {valid_fname}')
                    #print(train_df['Event'].sum()/train_df.shape[0])
                    #print(valid_df['Event'].sum()/valid_df.shape[0])


# ---
# # v7.1
# - normal distance 22 + stem distance 22 + stem closeness + concatenated pc1 (full 22 autosomes, s_per_chrom version) + age + gender
# - v6.1의 데이터에 full 22 autosome의 s_per_chrom pc1 concat 데이터를 이어붙이기

# In[8]:


save_version = '7.1'
load_version = '6.1'
import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{save_version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_data_fname), index_col = 0)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)

    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\n{target}')
        for j in range(num_folds):
            print(f"\nFold{j+1}")
            df_train = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_train_fold{j+1}.csv'), index_col = 0)
            df_valid = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_valid_fold{j+1}.csv'), index_col = 0)
            df_test = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_test_fold{j+1}.csv'), index_col = 0)
            
            train_samples = np.intersect1d(all_df.index.values, df_train.index.values)
            valid_samples = np.intersect1d(all_df.index.values, df_valid.index.values)
            test_samples = np.intersect1d(all_df.index.values, df_test.index.values)
            
            df_train_merged = df_train.merge(all_df.loc[train_samples].copy(), left_index = True, right_index = True)
            df_valid_merged = df_valid.merge(all_df.loc[valid_samples].copy(), left_index = True, right_index = True)
            df_test_merged = df_test.merge(all_df.loc[test_samples].copy(), left_index = True, right_index = True)
            
            df_train_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_train_fold{j+1}.csv')
            df_valid_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_valid_fold{j+1}.csv')
            df_test_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_test_fold{j+1}.csv')
            
            df_train_merged.to_csv(df_train_merged_fname)
            df_valid_merged.to_csv(df_valid_merged_fname)
            df_test_merged.to_csv(df_test_merged_fname)
            
            print(f'train: {df_train_merged_fname}')
            print(f'valid: {df_valid_merged_fname}')
            print(f'test: {df_test_merged_fname}')


# In[ ]:





# ---
# # v7.2
# - normal distance 22 + stem distance 22 + stem closeness + concatenated pc1 (full 22 autosomes, s_per_chrom version)

# In[9]:


save_version = '7.2'
load_version = '6.2'
import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{save_version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_data_fname), index_col = 0)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)

    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\n{target}')
        for j in range(num_folds):
            print(f"\nFold{j+1}")
            df_train = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_train_fold{j+1}.csv'), index_col = 0)
            df_valid = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_valid_fold{j+1}.csv'), index_col = 0)
            df_test = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_test_fold{j+1}.csv'), index_col = 0)
            
            train_samples = np.intersect1d(all_df.index.values, df_train.index.values)
            valid_samples = np.intersect1d(all_df.index.values, df_valid.index.values)
            test_samples = np.intersect1d(all_df.index.values, df_test.index.values)
            
            df_train_merged = df_train.merge(all_df.loc[train_samples].copy(), left_index = True, right_index = True)
            df_valid_merged = df_valid.merge(all_df.loc[valid_samples].copy(), left_index = True, right_index = True)
            df_test_merged = df_test.merge(all_df.loc[test_samples].copy(), left_index = True, right_index = True)
            
            df_train_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_train_fold{j+1}.csv')
            df_valid_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_valid_fold{j+1}.csv')
            df_test_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_test_fold{j+1}.csv')
            
            df_train_merged.to_csv(df_train_merged_fname)
            df_valid_merged.to_csv(df_valid_merged_fname)
            df_test_merged.to_csv(df_test_merged_fname)
            
            print(f'train: {df_train_merged_fname}')
            print(f'valid: {df_valid_merged_fname}')
            print(f'test: {df_test_merged_fname}')


# ---
# # v8.1
# - normal distance 22개 + stem distance 22개 + stem closeness + concatenated pc1 (full, 10-averaged, s_per_chrom) + age + gender
# - import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_avg10.csv'

# In[10]:


save_version = '8.1'
load_version = '6.1'
import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_avg10.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{save_version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_data_fname), index_col = 0)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)

    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\n{target}')
        for j in range(num_folds):
            print(f"\nFold{j+1}")
            df_train = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_train_fold{j+1}.csv'), index_col = 0)
            df_valid = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_valid_fold{j+1}.csv'), index_col = 0)
            df_test = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_test_fold{j+1}.csv'), index_col = 0)
            
            train_samples = np.intersect1d(all_df.index.values, df_train.index.values)
            valid_samples = np.intersect1d(all_df.index.values, df_valid.index.values)
            test_samples = np.intersect1d(all_df.index.values, df_test.index.values)
            
            df_train_merged = df_train.merge(all_df.loc[train_samples].copy(), left_index = True, right_index = True)
            df_valid_merged = df_valid.merge(all_df.loc[valid_samples].copy(), left_index = True, right_index = True)
            df_test_merged = df_test.merge(all_df.loc[test_samples].copy(), left_index = True, right_index = True)
            
            df_train_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_train_fold{j+1}.csv')
            df_valid_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_valid_fold{j+1}.csv')
            df_test_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_test_fold{j+1}.csv')
            
            df_train_merged.to_csv(df_train_merged_fname)
            df_valid_merged.to_csv(df_valid_merged_fname)
            df_test_merged.to_csv(df_test_merged_fname)
            
            print(f'train: {df_train_merged_fname}')
            print(f'valid: {df_valid_merged_fname}')
            print(f'test: {df_test_merged_fname}')


# ---
# # v8.2
# - normal distance 22개 + stem distance 22개 + stem closeness + concatenated pc1 (full, 10-averaged, s_per_chrom)
# - import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_avg10.csv'

# In[11]:


save_version = '8.2'
load_version = '6.2'
import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_avg10.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{save_version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_data_fname), index_col = 0)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)

    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\n{target}')
        for j in range(num_folds):
            print(f"\nFold{j+1}")
            df_train = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_train_fold{j+1}.csv'), index_col = 0)
            df_valid = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_valid_fold{j+1}.csv'), index_col = 0)
            df_test = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_test_fold{j+1}.csv'), index_col = 0)
            
            train_samples = np.intersect1d(all_df.index.values, df_train.index.values)
            valid_samples = np.intersect1d(all_df.index.values, df_valid.index.values)
            test_samples = np.intersect1d(all_df.index.values, df_test.index.values)
            
            df_train_merged = df_train.merge(all_df.loc[train_samples].copy(), left_index = True, right_index = True)
            df_valid_merged = df_valid.merge(all_df.loc[valid_samples].copy(), left_index = True, right_index = True)
            df_test_merged = df_test.merge(all_df.loc[test_samples].copy(), left_index = True, right_index = True)
            
            df_train_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_train_fold{j+1}.csv')
            df_valid_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_valid_fold{j+1}.csv')
            df_test_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_test_fold{j+1}.csv')
            
            df_train_merged.to_csv(df_train_merged_fname)
            df_valid_merged.to_csv(df_valid_merged_fname)
            df_test_merged.to_csv(df_test_merged_fname)
            
            print(f'train: {df_train_merged_fname}')
            print(f'valid: {df_valid_merged_fname}')
            print(f'test: {df_test_merged_fname}')


# ---
# # v9.1
# - normal distance 22개 + stem distance 22개 + stem closeness + concatenated pc1 (full, 20-averaged, s_per_chrom) + age + gender
# - v6.1 데이터에 full autosome의 20_averaged s_perchrom pc1 데이터를 이어붙이기
# - import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_avg20.csv'

# In[12]:


save_version = '9.1'
load_version = '6.1'
import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_avg20.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{save_version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_data_fname), index_col = 0)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)

    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\n{target}')
        for j in range(num_folds):
            print(f"\nFold{j+1}")
            df_train = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_train_fold{j+1}.csv'), index_col = 0)
            df_valid = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_valid_fold{j+1}.csv'), index_col = 0)
            df_test = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_test_fold{j+1}.csv'), index_col = 0)
            
            train_samples = np.intersect1d(all_df.index.values, df_train.index.values)
            valid_samples = np.intersect1d(all_df.index.values, df_valid.index.values)
            test_samples = np.intersect1d(all_df.index.values, df_test.index.values)
            
            df_train_merged = df_train.merge(all_df.loc[train_samples].copy(), left_index = True, right_index = True)
            df_valid_merged = df_valid.merge(all_df.loc[valid_samples].copy(), left_index = True, right_index = True)
            df_test_merged = df_test.merge(all_df.loc[test_samples].copy(), left_index = True, right_index = True)
            
            df_train_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_train_fold{j+1}.csv')
            df_valid_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_valid_fold{j+1}.csv')
            df_test_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_test_fold{j+1}.csv')
            
            df_train_merged.to_csv(df_train_merged_fname)
            df_valid_merged.to_csv(df_valid_merged_fname)
            df_test_merged.to_csv(df_test_merged_fname)
            
            print(f'train: {df_train_merged_fname}')
            print(f'valid: {df_valid_merged_fname}')
            print(f'test: {df_test_merged_fname}')


# ---
# # v9.2
# - normal distance 22개 + stem distance 22개 + stem closeness + concatenated pc1 (full, 20-averaged, s_per_chrom)
# - v6.2 데이터에 full autosome의 20_averaged s_perchrom pc1 데이터를 이어붙이기
# - import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_avg20.csv'

# In[13]:


save_version = '9.2'
load_version = '6.2'
import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_avg20.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{save_version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_data_fname), index_col = 0)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)

    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\n{target}')
        for j in range(num_folds):
            print(f"\nFold{j+1}")
            df_train = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_train_fold{j+1}.csv'), index_col = 0)
            df_valid = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_valid_fold{j+1}.csv'), index_col = 0)
            df_test = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_test_fold{j+1}.csv'), index_col = 0)
            
            train_samples = np.intersect1d(all_df.index.values, df_train.index.values)
            valid_samples = np.intersect1d(all_df.index.values, df_valid.index.values)
            test_samples = np.intersect1d(all_df.index.values, df_test.index.values)
            
            df_train_merged = df_train.merge(all_df.loc[train_samples].copy(), left_index = True, right_index = True)
            df_valid_merged = df_valid.merge(all_df.loc[valid_samples].copy(), left_index = True, right_index = True)
            df_test_merged = df_test.merge(all_df.loc[test_samples].copy(), left_index = True, right_index = True)
            
            df_train_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_train_fold{j+1}.csv')
            df_valid_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_valid_fold{j+1}.csv')
            df_test_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_test_fold{j+1}.csv')
            
            df_train_merged.to_csv(df_train_merged_fname)
            df_valid_merged.to_csv(df_valid_merged_fname)
            df_test_merged.to_csv(df_test_merged_fname)
            
            print(f'train: {df_train_merged_fname}')
            print(f'valid: {df_valid_merged_fname}')
            print(f'test: {df_test_merged_fname}')


# ---
# # v10.1
# - normal distance 22개 + stem distance 22개 + stem closeness + concatenated pc1 (first 8 autosomes, s_per_chrom) + age + gender
# - v6.1 데이터에 first 8 autosome의 s_per_chrom pc1 데이터를 이어붙이기
# - import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_8.csv'

# In[14]:


save_version = '10.1'
load_version = '6.1'
import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_8.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{save_version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_data_fname), index_col = 0)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)

    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\n{target}')
        for j in range(num_folds):
            print(f"\nFold{j+1}")
            df_train = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_train_fold{j+1}.csv'), index_col = 0)
            df_valid = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_valid_fold{j+1}.csv'), index_col = 0)
            df_test = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_test_fold{j+1}.csv'), index_col = 0)
            
            train_samples = np.intersect1d(all_df.index.values, df_train.index.values)
            valid_samples = np.intersect1d(all_df.index.values, df_valid.index.values)
            test_samples = np.intersect1d(all_df.index.values, df_test.index.values)
            
            df_train_merged = df_train.merge(all_df.loc[train_samples].copy(), left_index = True, right_index = True)
            df_valid_merged = df_valid.merge(all_df.loc[valid_samples].copy(), left_index = True, right_index = True)
            df_test_merged = df_test.merge(all_df.loc[test_samples].copy(), left_index = True, right_index = True)
            
            df_train_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_train_fold{j+1}.csv')
            df_valid_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_valid_fold{j+1}.csv')
            df_test_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_test_fold{j+1}.csv')
            
            df_train_merged.to_csv(df_train_merged_fname)
            df_valid_merged.to_csv(df_valid_merged_fname)
            df_test_merged.to_csv(df_test_merged_fname)
            
            print(f'train: {df_train_merged_fname}')
            print(f'valid: {df_valid_merged_fname}')
            print(f'test: {df_test_merged_fname}')


# ---
# # v10.2
# - normal distance 22개 + stem distance 22개 + stem closeness + concatenated pc1 (first 8 autosomes, s_per_chrom)
# - v6.2 데이터에 first 8 autosome의 s_per_chrom pc1 데이터를 이어붙이기
# - import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_8.csv'

# In[15]:


save_version = '10.2'
load_version = '6.2'
import_data_fname = 'all_samples_pc1_concat_bdm.s_per_chrom_8.csv'

for i in range(len(target_cohorts)):#all
#for i in range(1):#debug
    cohort = target_cohorts[i]
    cohort_type = cohort.split('TCGA-')[1]
    print(f"===\n{cohort}")
    cohort_data_dir = os.path.join(data_dir, cohort)
    cohort_save_dir = os.path.join(data_dir, cohort, f'v{save_version}')
    if not os.path.exists(cohort_save_dir):
        os.makedirs(cohort_save_dir)
        
    all_df = pd.read_csv(os.path.join(cohort_data_dir, import_data_fname), index_col = 0)
    short_barcode = [x[:12] for x in all_df.index.values]
    all_df['barcode'] = short_barcode
    all_df.drop_duplicates(subset = ['barcode'], inplace=True)
    all_df.index = all_df['barcode'].values.flatten()
    all_df.drop(['barcode'], axis = 1, inplace = True)

    clinical_cohort = clinical[clinical['type']==cohort_type].copy()
    intersecting_samples = np.intersect1d(all_df.index.values, clinical_cohort.index.values)
    all_df = all_df.loc[intersecting_samples].copy()
    
    for target in survival_targets:
        print(f'---\n{target}')
        for j in range(num_folds):
            print(f"\nFold{j+1}")
            df_train = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_train_fold{j+1}.csv'), index_col = 0)
            df_valid = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_valid_fold{j+1}.csv'), index_col = 0)
            df_test = pd.read_csv(os.path.join(cohort_data_dir, f'v{load_version}', f'{target}_test_fold{j+1}.csv'), index_col = 0)
            
            train_samples = np.intersect1d(all_df.index.values, df_train.index.values)
            valid_samples = np.intersect1d(all_df.index.values, df_valid.index.values)
            test_samples = np.intersect1d(all_df.index.values, df_test.index.values)
            
            df_train_merged = df_train.merge(all_df.loc[train_samples].copy(), left_index = True, right_index = True)
            df_valid_merged = df_valid.merge(all_df.loc[valid_samples].copy(), left_index = True, right_index = True)
            df_test_merged = df_test.merge(all_df.loc[test_samples].copy(), left_index = True, right_index = True)
            
            df_train_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_train_fold{j+1}.csv')
            df_valid_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_valid_fold{j+1}.csv')
            df_test_merged_fname = os.path.join(cohort_data_dir, f'v{save_version}', f'{target}_test_fold{j+1}.csv')
            
            df_train_merged.to_csv(df_train_merged_fname)
            df_valid_merged.to_csv(df_valid_merged_fname)
            df_test_merged.to_csv(df_test_merged_fname)
            
            print(f'train: {df_train_merged_fname}')
            print(f'valid: {df_valid_merged_fname}')
            print(f'test: {df_test_merged_fname}')


# ---
# # v11.1
# - final normal/stem distance & stem closeness + age + gender

# In[16]:


version = '11.1'
import_filename = 'final_distance_score.csv'

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


# ---
# # v11.2
# - final normal/stem distance & stem closeness 

# In[17]:


version = '11.2'
import_filename = 'final_distance_score.csv'

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
    
    for target in survival_targets:
        print(f'---\nTarget: {target}')
        current_target_col = [target, f'{target}.time']
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
            test_df.dropna().to_csv(test_fname)
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
                    
                    train_df.dropna().to_csv(train_fname)
                    valid_df.dropna().to_csv(valid_fname)
                    print(f'train fold{j+1}: {train_fname}')
                    print(f'valid fold{j+1}: {valid_fname}')
                    #print(train_df['Event'].sum()/train_df.shape[0])
                    #print(valid_df['Event'].sum()/valid_df.shape[0])


# ---
# # Baseline
# - age + gender

# In[18]:


version = '0'
import_filename = 'final_distance_score.csv'

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
    #all_df = all_df.loc[intersecting_samples].copy()
    clinical_cohort = clinical_cohort.loc[intersecting_samples].copy()
    
    clinical_cohort.rename(columns = {'age_at_initial_pathologic_diagnosis':'age'}, inplace = True)
    clinical_cohort = clinical_cohort[clinical_cohort.age.notnull() & clinical_cohort.gender.notnull()].copy()
    #clinical_cohort['age'] /= max_age
    clinical_cohort['gender'] = clinical_cohort['gender'].map(gender_dict) 
    clinical_cohort.rename(columns = {'gender':'gender=FEMALE'},inplace=True) 
    
    for target in survival_targets:
        print(f'---\nTarget: {target}')
        current_target_col = [target, f'{target}.time', 'age', 'gender=FEMALE']
        clinical_tmp = clinical_cohort[current_target_col].copy()
        #feature_clinical = all_df.merge(clinical_tmp, left_index=True, right_index=True)
        #feature_clinical = feature_clinical[feature_clinical[target].notnull() & feature_clinical[f'{target}.time'].notnull()].copy()
        feature_clinical = clinical_tmp[clinical_tmp[target].notnull() & clinical_tmp[f'{target}.time'].notnull()].copy()
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
            train_valid_fname = os.path.join(cohort_save_dir, f'{target}_train_valid_fold{j+1}.csv')
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


# ---
# # v12.1
# - average opensea beta value of each autosome + age + gender

# In[8]:


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


# ---
# # v12.2
# - average opensea beta value per autosome

# In[9]:


version = '12.2'
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
    
    for target in survival_targets:
        print(f'---\nTarget: {target}')
        current_target_col = [target, f'{target}.time']
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
            test_df.dropna().to_csv(test_fname)
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
                    
                    train_df.dropna().to_csv(train_fname)
                    valid_df.dropna().to_csv(valid_fname)
                    print(f'train fold{j+1}: {train_fname}')
                    print(f'valid fold{j+1}: {valid_fname}')
                    #print(train_df['Event'].sum()/train_df.shape[0])
                    #print(valid_df['Event'].sum()/valid_df.shape[0])


# In[ ]:




