#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import glob


# In[2]:


all_lrs = [0.0005, 0.0003, 0.0001, 0.00007, 0.00005, 0.00003, 0.00001]
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']


# In[3]:


d = '/data/project/3dith/result/opensea/TCGA-BLCA/v9.1'
files = glob.glob(os.path.join(d, f'all_versions*'))


# In[4]:


files


# In[5]:


server = 'bhi-gpu3'


# # parse logs of v6 ~ 
# ---

# ## version 9.1

# In[6]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 9.1

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print(f'significant figure files log: {figure_log_fname}')
print(f'significant figure files download: {figure_download_fname}')

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            for lr in all_lrs:
                print(f'---\nlr: {lr}')

                cohort_type = cohort.split('-')[1]
                cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
                log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
                #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
                if len(log_list) ==1:
                    cohort_log = log_list[0]

                    globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
                    globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
                    globals()[f'{cohort_type}'].rename({"time":"time (sec)"}, axis = 'columns', inplace = True)
                    tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
                    tmp = tmp[tmp['train_c']>=0.5].copy()
                    tmp = tmp[tmp['valid_c']>=0.5].copy()
                    tmp = tmp[tmp['test_c']>=0.5].copy()
                    if tmp.shape[0] >= 1:
                        if cohort not in significant_cohorts:
                            significant_cohorts.append(cohort)
                            f2.write(f'mkdir {cohort}')
                            f2.write('\n')
                        display(tmp)
                        #print(tmp.fname.values)
                        for i in range(tmp.shape[0]):
                            print(tmp.fname.values[i])
                            print(tmp.train_valid_nll.values[i])
                            print(tmp.train_valid_c_plot.values[i])
                            print(tmp.LogRank_plot.values[i])
                            f.write(tmp.train_valid_nll.values[i])
                            f.write('\n')
                            f.write(tmp.train_valid_c_plot.values[i])
                            f.write('\n')
                            f.write(tmp.LogRank_plot.values[i])
                            f.write('\n')
                            f2.write(f'scp jwyang@{server}.snu.ac.kr:{tmp.train_valid_nll.values[i]} ./{cohort}')
                            f2.write('\n')
                            f2.write(f'scp jwyang@{server}.snu.ac.kr:{tmp.train_valid_c_plot.values[i]} ./{cohort}')
                            f2.write('\n')
                            f2.write(f'scp jwyang@{server}.snu.ac.kr:{tmp.LogRank_plot.values[i]} ./{cohort}')
                            f2.write('\n')
                            
                    else:
                        print(f"All results are insignificant")
                elif len(log_list) == 0:
                    print(f"{cohort} has no log file for lr {lr}")
                else:
                    raise ValueError
f.close()
f2.close()


# ---
# # v9.2

# In[7]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 9.2
figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print(f'significant figure files log: {figure_log_fname}')
print(f'significant figure files download: {figure_download_fname}')

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            for lr in all_lrs:
                print(f'---\nlr: {lr}')

                cohort_type = cohort.split('-')[1]
                cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
                log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result_lr_{lr}.csv'))
                #log_list = glob.glob(os.path.join(cohort_dir, f'all_versions_deepsurv_logrank_result.csv'))
                if len(log_list) ==1:
                    cohort_log = log_list[0]

                    globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
                    globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
                    globals()[f'{cohort_type}'].rename({"time":"time (sec)"}, axis = 'columns', inplace = True)
                    tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
                    tmp = tmp[tmp['train_c']>=0.5].copy()
                    tmp = tmp[tmp['valid_c']>=0.5].copy()
                    tmp = tmp[tmp['test_c']>=0.5].copy()
                    if tmp.shape[0] >= 1:
                        if cohort not in significant_cohorts:
                            significant_cohorts.append(cohort)
                            f2.write(f'mkdir {cohort}')
                            f2.write('\n')
                        display(tmp)
                        #print(tmp.fname.values)
                        for i in range(tmp.shape[0]):
                            print(tmp.fname.values[i])
                            print(tmp.train_valid_nll.values[i])
                            print(tmp.train_valid_c_plot.values[i])
                            print(tmp.LogRank_plot.values[i])
                            f.write(tmp.train_valid_nll.values[i])
                            f.write('\n')
                            f.write(tmp.train_valid_c_plot.values[i])
                            f.write('\n')
                            f.write(tmp.LogRank_plot.values[i])
                            f.write('\n')
                            f2.write(f'scp jwyang@{server}.snu.ac.kr:{tmp.train_valid_nll.values[i]} ./{cohort}')
                            f2.write('\n')
                            f2.write(f'scp jwyang@{server}.snu.ac.kr:{tmp.train_valid_c_plot.values[i]} ./{cohort}')
                            f2.write('\n')
                            f2.write(f'scp jwyang@{server}.snu.ac.kr:{tmp.LogRank_plot.values[i]} ./{cohort}')
                            f2.write('\n')
                            
                    else:
                        print(f"All results are insignificant")
                elif len(log_list) == 0:
                    print(f"{cohort} has no log file for lr {lr}")
                else:
                    raise ValueError
f.close()
f2.close()


# In[ ]:





# ---

# In[8]:


'''
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 1.2
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']
for cohort in all_cohorts:
    print(f"===\n{cohort}")
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}_v{version}')
    cohort_log = os.path.join(cohort_dir, 'all_versions_deepsurv_logrank_result.csv')
    if not os.path.exists(cohort_log):
        print(f"{cohort} has no log file")
    else:
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
        globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
        tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
        tmp = tmp[tmp['train_c']>=0.5].copy()
        tmp = tmp[tmp['valid_c']>=0.5].copy()
        tmp = tmp[tmp['test_c']>=0.5].copy()
        if tmp.shape[0] >= 1:
            display(tmp)
            print(tmp.fname.values)
        else:
            print(f"All results are insignificant")
'''


# ---

# In[9]:


'''
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 2
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']
for cohort in all_cohorts:
    print(f"===\n{cohort}")
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}_v{version}')
    cohort_log = os.path.join(cohort_dir, 'all_versions_deepsurv_logrank_result.csv')
    if not os.path.exists(cohort_log):
        print(f"{cohort} has no log file")
    else:
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
        globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
        tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
        tmp = tmp[tmp['train_c']>=0.5].copy()
        tmp = tmp[tmp['valid_c']>=0.5].copy()
        tmp = tmp[tmp['test_c']>=0.5].copy()
        if tmp.shape[0] >= 1:
            display(tmp)
            print(tmp.fname.values)
        else:
            print(f"All results are insignificant")
'''


# ---

# In[10]:


'''
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 2.2
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']
for cohort in all_cohorts:
    print(f"===\n{cohort}")
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}_v{version}')
    cohort_log = os.path.join(cohort_dir, 'all_versions_deepsurv_logrank_result.csv')
    if not os.path.exists(cohort_log):
        print(f"{cohort} has no log file")
    else:
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
        globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
        tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
        tmp = tmp[tmp['train_c']>=0.5].copy()
        tmp = tmp[tmp['valid_c']>=0.5].copy()
        tmp = tmp[tmp['test_c']>=0.5].copy()
        if tmp.shape[0] >= 1:
            display(tmp)
            print(tmp.fname.values)
        else:
            print(f"All results are insignificant")
'''


# ---

# In[11]:


'''
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 3
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']
for cohort in all_cohorts:
    print(f"===\n{cohort}")
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}_v{version}')
    cohort_log = os.path.join(cohort_dir, 'all_versions_deepsurv_logrank_result.csv')
    if not os.path.exists(cohort_log):
        print(f"{cohort} has no log file")
    else:
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
        globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
        tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
        tmp = tmp[tmp['train_c']>=0.5].copy()
        tmp = tmp[tmp['valid_c']>=0.5].copy()
        tmp = tmp[tmp['test_c']>=0.5].copy()
        if tmp.shape[0] >= 1:
            display(tmp)
            print(tmp.fname.values)
        else:
            print(f"All results are insignificant")
'''


# ---

# In[12]:


'''
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 3.2
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']
for cohort in all_cohorts:
    print(f"===\n{cohort}")
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}_v{version}')
    cohort_log = os.path.join(cohort_dir, 'all_versions_deepsurv_logrank_result.csv')
    if not os.path.exists(cohort_log):
        print(f"{cohort} has no log file")
    else:
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
        globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
        tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
        tmp = tmp[tmp['train_c']>=0.5].copy()
        tmp = tmp[tmp['valid_c']>=0.5].copy()
        tmp = tmp[tmp['test_c']>=0.5].copy()
        if tmp.shape[0] >= 1:
            display(tmp)
            print(tmp.fname.values)
        else:
            print(f"All results are insignificant")
'''


# ---

# In[13]:


'''
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 4
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']
for cohort in all_cohorts:
    print(f"===\n{cohort}")
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}_v{version}')
    cohort_log = os.path.join(cohort_dir, 'all_versions_deepsurv_logrank_result.csv')
    if not os.path.exists(cohort_log):
        print(f"{cohort} has no log file")
    else:
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
        globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
        tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
        tmp = tmp[tmp['train_c']>=0.5].copy()
        tmp = tmp[tmp['valid_c']>=0.5].copy()
        tmp = tmp[tmp['test_c']>=0.5].copy()
        if tmp.shape[0] >= 1:
            display(tmp)
            print(tmp.fname.values)
        else:
            print(f"All results are insignificant")
'''


# ---

# In[14]:


'''
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 4.2
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']
for cohort in all_cohorts:
    print(f"===\n{cohort}")
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}_v{version}')
    cohort_log = os.path.join(cohort_dir, 'all_versions_deepsurv_logrank_result.csv')
    if not os.path.exists(cohort_log):
        print(f"{cohort} has no log file")
    else:
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
        globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
        tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
        tmp = tmp[tmp['train_c']>=0.5].copy()
        tmp = tmp[tmp['valid_c']>=0.5].copy()
        tmp = tmp[tmp['test_c']>=0.5].copy()
        if tmp.shape[0] >= 1:
            display(tmp)
            print(tmp.fname.values)
        else:
            print(f"All results are insignificant")
'''


# ---

# In[15]:


'''
cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 5
all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']
for cohort in all_cohorts:
    print(f"===\n{cohort}")
    cohort_type = cohort.split('-')[1]
    cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}_v{version}')
    cohort_log = os.path.join(cohort_dir, 'all_versions_deepsurv_logrank_result.csv')
    if not os.path.exists(cohort_log):
        print(f"{cohort} has no log file")
    else:
        globals()[f'{cohort_type}'] = pd.read_csv(cohort_log)
        globals()[f'{cohort_type}'].rename({"Unnamed: 0":"fname"}, axis = 'columns', inplace = True)
        tmp = globals()[f'{cohort_type}'][globals()[f'{cohort_type}']['LogRank_p']<=5e-2].copy()
        tmp = tmp[tmp['train_c']>=0.5].copy()
        tmp = tmp[tmp['valid_c']>=0.5].copy()
        tmp = tmp[tmp['test_c']>=0.5].copy()
        if tmp.shape[0] >= 1:
            display(tmp)
            print(tmp.fname.values)
        else:
            print(f"All results are insignificant")
'''


# In[ ]:




