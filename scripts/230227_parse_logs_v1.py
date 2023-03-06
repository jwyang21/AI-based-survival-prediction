#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os


# In[3]:


all_cohorts = ['TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-THCA', 'TCGA-UCEC']


# In[4]:


server = 'bhi-gpu3'


# # Parse logs of v1 ~ v5

# ## version 1

# In[19]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 1

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# ---
# # v1.2

# In[20]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 1.2

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# In[ ]:





# ---

# # v2

# In[22]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 2

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# ---
# # v2.2

# In[23]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 2.2

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# ---

# # v3

# - running in bhi-gpu1 (23/03/03)

# In[24]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 3

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# ---
# # v3.2

# In[29]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 3.2

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# ---

# # v4

# In[26]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 4

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# ---

# # v4.2

# In[27]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 4.2

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# ---

# # v5

# In[28]:


cpg_type = 'opensea'
result_dir = '/data/project/3dith/result/'
version = 5

figure_log_fname = os.path.join(result_dir, f'v{version}_significant_figures.txt')
figure_download_fname = os.path.join(result_dir, f'download_v{version}_significant_figures.sh')
print('figure_log_fname:',figure_log_fname)
print('figure_download_fname:',figure_download_fname)

significant_cohorts = []

with open(figure_log_fname, 'w') as f:
    with open(figure_download_fname, 'w') as f2:
        for cohort in all_cohorts:
            print(f"===\n{cohort}")
            cohort_type = cohort.split('-')[1]
            cohort_dir = os.path.join(result_dir, cpg_type, f'{cohort}', f'v{version}')
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
                    if cohort not in significant_cohorts:
                        significant_cohorts.append(cohort)
                        f2.write(f'mkdir {cohort}')
                        f2.write('\n')
                    display(tmp)
                    print(tmp.fname.values)
                    for i in range(tmp.shape[0]):
                        current_csv_file = tmp.fname.values[i]
                        current_item = current_csv_file.split('.csv')[0]
                        logrank_fname = os.path.join(cohort_dir, f'LogRank_{current_item}.png')
                        nll_fname = os.path.join(cohort_dir, f'{current_item}_train_log_nll.png')
                        c_index_fname = os.path.join(cohort_dir, f'{current_item}_train_log_c-index.png')
                        assert os.path.exists(logrank_fname)
                        assert os.path.exists(nll_fname)
                        assert os.path.exists(c_index_fname)
                        f.write(logrank_fname)
                        f.write('\n')
                        f.write(nll_fname)
                        f.write('\n')
                        f.write(c_index_fname)
                        f.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{logrank_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{nll_fname} ./{cohort}')
                        f2.write('\n')
                        f2.write(f'scp jwyang@{server}.snu.ac.kr:{c_index_fname} ./{cohort}')
                        f2.write('\n')
                else:
                    print(f"All results are insignificant")
f2.close()
f.close()


# In[ ]:




