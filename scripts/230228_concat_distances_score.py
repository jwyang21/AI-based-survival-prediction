#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os


# In[2]:


scores = pd.read_csv('/data/project/3dith/data/cohort-1-best-score-km.csv', index_col = 0)


# In[14]:


scores.head(3)


# In[15]:


cpg_type = 'opensea'
cohort = 'TCGA-BLCA'


# In[16]:


score_fname = scores.loc[cohort][f'filename_{cpg_type}']


# In[17]:


score_fname


# In[19]:


pd.read_csv(score_fname, index_col = 0)


# In[12]:


cohort_pc1_dir =


# In[13]:


pd.read_csv(f'/data/project/3dith/pipelines/{cpg_type}-pipeline/1_compute-score-{cpg_type}/result/{cohort}/normal-distance_cosine-sim_bdm_pc1-avg_simple-avg_half_standardized_to-chrom-8.csv')


# In[ ]:




