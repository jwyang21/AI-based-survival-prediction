import pandas as pd
import numpy as np
import os

def standardize(v):
    return (v-np.mean(v))/np.std(v)

cohort_fname = '/data/project/3dith/data/cohorts.manifest'
cohorts = pd.read_csv(cohort_fname).cohort.values
tcga_cohorts = []
for c in cohorts:
    if 'TCGA' in c:
        tcga_cohorts.append(c)
        
data_dir = '/data/project/3dith/data/'
samplenames_fname = '/data/project/3dith/data/samplenames.npz'
AUT_LIST = [f'chr{i}' for i in np.arange(1, 23)]
all_cpg_types = ['opensea','island', 'shelf_shore']

clinical_fname = '/data/project/3dith/data/TCGA-CDR-SupplementalTableS1.xlsx'
if not os.path.exists(clinical_fname):
    cmd_ = f'wget https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81 -O {clinical_fname}'
    os.system(cmd_)
    print(f"TCGA clinical data is downloaded to {clinical_fname}")
clinical = pd.read_excel(clinical_fname, sheet_name = 'TCGA-CDR', index_col = 0)


clinical_all_columns = ['bcr_patient_barcode', 'type', 'age_at_initial_pathologic_diagnosis', 'gender', 'race', 'ajcc_pathologic_tumor_stage', 'clinical_stage', \
                        'histological_type', 'histological_grade', 'initial_pathologic_dx_year', 'menopause_status', 'birth_days_to', 'vital_status', 'tumor_status', \
                        'last_contact_days_to', 'death_days_to', 'cause_of_death', 'new_tumor_event_type', 'new_tumor_event_site', 'new_tumor_event_site_other', 'new_tumor_event_dx_days_to', \
                        'treatment_outcome_first_course', 'margin_status', 'residual_tumor', 'OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time', 'Redaction']

clinical_target_columns = ['OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time']


