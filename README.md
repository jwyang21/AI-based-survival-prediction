# AI-based-survival-prediction

## 1. Introduction    
This repository contains codes for an AI-based survival prediction, utilizing stem closeness and its related features as input.         
For stem closeness, please refer to [stem-closeness](https://github.com/jwyang21/stem-closeness).     
For the AI model, modified codes from [DeepSurv](https://github.com/jaredleekatzman/DeepSurv) was used.     

## 2. Installation
```python
conda env create -f surv-prediction.yaml
```

## 3. Preprocess
Before running the model, files containing input features should be made first.     
```python
python3 230223_concat_autosome_pc1-8.py
python3 230223_concat_autosome_pc1-ALL-avg10_20.py
python3 230223_concat_autosome_pc1-ALL.py
python3 230301_concat_10_20_averaged_pc1.py
python3 230228_concat_distances_score.py
python3 230223_merge_pc1_clinical.py
python3 230306_compute_avg_opensea_beta_per_autosome.py
python3 230301_prepare_data.py
```

## 4. Run model
```shell
bash 230301_DeepSurv_v6_to_v11.sh
bash 230303_DeepSurv_MLP_optimization.sh
```
