cpg_type=opensea
result_dir=/data/project/3dith/result
data_dir=/data/project/3dith/data/
version=9.1
script=230301_DeepSurv.py
num_folds=5
lr=0.0005
cohort=TCGA-PRAD
echo "python3 $script --version $version --cpg_type $cpg_type --cohort $cohort --result_dir $result_dir --data_dir $data_dir --num_folds $num_folds --lr $lr"
python3 $script --version $version --cpg_type $cpg_type --cohort $cohort --result_dir $result_dir --data_dir $data_dir --num_folds $num_folds --lr $lr
