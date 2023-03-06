cpg_type=opensea
result_dir=/data/project/3dith/result
data_dir=/data/project/3dith/data/
version=9.2
script=230301_DeepSurv.py
num_folds=5

for lr in 0.0005 0.0003 0.0001 0.00007 0.00005 0.00003 0.00001
do  
    for cohort in TCGA-BLCA TCGA-BRCA TCGA-CHOL TCGA-COAD TCGA-KIRC TCGA-KIRP TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PAAD TCGA-PRAD TCGA-THCA TCGA-UCEC
    do
        echo "python3 $script --version $version --cpg_type $cpg_type --cohort $cohort --result_dir $result_dir --data_dir $data_dir --num_folds $num_folds --lr $lr"
        python3 $script --version $version --cpg_type $cpg_type --cohort $cohort --result_dir $result_dir --data_dir $data_dir --num_folds $num_folds --lr $lr
    done
done