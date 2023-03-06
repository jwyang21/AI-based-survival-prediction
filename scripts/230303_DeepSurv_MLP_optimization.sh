cpg_type=opensea
result_dir=/data/project/3dith/result
data_dir=/data/project/3dith/data/
data_version=11.1
script=230303_DeepSurv_MLP_optimization.py 
num_folds=5

for save_version in 1 1.2 2 2.2 3 3.2 4 4.2 5 5.2
do
	for cohort in TCGA-BLCA TCGA-BRCA TCGA-CHOL TCGA-COAD TCGA-KIRC TCGA-KIRP TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PAAD TCGA-PRAD TCGA-THCA TCGA-UCEC
	do
		echo "python3 $script --save_version $save_version --data_version $data_version --cpg_type $cpg_type --cohort $cohort --result_dir $result_dir --data_dir $data_dir --num_folds $num_folds > ../log/230303_DeepSurv_MLP_optimization_v$save_version.log"
        	python3 $script --save_version $save_version --data_version $data_version --cpg_type $cpg_type --cohort $cohort --result_dir $result_dir --data_dir $data_dir --num_folds $num_folds > ../log/230303_DeepSurv_MLP_optimization_v$save_version.log
	done
done
