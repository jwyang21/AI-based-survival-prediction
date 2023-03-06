#for version in 9.1 9.2 10.1 10.2 11.2 11.1 0 6.1 6.2 7.1 7.2 8.1 8.2 
for version in 9.2 10.1 10.2 11.2 11.1 0 6.1 6.2 7.1 7.2 8.1 8.2
do
    echo "bash 230301_DeepSurv-all-cohorts-v$version.sh > ../log/230301_DeepSurv-all-cohorts-v$version.log"
    bash 230301_DeepSurv-all-cohorts-v$version.sh > ../log/230301_DeepSurv-all-cohorts-v$version.log
done
