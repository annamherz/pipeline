#! /bin/bash

slurmlog=$2
slurm_dir=$1

echo checking $slurm_logs in $slurm_dir

file_base="$slurm_dir/prod_$slurmlog\_*"


for f in $file_base; do

declare regex="\s+STOP PMEMD Terminated Abnormally!\s+"
declare file_content=$( cat "${f}" )
if [[ " $file_content " =~ $regex ]] # please note the space before and after the file content
    then
        echo "$f is a failed run, pmemd terminated"
fi

declare regex="\s+unsigned long long of length = 42Failed an illegal memory access was encountered\s+"
declare file_content=$( cat "${f}" )
if [[ " $file_content " =~ $regex ]] # please note the space before and after the file content
    then
	echo "$f is a failed run, illegal memory access"
	   
	free=$( sed -n 4p $f )
	bound=$( sed -n 10p $f )
	for leg in "$free" "$bound" ; do

	while IFS=' ' read -ra addr; do
	folder=${addr[2]}
	done <<< "$leg"
	
	dirs=( lambda_0.0000 lambda_0.0909 lambda_0.1818 lambda_0.2727 lambda_0.3636 lambda_0.4546 lambda_0.5454 lambda_0.6364 lambda_0.7273 lambda_0.8182 lambda_0.9091 lambda_1.0000 )
	index=$( echo $f | tail -c 1 )
	amber_file="$folder/${dirs[index]}/amber.out"
	#echo $leg is $( grep -F "CUDA_VISIBLE_DEVICES" $amber_file )

	done
fi	

done
