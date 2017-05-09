#/bin/bash

script=/scratch/scripts/chunk.R
out_loc=s3://cory-dbtest/cer-trn/
start=1
end=15225
int=500
count=1

for inc in $(eval  echo {$start..$end..$int})
#for inc in {1..15225..500}
#for inc in {1..20..20}
do
	first=$(($inc-1))
	
	if [ $first -gt 10 ]; then
		if [ $(($first + $int)) -gt $end ]; then
			echo "R -f $script "$inc" "$end" temp.$count.rds"
		else
			echo "R -f $script "$inc" "$(($first + $int))" temp.$count.rds"
		fi
	else
	echo "R -f $script "$start" "$int" temp.$count.rds"
	fi
	
	echo "aws s3 cp /scratch/scripts/temp.$count.rds $out_loc"
	echo -e "\n"
	count=$(($count + 1))
done
