#create input files and submit array job
for i in $(seq 0 9)
do
	x=$(($i+1))
	printf "$(($i*2500))\n$(($x*2500))\n" > "filter_in.${x}"

done

bsub -J "job_Array[1-10]" -i "filter_in.%I" -o "filter_out.%I" -e "filter_out.%I" < filter-submit-array.bsub
