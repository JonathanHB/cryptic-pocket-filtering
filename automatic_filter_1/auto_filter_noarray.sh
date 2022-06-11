#create .bsub files and submit them

jobs=9 #number of jobs to run
increment=2500 #spacing between jobs
restart=false #whether to use the starts array below or start at the increments above
starts=(137 2827 5200 7633 10241 13368 15915 18008 20301) #indices at which to start
directory="/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/new_pockets_2"

for i in $(seq 0 $(($jobs-1)))
do
	x=$(($i+1))
	if $restart
	then
		j=${starts[i]}
	else
		j=$(($i*$increment))
	fi

	#note that loading the original filter-submit-noarray.bsub file into a bash variable,
	#processing the variable with sed, and then saving it using printf loses all the returns
	sed "s/index_arg1/${j}/" "${directory}/cryptic-pocket-filtering-2/automatic_filter_1/filter-submit-noarray.bsub" > "${directory}/iofiles_na/submission_pairing_scripts/filter-submit-noarray-${i}.bsub"
	sed -i "s/index_arg2/$(($x*$increment))/" "${directory}/iofiles_na/submission_pairing_scripts/filter-submit-noarray-${i}.bsub"

	bsub < "${directory}/iofiles_na/submission_pairing_scripts/filter-submit-noarray-${i}.bsub"

done
