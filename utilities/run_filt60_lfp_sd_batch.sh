



#for i in {5,19,29,40}
for i in {1..79}
do

	# Added runtime submission directive
	qsub -l h_rt=576:00:00 single_node_batch "filter60_lfp_sample_delay(${i})" localOutput

done

