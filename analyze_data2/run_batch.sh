# Run this file. run_batch2.sh is dummy file for permutation case

floor()
{
echo ${1/.*}
}

	
function run_qsub {

mode=$1
stage=$2

#mode_prefix=`perl -w -e "use POSIX; print floor($mode), qq{\n}"` 
mode_firstdecimal=`perl -w -e "use POSIX; print floor(${mode}*10) - 10*floor($mode), qq{\n}"`
mode_sixthdecimal=`perl -w -e "use POSIX; print floor(${mode}*1000000) - 10*floor(${mode}*100000), qq{\n}"`
#perl -w -e "use POSIX; print floor(${mode}*1000000) - 10*floor($mode*100000), qq{\n}"
#perl -w -e "use POSIX; print floor(${mode}*1000000.0) , qq{\n}"
#perl -w -e "use POSIX; print  10*floor(${mode}*100000.0), qq{\n}"
#perl -w -e "use POSIX; print  floor(${mode}*100000.0), qq{\n}"
#perl -w -e "use POSIX; print  ${mode}*100000.0, qq{\n}"
#perl -w -e "use POSIX; print  ${mode}, qq{\n}"


temp=`echo "$mode * 10^6" | bc -l`
temp=$(floor $temp 0)
temp2=`echo "$mode * 10^5" | bc -l`
temp2=$(floor $temp2 0)
temp2=`echo "$temp2 * 10" | bc -l`
mode_sixthdecimal=`echo "$temp - $temp2" | bc -l`



echo mode_sixthdecimal $mode_sixthdecimal

source Nunits.txt
source Nunits_field.txt
source Nfield_field.txt
source Nunits_units.txt
source Nlfp.txt

#for fi in {74,75,76}
#for fi in {1..79}
for fi in {1..2}
do
	#for i in {1..11}
	for i in 1
	do
		if [ $mode_firstdecimal -le 2 ];
		then
			curr_Ncells=${Nunits_var[${fi}-1]}
		elif [ $mode_firstdecimal -eq 3 ];
		then
			curr_Ncells=${Nunits_field_var[${fi}-1]}
		elif [ $mode_firstdecimal -eq 4 ];
                then
                        curr_Ncells=${Nfield_field_var[${fi}-1]}
		elif [ $mode_firstdecimal -eq 5 ];
                then
                        curr_Ncells=${Nunits_units_var[${fi}-1]}
		elif [ $mode_firstdecimal -eq 6 ];
                then
                        curr_Ncells=${Nlfp_var[${fi}-1]}
		fi

		echo Curr Ncells = $curr_Ncells


		for j in `seq 1 $curr_Ncells`
		#for j in `seq $curr_Ncells $curr_Ncells`
		do
			if [ $mode_sixthdecimal -lt 1 ]
			then
	#			qsub -l h_rt=576:00:00 matlab_single_node_batch.sh "run_analysis_wrapp(${mode},${stage},${fi},1:11,${j},1)" localOutput
	#			qsub -l h_rt=0:30:00 matlab_single_node_batch.sh "run_analysis_wrapp(${mode},${stage},${fi},[1:12],${j},1)" localOutput
				qsub -l h_rt=0:30:00 run_standalone_job.sh run_analysis_binwrapp ${mode} ${stage} ${fi} "[1:12]" ${j} 1

			else
	#			qsub -l h_rt=24:00:00 matlab_single_node_batch.sh "run_analysis_wrapp(${mode},${stage},${fi},[1 2 5 6],${j},1)" localOutput
				qsub -l h_rt=4:00:00 run_standalone_job.sh run_analysis_binwrapp ${mode} ${stage} ${fi} "[1:12]" ${j} 1 	#For permutation or bootstrap test
	#			qsub -pe omp 3 run_standalone_job.sh run_analysis_binwrapp ${mode} ${stage} ${fi} 9:10 ${j} 1 	#For multithreaded. Shouldn't need this now since binary is forced to be single threaded
	#			./run_standalone_job.sh ${mode} ${stage} ${fi} "[1 2 5 6]" ${j} 1 
			fi


		done
	done
done

}


function run_allstages {

#matlab -nodisplay -nosplash -r "calc_ncells; exit" >! nfiles_out
#matlab -nodisplay -nosplash -r "calc_nlfp; exit" >! nfiles_out

mode=$1
#run_qsub $mode -2
#run_qsub $mode -1
#run_qsub $mode 0
run_qsub $mode 2
#run_qsub $mode 3
#run_qsub $mode 6

}


#run_allstages 1.2

#run_allstages 2.0
#run_allstages 2.21
#run_allstages 2.00001
#run_allstages 2.00011
#run_allstages 2.20001
#run_allstages 2.20011
#run_allstages 2.201
#run_allstages 2.20101
#run_allstages 2.20111
#run_allstages 2.21101
#run_allstages 2.2012
#run_allstages 2.202

#run_allstages 2.201112
#run_allstages 2.201511
#run_allstages 22.401311

#run_allstages 2.201310
#run_allstages 2.201410
#run_allstages 2.201510
#run_allstages 2.201110


#run_allstages 2.301410
#run_allstages 22.401310
#run_allstages 32.501110

#run_qsub 3.2 5
#run_qsub 3.201 5
#run_qsub 3.20001 5
#run_qsub 3.20101 5

#run_allstages 5
#run_qsub 5.1 5

#run_qsub 8 5
#run_qsub 9 5

#run_qsub 8.00001 5
#run_qsub 9.00001 5

#run_allstages 12.2
#run_allstages 15.2

#run_allstages 40.601310


#run_allstages 43.601311
#run_allstages 43.601210

#run_allstages 42.601310
#run_allstages 44.601311

#run_allstages 42.601210
#run_allstages 44.601210


#run_allstages 41.601310
#run_allstages 22.401310
#run_allstages 2.221510


#run_allstages 41.611310
#run_allstages 41.641311

#run_allstages 22.411310
#run_allstages 22.441311

#run_allstages 22.4113100
#run_allstages 22.4113101

#run_allstages 22.4313110

#run_allstages 2.2415100
#run_allstages 2.2415101

#run_allstages 41.6213101
#run_allstages 41.6313101
#run_allstages 41.6413101


#run_allstages 2.2511101
#run_allstages 2.2511111
#run_allstages 41.6513101
#run_allstages 41.6513111
#run_allstages 22.4513101
#run_allstages 22.4513111

#run_allstages 2.2611101
#run_allstages 2.2611111
#run_allstages 41.6613101
#run_allstages 41.6613111
#run_allstages 22.4613101
run_allstages 41.6014101
#run_allstages 22.4615111



#run_allstages 41.6213111
#run_allstages 22.4213111

#run_allstages 41.6313111
#run_allstages 41.6413111






