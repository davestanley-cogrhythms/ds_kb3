echo Warning: This code was only tested in CentOS Linux

for f in  `ls plot_*`
do
	fnew=$(echo $f | sed s/plot_/plott_/)
	echo git mv $f $fnew
	git mv $f $fnew

done

