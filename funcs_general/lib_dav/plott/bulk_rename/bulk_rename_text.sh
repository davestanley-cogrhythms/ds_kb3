

function rename_plott()
{
	echo find -name "$1" -exec sed -i 's/plot'$2'/plott'$2'/g' '{}' \;
	find -name "$1" -exec sed -i 's/plot'$2'/plott'$2'/g' '{}' \;
}


function rename_all_files()
{
	rename_plott "$1" _all
	rename_plott "$1" _ani_pairs
	rename_plott "$1" _fs
	rename_plott "$1" _matrix3D
	rename_plott "$1" _spect
	rename_plott "$1" _ani
	rename_plott "$1" _coherency
	rename_plott "$1" _handles
	rename_plott "$1" _psd
}

echo Warning: This code	was only tested	in CentOS Linux. The find and sed functions will
echo likely choke on Mac.

rename_all_files '*.m'
rename_all_files '*.html'


