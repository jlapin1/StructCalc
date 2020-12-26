#/bin/bash

if [ $# -lt 1 ]
then
	dir="untitled_project"
else
	dir=$1
fi

mkdir $dir
cp ~/structcalc/prgms/src/config.conf $dir/
cd $dir
mkdir raw input output analysis

printf "\nCreated new project in directory: $dir/\n\n" 
