#!/bin/bash

#where you want to look at the data
analysis_dir=10.12.22_tree_analysis

#Dir with tree file folders
#file_dir=

#list file option #1
files=$(ls | grep sws_tree | xargs)

for file in ${files}
do  
    cp -f ${file}/tree/eiganstrat_input.min4.tree ${analysis_dir}/${file}.min4.tree

done