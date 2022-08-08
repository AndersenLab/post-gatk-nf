#!/bin/bash
while getopts f:n: flag
do
    case "${flag}" in
        f) file=${OPTARG};;
        n) num_out=${OPTARG};;
    esac
done

@num_out=20

#Must keep double quotes so that shell will expand the variable
cat ${file} | sed -e "s/numoutlieriter:  ../numoutlieriter:  ${num_out}/g" > outlier_eigpar 