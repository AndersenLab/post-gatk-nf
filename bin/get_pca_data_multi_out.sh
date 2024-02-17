#!/bin/bash

analysis_dir=/projects/b1059/projects/Ryan/telomere_variation/cb_optimal_mapping/data/10.14.22_sws_tree_analysis
mkdir $analysis_dir


ce_w=/projects/b1059/projects/Ryan/post-gatk-nf/10.6.22_ce_singleton_sws
ce_wo=/projects/b1059/projects/Ryan/post-gatk-nf/10.6.22_ce_wo_singleton_sws

cb_w=/projects/b1059/projects/Ryan/post-gatk-nf/10.6.22_cb_singleton_sws
cb_wo=/projects/b1059/projects/Ryan/post-gatk-nf/10.6.22_cb_wo_singleton_sws

ct_w=/projects/b1059/projects/Ryan/post-gatk-nf/10.6.22_ct_singleton_sws
ct_wo=/projects/b1059/projects/Ryan/post-gatk-nf/10.6.22_ct_wo_singleton_sws

#for dir in ${ce_dir} ${cb_dir} ${ct_dir}
for dir in ${ce_w} ${ce_wo} ${cb_w} ${cb_wo} ${ct_w} ${ct_wo}
do 
#Get the species name
sp=$(basename ${dir})

cd ${analysis_dir}

mkdir ${sp}

cd ${dir}/EIGESTRAT

#Create a variable for each output dir for LD params
lds=$(ls) 

    for ld in ${lds}
    do

    #Pull without outlier removeal
    mkdir ${analysis_dir}/${sp}/${ld}

    #Pull the vcf into the analysis dir
    #cp ${ld}/INPUTFILES/*.vcf.gz

    cp ${ld}/NO_REMOVAL/*.tsv ${analysis_dir}/${sp}/${ld}/
    cp ${ld}/NO_REMOVAL/*.evac ${analysis_dir}/${sp}/${ld}/
    
    cd ${ld}/OUTLIER_REMOVAL
    out_its=$(ls)
        for it in ${out_its}
        do
        #Make an output folder in the analysis directory
        cd ${analysis_dir}/${sp}/${ld}
        mkdir ${it}
        
        cd ${dir}/EIGESTRAT
        #With outlier removal
        cp ${ld}/OUTLIER_REMOVAL/${it}/*.tsv ${analysis_dir}/${sp}/${ld}/${it}/
        cp ${ld}/OUTLIER_REMOVAL/${it}/*.evac ${analysis_dir}/${sp}/${ld}/${it}/
        done

    done

cd ${analysis_dir}

done
