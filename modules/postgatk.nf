/* 
    ==================
    Subset isotype reference strains
    ==================
*/
process subset_iso_ref_strains {

    label 'postgatk'

    // conda "/data/eande106/software/conda_envs/popgen-nf_env"

    memory 20.GB
    cpus 4

    publishDir "${params.output}/variation", mode: 'copy'
    publishDir "${params.output}/NemaScan", pattern: "div_isotype_list.txt", mode: 'copy'

    input: 
        tuple file(vcf), file(vcf_index), file("ref_strain_isotype.tsv")

    output: 
        tuple file("*.isotype.vcf.gz"), file("*.isotype.vcf.gz.tbi"), emit: vcf
        path 'div_isotype_list.txt', emit: pop_strains


    """
    output=`echo $vcf | sed 's/.vcf.gz/.isotype.vcf.gz/'`
    cut -f1 ${params.sample_sheet} > strain_list.txt

    bcftools view -S strain_list.txt -O u ${vcf} | \\
      bcftools view -O v --min-af 0.000001 --max-af 0.999999 | \\
      vcffixup - | \\
      bcftools view --threads ${task.cpus} -O z > WI.ref_strain.vcf.gz

    bcftools index WI.ref_strain.vcf.gz
    bcftools index --tbi WI.ref_strain.vcf.gz

    # If c.e., convert to isotype names
    if [ ${params.species} = "c_elegans" ]
    then
        bcftools reheader -s ref_strain_isotype.tsv WI.ref_strain.vcf.gz | bcftools view -O z > \${output}
        bcftools index --tbi \${output}
    else
        mv WI.ref_strain.vcf.gz \${output}
        mv WI.ref_strain.vcf.gz.tbi \${output}.tbi
    fi

    # output list of strains for divergent
    bcftools query -l \${output} > div_isotype_list.txt

    """

}

// i know there has to be a better way to do this, but this should work. subset iso ref strains for soft filter vcf
process subset_iso_ref_soft {

    label 'postgatk'

    // conda "/data/eande106/software/conda_envs/popgen-nf_env"

    memory 20.GB
    cpus 4

    publishDir "${params.output}/variation", mode: 'copy'

    input: 
        tuple file(vcf), file(vcf_index), file("ref_strain_isotype.tsv")

    output: 
        tuple file("*.isotype.vcf.gz"), file("*.isotype.vcf.gz.tbi")


    """
    output=`echo $vcf | sed 's/.vcf.gz/.isotype.vcf.gz/'`
    cut -f1 ${params.sample_sheet} > strain_list.txt

    bcftools view -S strain_list.txt -O u ${vcf} | \\
      bcftools view -O v --min-af 0.000001 --max-af 0.999999 | \\
      vcffixup - | \\
      bcftools view --threads ${task.cpus} -O z > WI.ref_strain.vcf.gz

    bcftools index WI.ref_strain.vcf.gz
    bcftools index --tbi WI.ref_strain.vcf.gz

    # If c.e., convert to isotype names
    if [ ${params.species} = "c_elegans" ]
    then
        bcftools reheader -s ref_strain_isotype.tsv WI.ref_strain.vcf.gz | bcftools view -O z > \${output}
        bcftools index --tbi \${output}
    else
        mv WI.ref_strain.vcf.gz \${output}
        mv WI.ref_strain.vcf.gz.tbi \${output}.tbi
    fi

    """

}

process subset_snv {

    label 'postgatk'

    // conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    publishDir "${params.output}/variation", mode: 'copy'

    input:
        tuple file(hardvcf), file(hardvcf_index)

    output:
        tuple file("WI.*.hard-filter.isotype.SNV.vcf.gz"), file("WI.*.hard-filter.isotype.SNV.vcf.gz.tbi")

    """
    output=`echo ${hardvcf} | sed 's/.isotype.vcf.gz/.isotype.SNV.vcf.gz/'`
    bcftools view -O u ${hardvcf} | \
    bcftools view -O v --types snps --min-af 0.000001 --max-af 0.999999 | \
    vcffixup - | \
    bcftools view --threads=3 -O z > \$output
    
    bcftools index --tbi \$output

    """

}

// make a small genotype only vcf for download from cendr for nemascan
process make_small_vcf {

    label 'postgatk'

    // conda "/data/eande106/software/conda_envs/popgen-nf_env"
    publishDir "${params.output}/variation", mode: 'copy'

    input:
        tuple file(vcf), file(vcf_index)

    output:
        tuple file("WI.*.small.hard-filter.isotype.vcf.gz"), file("WI.*.small.hard-filter.isotype.vcf.gz.tbi")

"""
    output=`echo ${vcf} | sed 's/.hard-filter/.small.hard-filter/'`

    bcftools annotate -x INFO,^FORMAT/GT ${vcf} | bgzip > \$output
    bcftools index -t \$output
"""

}


/* 
    ==================
    Tree (only uses SNVs)
    ==================
*/

process convert_tree {

    // label 'tree'
    // container 'shub://bioconvert/bioconvert:latest'

     conda "/data/eande106/software/conda_envs/popgen-nf_env"

    memory { 24.GB + 10.GB * task.attempt }
    errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }

    input:
        tuple file(vcf), file(vcf_index)

    output:
        file("*.stockholm")

"""

    output_phylip=`echo ${vcf} | sed 's/.vcf.gz/.min4.phy/'`

    output_stockholm=`echo \${output_phylip} | sed 's/.phy/.stockholm/'`

    python ${workflow.projectDir}/bin/vcf2phylip.py -i ${vcf}

    bioconvert phylip2stockholm \${output_phylip} \${output_stockholm}

"""

}

process quick_tree {

    label 'tree'

    memory { 44.GB + 10.GB * task.attempt }
    errorStrategy { task.attempt < 5 ? 'retry' : 'ignore' }
    publishDir "${params.output}/tree", mode: 'copy'

    input:
        file(stockholm)

    output:
        file("*.tree")

    """
    output_tree=`echo ${stockholm} | sed 's/.stockholm/.tree/'`

    quicktree -in a -out t ${stockholm} > \${output_tree}

    """

}


process plot_tree {

    label 'R'
    memory { 24.GB + 10.GB * task.attempt }

    // conda "/data/eande106/software/conda_envs/popgen-nf-r_env"

    publishDir "${params.output}/tree", mode: 'copy'

    input:
        file(tree)

    output:
        file("*.pdf")

"""
    Rscript --vanilla ${workflow.projectDir}/bin/plot_tree.R ${tree}
"""



}

/* 
    ==================
    Haplotype and sweep (use both SNVs and INDELs)
    ==================
*/

process haplotype_sweep_IBD {

    label 'postgatk'

    // conda "/data/eande106/software/conda_envs/popgen-nf_env"

    publishDir "${params.output}/haplotype", mode: 'copy'

    memory 17.GB
    cpus 4

    input:
        tuple file(vcf), file(vcf_index), val(contig)

    output: file("*.ibd")

    """
    java -Xmx16G -jar ${workflow.projectDir}/bin/ibdseq.r1206.jar gt=${vcf} minalleles=4 r2max=0.3 ibdlod=3 r2window=1500 nthreads=4 chrom=${contig} out=${contig}


    """

}



process haplotype_sweep_plot {

    label 'R'

    // conda "/data/eande106/software/conda_envs/popgen-nf-r_env"

    memory { 20.GB + 20.GB * task.attempt }
    errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }

    publishDir "${params.output}/haplotype", mode: 'copy'
    publishDir "${params.output}/NemaScan", pattern: "haplotype_df_isotype.bed", mode: 'copy'

    input: 
        path("*")

    output: 
        path("*")

    """

    cat *.ibd > haplotype.tsv

    Rscript --vanilla ${workflow.projectDir}/bin/process_ibd_nf_final.R haplotype.tsv

    """
}


/* 
    ==================
    Divergent regions
    ==================
*/


process prep_variant_coverage {

    label 'postgatk'

    // conda "/data/eande106/software/conda_envs/popgen-nf_env"

    // publishDir "${params.output}/divergent_regions/Mask_DF", mode: 'copy'

    input:
        tuple path(vcf), path(index), val(strain), val(bam), path("ref_strain_isotype.tsv")
        // tuple val(strain), path(bam), path("ref_strain_isotype.tsv")

    output:
        file("*.csv")
        // val(isotype_name)

    """
    # if c.e and the strain is one of those that isotype name differs from strain name, convert isotype name to isotype ref name
    if [ ${params.species} = "c_elegans" ] && grep -q ${strain} ref_strain_isotype.tsv
    then
        ref_iso_pair=`grep ${strain} ref_strain_isotype.tsv`
        ref_iso_pair_array=(\$ref_iso_pair)
        ref_strain_name=\${ref_iso_pair_array[0]}
        iso_name=\${ref_iso_pair_array[1]}
        isotype_name=`echo ${strain} | sed "s/\$ref_strain_name/\$iso_name/"`
    else
        isotype_name=${strain}
    fi

    echo strain,bam,bai > \$isotype_name.csv
    echo \$isotype_name,${params.bam_folder}/${strain}.bam,${params.bam_folder}/${strain}.bam.bai >> \$isotype_name.csv
    """
}

process count_variant_coverage {

    label 'postgatk'

    // conda "/data/eande106/software/conda_envs/popgen-nf_env"

    publishDir "${params.output}/divergent_regions/Mask_DF", mode: 'copy'

    input:
        tuple val(isotype_name), path(bam), path(bai), path(vcf), path(index)

    output:
        file("*_Mask_DF.tsv")

    """
    # although calling the isotype strain, it is actually using data from ref strain.
    bcftools view -s ${isotype_name} ${vcf} | bcftools filter -i 'GT="alt"' -Ov | \\
    bedtools coverage -a ${params.bin_bed} -b stdin -counts > ${isotype_name}_variant_counts.txt

    # need to use the isotype ref strain here, not the isotype strain!
    mosdepth -b ${params.bin_bed} ${isotype_name} ${bam}

    gunzip ${isotype_name}.regions.bed.gz

    echo 'CHROM\tSTART_BIN\tEND_BIN\tVAR_COUNT\tCOVERAGE' > ${isotype_name}_Mask_DF.tsv

    paste ${isotype_name}_variant_counts.txt ${isotype_name}.regions.bed | cut -f 1-4,8 >> ${isotype_name}_Mask_DF.tsv

    """
}


process define_divergent_region {

    label 'R'

    // conda "/data/eande106/software/conda_envs/popgen-nf-r_env"

    publishDir "${params.output}/divergent_regions", mode: 'copy'
    publishDir "${params.output}/NemaScan", pattern: 'divergent_bins.bed', mode: 'copy'
    publishDir "${params.output}/NemaScan", pattern: 'divergent_df_isotype.bed', mode: 'copy'

    memory { 128.GB + 20.GB * task.attempt }
    // errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }
    errorStrategy 'ignore'

    input:
        file("*")

    output:
        tuple file("divergent_regions_strain.bed"), file("divergent_regions_all.bed"), file("divergent_regions.png")

    """
    cp ${workflow.projectDir}/bin/reoptimzied_divergent_region_characterization.Rmd reoptimzied_divergent_region_characterization.Rmd
    
    cp ${workflow.projectDir}/bin/${params.species}_chr_lengths.tsv ./df_chr_length.tsv
    Rscript -e "rmarkdown::render('reoptimzied_divergent_region_characterization.Rmd')"

    # files for NemaScan
    cp ${params.bin_bed} ./divergent_bins.bed
    cp divergent_regions_strain.bed ./divergent_df_isotype.bed

    # gzip divergent files
    cp All_divergent_regions.tsv divergent_regions_all.bed

    """
}


process get_species_sheet {

    label 'R'
    
    publishDir "${params.output}/NemaScan/", mode: 'copy'

    // conda "/data/eande106/software/conda_envs/popgen-nf-r_env"

    output:
        file("strain_isotype_lookup.tsv")

    """
    Rscript --vanilla ${workflow.projectDir}/bin/download_google_sheet.R ${params.species}
        
    """

}
