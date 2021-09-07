#!/usr/bin/env nextflow 
/*
    Authors:
    - Dan Lu <dan.lu@northwestern.edu>
*/

nextflow.preview.dsl=2
// NXF_VER=20.01.0" Require later version of nextflow
//assert System.getenv("NXF_VER") == "20.01.0"

date = new Date().format( 'yyyyMMdd' )
contigs = Channel.from("I","II","III","IV","V","X")
// 1kb bins for all chromosomes
params.bin_bed = "${workflow.projectDir}/bin/bins_1kb_${params.species}.bed"

if (params.debug) {
    params.vcf = "${workflow.projectDir}/test_data/WI.20201230.hard-filter.vcf.gz"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bam_folder = "${workflow.projectDir}/test_data/bam"
    params.output = "popgen-${date}-debug"

} else {
    // Read input
    params.vcf = ""
    params.sample_sheet = ""

    // folder for the bam files. currently need to put all bam in the same folder
    params.bam_folder = "/projects/b1059/data/${params.species}/WI/alignments/"
    params.output = "popgen-${date}"
}


// Note that params.species is set in the config to be c_elegans (default)
if ( params.vcf==null ) error "Parameter --vcf is required. Specify path to the full vcf."
if ( params.sample_sheet==null ) error "Parameter --sample_sheet is required. It should contain a column of strain names, a column of bam file names and a column of bai file names WITH NO HEADERS. If the bam and bai column do not contain full path to the files, specify that path with --bam_folder."
if ( params.species==null ) error "Parameter --species is required. Please select c_elegans, c_briggsae, or c_tropicalis."



def log_summary() {
/*
    Generates a log
*/

out = '''

Build trees, define haplotypes and divergent regions.

''' + """

      * * * *                    **           * * * *    * * *    * * * *    *   *                         *
     *       *                * * * * *     *        *  *     *      *       *  *                         * *
    *        *                   **         *           *     *      *       * *                         * *
   *        *   * * * * * *      **    ***  *           * * * *      *       * *      ***      *          *
  * * * * *    *   * *   *  *    **         *    * * *  *     *      *       *  *             * * *      *
 *            *     *   *   *   *  *        *        *  *     *      *       *   *           *     *    *   *
*              * * *   * * * * *    *        * * * * *  *     *      *       *    *         *      * * * * *  
                                                                                                      **
                                                                                                     * * 
                                                                                                    *  *
                                                                                                   *  *
                                                                                                    *

nextflow main.nf -profile quest --debug

nextflow main.nf -profile quest --vcf=hard-filtered.vcf --sample_sheet=sample_sheet.tsv --bam_folder=/path/bam_folder --species=c_elegans 

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --debug              Set to 'true' to test                                    ${params.debug}
    --species            Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     ${params.species}
    --vcf                hard filtered, isotype-filtered vcf                      ${params.vcf}
    --sample_sheet       TSV with column iso-ref strain, bam, bai. no header      ${params.sample_sheet}
    --output             (Optional) output folder name                            ${params.output}
 
    username                                                                      ${"whoami".execute().in.text}

    HELP: http://andersenlab.org/dry-guide/pipeline-postGATK   
    ----------------------------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId] 
"""
out
}


log.info(log_summary())


if (params.help) {
    exit 1
}


workflow { 

    // Read input
    input_vcf = Channel.fromPath("${params.vcf}")
    input_vcf_index = Channel.fromPath("${params.vcf}.tbi")

    // To convert ref strain to isotype names. Note only strains that match the first col will get changed, so won't impact ct and cb.
    isotype_convert_table = Channel.fromPath("${workflow.projectDir}/bin/ref_strain_isotype.tsv") 

    // Read sample sheet
    // No header. Columns are strain name, bam and bai
    sample_sheet = Channel.fromPath(params.sample_sheet)
                           .splitCsv(sep: "\t")
                          .map { row -> [ row[0], row[1] ]}

    //subset isotype ref strains
    input_vcf.combine(input_vcf_index).combine(isotype_convert_table) | subset_iso_ref_strains

    // build tree
    input_vcf.combine(input_vcf_index).concat(subset_iso_ref_strains.out) | build_tree

    // haplotype
    subset_iso_ref_strains.out.combine(contigs) | haplotype_sweep_IBD
    subset_iso_ref_strains.out.combine(sample_sheet).combine(isotype_convert_table) | count_variant_coverage

    // haplotype_sweep_plot and define_divergent_region always give error during debugging run prob b/c the debug dataset is too small. so turn it off when debugging
    if (!params.debug) {
      haplotype_sweep_IBD.out.collect() | haplotype_sweep_plot
      count_variant_coverage.out.collect() | define_divergent_region
    }

}


/* 
    ==================
    Subset isotype reference strains
    ==================
*/

process subset_iso_ref_strains {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    memory 16.GB
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



/* 
    ==================
    Tree (only uses SNVs)
    ==================
*/

process build_tree {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    memory 20.GB

    publishDir "${params.output}/tree", mode: 'copy'

    input:
        tuple file(vcf), file(vcf_index)

    output:
        file("*.tree")

"""

    output_phylip=`echo ${vcf} | sed 's/.vcf.gz/.min4.phy/'`

    output_stockholm=`echo \${output_phylip} | sed 's/.phy/.stockholm/'`

    output_tree=`echo \${output_phylip} | sed 's/.phy/.tree/'`


    python ${workflow.projectDir}/bin/vcf2phylip.py -i ${vcf}

    bioconvert phylip2stockholm \${output_phylip}

    quicktree -in a -out t \${output_stockholm} > \${output_tree}

"""

}


/* 
    ==================
    Haplotype and sweep (use both SNVs and INDELs)
    ==================
*/

process haplotype_sweep_IBD {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

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

    conda "/projects/b1059/software/conda_envs/popgen-nf-r_env"

    memory { 10.GB + 10.GB * task.attempt }
    errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }

    publishDir "${params.output}/haplotype", mode: 'copy'

    input: path("*")

    output: path("*")

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



process count_variant_coverage {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    publishDir "${params.output}/divergent_regions/Mask_DF", mode: 'copy'

    input:
        tuple path(vcf), path(index), val(strain), val(bam), path("ref_strain_isotype.tsv")

    output:
        file("*_Mask_DF.tsv")

    """
    # if c.e and the strain is one of those that isotype name differs from strain name, convert isotype ref strain name to isotype name
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

    bcftools view -s \${isotype_name} ${vcf} | bcftools filter -i 'GT="alt"' -Ov | \\

    bedtools coverage -a ${params.bin_bed} -b stdin -counts > \${isotype_name}_variant_counts.txt

    mosdepth -b ${params.bin_bed} \${isotype_name} ${params.bam_folder}/${bam}

    gunzip \${isotype_name}.regions.bed.gz

    echo 'CHROM\tSTART_BIN\tEND_BIN\tVAR_COUNT\tCOVERAGE' > \${isotype_name}_Mask_DF.tsv

    paste \${isotype_name}_variant_counts.txt \${isotype_name}.regions.bed | cut -f 1-4,8 >> \${isotype_name}_Mask_DF.tsv

    """
}




process define_divergent_region {

    conda "/projects/b1059/software/conda_envs/popgen-nf-r_env"

    publishDir "${params.output}/divergent_regions", mode: 'copy'

    memory { 128.GB + 20.GB * task.attempt }
    errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }

    input:
        file("*")

    output:
        file("divergent_regions_strain.bed")

    """
    cp ${workflow.projectDir}/bin/reoptimzied_divergent_region_characterization.Rmd reoptimzied_divergent_region_characterization.Rmd
    Rscript -e "rmarkdown::render('reoptimzied_divergent_region_characterization.Rmd')"

    """
}


