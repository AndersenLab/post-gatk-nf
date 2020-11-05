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
params.species = ""

if (params.debug.toString() == "true") {

    params.vcf = "${workflow.projectDir}/test_data/WI.20200811.hard-filter.vcf.gz"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bin_bed = "${workflow.projectDir}/bin/bins_1kb_ce.bed"
    params.bam_folder = "${workflow.projectDir}/test_data/bam"
    params.output = "popgen-${date}-debug"

} else {
    // Read input
    params.vcf = ""
    params.sample_sheet = ""

    // 1kb bins for all chromosomes
    params.bin_bed = "${workflow.projectDir}/bin/bins_1kb_${params.species}.bed"

    // folder for the bam files. currently need to put all bam in the same folder
    params.bam_folder = ""

    params.output = "popgen-${date}"
}


if ( !params.species ) error "Parameter --species is required, options are: ce, ct or cb."
if ( !params.vcf ) error "Parameter --vcf is required. Specify path to the full vcf."
if ( !params.sample_sheet ) error "Parameter --sample_sheet is required. It should contain a column of strain names, a column of bam file names and a column of bai file names WITH NO HEADERS. If the bam and bai column do not contain full path to the files, specify that path with --bam_folder."


// Read input
vcf = Channel.fromPath("${params.vcf}")
vcf_index = Channel.fromPath("${params.vcf}.tbi")

// Read sample sheet
sample_sheet = Channel.fromPath(params.sample_sheet)
                       .splitCsv(sep: "\t")
                      .map { row -> [ row[0], row[1] ]}


def log_summary() {
/*
    Generates a log
*/

out = '''

Subset isotype reference strains from hard-filter vcf, build trees, define haplotypes and divergent regions.
Should work for c.e, c.b, c.t since they all have the same number of chromosomes.

''' + """

nextflow main.nf -profile quest --debug=true

nextflow main.nf -profile quest --vcf=hard-filtered.vcf --sample_sheet=sample_sheet.tsv --bam_folder=/path/bam_folder --species=ce 

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --debug              Set to 'true' to test                                    ${params.debug}
    --species            Species: 'ce', 'ct' or 'cb'                              ${params.species}
    --vcf                hard filtered vcf to calculate variant density           ${params.vcf}
    --sample_sheet       Columns are strain, bam, bai without headers             ${params.sample_sheet}
    --bam_folder         (Optional) path to prefix the bam file. No end slash.    ${params.bam_folder}
    --output             (Optional) output folder name                            ${params.output}
 
    username                                                                      ${"whoami".execute().in.text}

    HELP: http://andersenlab.org/dry-guide/
"""
out
}


log.info(log_summary())


if (params.help) {
    exit 1
}



workflow { 

    vcf.combine(vcf_index) | subset_iso_ref_strains

    vcf.combine(vcf_index).concat(subset_iso_ref_strains.out.ref_strain_vcf) | build_tree

    subset_iso_ref_strains.out.ref_strain_vcf.combine(contigs) | haplotype_sweep_IBD

    haplotype_sweep_IBD.out.collect() | haplotype_sweep_plot

    subset_iso_ref_strains.out.ref_strain_vcf.combine(sample_sheet) | count_variant_coverage

    count_variant_coverage.out.collect() | define_divergent_region

}


/* 
    ==================
    Subset isotype reference strains
    ==================
*/

process subset_iso_ref_strains {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    memory 16.GB

    publishDir "${params.output}", mode: 'copy'

    input: 
        tuple file(vcf), file(vcf_index)

    output: 
        tuple file("*hard-filter.ref_strain.vcf.gz"), file("*hard-filter.ref_strain.vcf.gz.tbi"), emit: ref_strain_vcf
        file("*.stats.txt")


    """
    cp /projects/b1059/projects/Dan/post-gatk-nf/20201101_release20201015/popgen-20201101/WI.20201015.hard-filter.ref_strain.vcf.gz WI.20201015.hard-filter.ref_strain.vcf.gz

    cp /projects/b1059/projects/Dan/post-gatk-nf/20201101_release20201015/popgen-20201101/WI.20201015.hard-filter.ref_strain.vcf.gz.tbi WI.20201015.hard-filter.ref_strain.vcf.gz.tbi

    cp /projects/b1059/projects/Dan/post-gatk-nf/20201101_release20201015/popgen-20201101/WI.20201015.hard-filter.ref_strain.vcf.gz.stats.txt WI.20201015.hard-filter.ref_strain.vcf.gz.stats.txt


#    cut -f1 ${params.sample_sheet} > strain_list.txt

#    output_vcf=`echo ${vcf} | sed 's/hard-filter.vcf.gz/hard-filter.ref_strain.vcf.gz/'`

#    bcftools view -S strain_list.txt -O u ${vcf} | bcftools view -O v --min-af 0.000001 --max-af 0.999999 | vcffixup - | bcftools view --threads 3 -O z > \${output_vcf}

#    bcftools index \${output_vcf}

#    bcftools index --tbi \${output_vcf}

#    bcftools stats -s- \${output_vcf} > \${output_vcf}.stats.txt
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
        tuple path(vcf), path(index), val(strain), val(bam)

    output:
        file("*_Mask_DF.tsv")

    """

    bcftools view -s ${strain} ${vcf} | bcftools filter -i 'GT="alt"' -Ov | \\

    bedtools coverage -a ${params.bin_bed} -b stdin -counts > ${strain}_variant_counts.txt

    mosdepth -b ${params.bin_bed} ${strain} ${params.bam_folder}/${bam}

    gunzip ${strain}.regions.bed.gz

    echo 'CHROM\tSTART_BIN\tEND_BIN\tVAR_COUNT\tCOVERAGE' > ${strain}_Mask_DF.tsv

    paste ${strain}_variant_counts.txt ${strain}.regions.bed | cut -f 1-4,8 >> ${strain}_Mask_DF.tsv

    """
}




process define_divergent_region {

    conda "/projects/b1059/software/conda_envs/popgen-nf-r_env"

    publishDir "${params.output}/divergent_regions", mode: 'copy'

    memory { 30.GB + 10.GB * task.attempt }
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
