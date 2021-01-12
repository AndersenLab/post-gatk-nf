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
params.cendr = false

if (params.debug) {

    params.hard_filtered_vcf = "${workflow.projectDir}/test_data/WI.20201230.hard-filter.vcf.gz"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bin_bed = "${workflow.projectDir}/bin/bins_1kb_ce.bed"
    params.bam_folder = "${workflow.projectDir}/test_data/bam"
    params.output = "popgen-${date}-debug"

} else {
    // Read input
    params.hard_filtered_vcf = ""
    params.sample_sheet = ""

    // 1kb bins for all chromosomes
    params.bin_bed = "${workflow.projectDir}/bin/bins_1kb_${params.species}.bed"

    // folder for the bam files. currently need to put all bam in the same folder
    params.bam_folder = ""

    params.output = "popgen-${date}"
}


// Variant annotation files. The same for debug or normal run. 

// gff. below is the one DL made.
// need to be replaced by Ryan reformatted csq.
// "${params.reference_dir}/csq/${species}.${project}.${ws_build}.csq.gff3.gz" from DEC genomes-nf gives some error for transposons
params.snpeff_vcfanno_config = "${workflow.projectDir}/bin/vcfanno_snpeff.toml"
params.bcsq_vcfanno_config = "${workflow.projectDir}/bin/vcfanno.toml"


if ( params.species==null ) error "Parameter --species is required, options are: ce, ct or cb."
if ( params.hard_filtered_vcf==null ) error "Parameter --vcf is required. Specify path to the full vcf."
if ( params.sample_sheet==null ) error "Parameter --sample_sheet is required. It should contain a column of strain names, a column of bam file names and a column of bai file names WITH NO HEADERS. If the bam and bai column do not contain full path to the files, specify that path with --bam_folder."


// Read input
input_vcf = Channel.fromPath("${params.hard_filtered_vcf}")
input_vcf_index = Channel.fromPath("${params.hard_filtered_vcf}.tbi")

// To convert ref strain to isotype names. Note only strains that match the first col will get changed, so won't impact ct and cb.
isotype_convert_table = Channel.fromPath("${workflow.projectDir}/bin/ref_strain_isotype.tsv") 

// Read sample sheet
// No header. Columns are strain name, bam and bai
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
    --vcf                hard filtered vcf to calculate variant density           ${params.hard_filtered_vcf}
    --sample_sheet       TSV with column iso-ref strain, bam, bai. no header      ${params.sample_sheet}
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

    input_vcf.combine(input_vcf_index).combine(isotype_convert_table) | subset_iso_ref_strains

    // snpeff annotation
    subset_iso_ref_strains.out
      .combine(Channel.fromPath(params.snpeff_vcfanno_config))
      .combine(Channel.fromPath(params.dust_bed))
      .combine(Channel.fromPath(params.dust_bed + ".tbi"))
      .combine(Channel.fromPath(params.repeat_masker_bed))
      .combine(Channel.fromPath(params.repeat_masker_bed + ".tbi")) | snpeff_annotate_vcf

    // bcsq annotation
    subset_iso_ref_strains.out
      .combine(Channel.fromPath(params.csq_gff)) | bcsq_annotate_vcf

    bcsq_annotate_vcf.out.bcsq_tsv
      .combine(Channel.fromPath(params.AA_length))
      .combine(Channel.fromPath(params.AA_score)) | prep_other_annotation

    bcsq_annotate_vcf.out.bcsq_vcf
      .combine(prep_other_annotation.out)
      .combine(Channel.fromPath(params.bcsq_vcfanno_config))
      .combine(Channel.fromPath(params.dust_bed))
      .combine(Channel.fromPath(params.dust_bed + ".tbi"))
      .combine(Channel.fromPath(params.repeat_masker_bed))
      .combine(Channel.fromPath(params.repeat_masker_bed + ".tbi")) | AA_annotate_vcf

/*
    if (params.cendr == true) {

      /*
      // Generate Strain-level TSV and VCFs. In 20200815 release this step was removed. but I'm keeping the code here in case having single-strain vcf or tsv helps with displaying annotation on cendr
      // note 1: since annotation is only done for isotype-ref strains, here also only include isotype ref strains.
      // note 2: this step used vcf with only bcsq annotation b/c I'm not sure how to split into single sample vcf with protein length and amino acid score. Ryan has the code to split out single strain annotation for bcsq.
      
      bcsq_annotate_vcf.out | strain_list
      strain_set = strain_list.out.splitText( it.strip() )

      strain_set.combine( bcsq_annotate_vcf.out.anno_vcf ) | generate_strain_tsv
      

      // Extract SNPeff severity tracks
      mod_tracks = Channel.from(["LOW", "MODERATE", "HIGH", "MODIFIER"])
      bcsq_annotate_vcf.out.anno_vcf.spread(mod_tracks) | generate_severity_tracks

      }
*/


    input_vcf.combine(input_vcf_index).concat(subset_iso_ref_strains.out) | build_tree

    subset_iso_ref_strains.out.combine(contigs) | haplotype_sweep_IBD

    haplotype_sweep_IBD.out.collect() | haplotype_sweep_plot

    subset_iso_ref_strains.out.combine(sample_sheet) | count_variant_coverage

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
    cpus 4

    publishDir "${params.output}", mode: 'copy'

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

    # convert to isotype names
    bcftools reheader -s ref_strain_isotype.tsv WI.ref_strain.vcf.gz | bcftools view -O z > \${output}
    bcftools index --tbi \${output}
    """

}


/* 
    =========================
    Annotate isotype-only vcf
    =========================
*/


process snpeff_annotate_vcf {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    publishDir "${params.output}/snpeff", mode: 'copy'

    input:
        tuple file(vcf), file(vcf_index), \
              path(vcfanno), \
              path("dust.bed.gz"), \
              path("dust.bed.gz.tbi"), \
              path("repeat_masker.bed.gz"), \
              path("repeat_masker.bed.gz.tbi")

    output:
        tuple path("*.snpeff.vcf.gz"), path("*.snpeff.vcf.gz.tbi"), path("snpeff.stats.csv")


    script:
    """
        output=`echo ${vcf} | sed 's/.vcf.gz/.snpeff.vcf.gz/'`

        bcftools view -O v ${vcf} | \\
        snpEff eff -csvStats snpeff.stats.csv \\
                   -no-downstream \\
                   -no-intergenic \\
                   -no-upstream \\
                   -nodownload \\
        -dataDir ${params.snpeff_dir} \\
        -config ${params.snpeff_dir}/snpEff.config \\
        ${params.snpeff_reference} | \\
        bcftools view -O z > out.vcf.gz

        vcfanno ${vcfanno} out.vcf.gz | bcftools view -O z > \${output}
        bcftools index --tbi \${output}
    """

}



process bcsq_annotate_vcf {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    input:
        tuple file(vcf), file(vcf_index), file(gff)


    output:
        tuple file("bcsq.vcf.gz"), file("bcsq.vcf.gz.tbi"), emit:bcsq_vcf
        path "BCSQ.tsv", emit: bcsq_tsv

    script:
    """
        gzip -dc $gff | grep -v 'transposon' | bgzip > csq.gff.gz
        tabix -p gff csq.gff.gz

        bcftools csq -O z --fasta-ref ${params.reference} \\
                     --gff-annot csq.gff.gz \\
                     --phase a $vcf > bcsq.vcf.gz

        bcftools index --tbi bcsq.vcf.gz

        bcftools query -e 'INFO/BCSQ="."' -f '%CHROM %POS %BCSQ\\n' bcsq.vcf.gz > BCSQ.tsv
    """

}


process prep_other_annotation {

    conda "/projects/b1059/software/conda_envs/popgen-nf-r_env"
    echo true

    input:
        tuple path("BCSQ.tsv"), path("gff_AA_Length.tsv"), path("AA_Scores.tsv")

    output:
        path("BCSQ_bed.bed")

    """
      Rscript --vanilla ${workflow.projectDir}/bin/Parse_and_Score.R BCSQ.tsv AA_Scores.tsv gff_AA_Length.tsv
    """

}


process AA_annotate_vcf {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    publishDir "${params.output}/bcsq", mode: 'move'

    input:
        tuple file(vcf), file(vcf_index), file("BCSQ_bed.bed"), file(vcfanno), file("dust.bed.gz"), file("dust.bed.gz.tbi"), file("repeat_masker.bed.gz"), file("repeat_masker.bed.gz.tbi")

    output:
        tuple file("*hard-filter.ref_strain.vcf.gz"), file("*hard-filter.ref_strain.vcf.gz.tbi")
        file("*.stats.txt")


    """
        output_vcf=`basename ${params.hard_filtered_vcf} | sed 's/hard-filter.vcf.gz/hard-filter.ref_strain.vcf.gz/'`

        bgzip BCSQ_bed.bed
        tabix BCSQ_bed.bed.gz

        vcfanno ${vcfanno} $vcf | bcftools view -O z > \${output_vcf}

        bcftools index --tbi \${output_vcf}

        bcftools stats -s- \${output_vcf} > \${output_vcf}.stats.txt

    """
}




/* 
    ==============================
    Generate individual strain vcf. This is copied from wi-gatk in case it's needed for cendr browser.
    ==============================
*/

/*
process strain_list {

    input:
        tuple path(vcf), path(vcf_index)
    
    output:
        file("samples.txt")

    """
        bcftools query --list-samples ${vcf} > samples.txt
    """
}


process generate_strain_tsv {
    // Generate a single TSV for every strain.

    tag { strain }

    publishDir "${params.output}/strain/tsv", mode: 'copy', pattern: "*.tsv.gz*"

    input:
        tuple val(strain), path(vcf), file(vcf_index)

    output:
        tuple path("${strain}.${date}.tsv.gz"),  path("${strain}.${date}.tsv.gz.tbi")

    """
        # Generate TSV
        {
            echo -e 'CHROM\\tPOS\\tREF\\tALT\\tFILTER\\tFT\\tGT';
            bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%FT\t%TGT]\n' --samples ${strain} ${vcf};
        } | awk -F'\\t' -vOFS='\\t' '{ gsub("\\\\.", "PASS", \$6) ; print }' > ${strain}.${date}.tsv

        bgzip ${strain}.${date}.tsv
        tabix -S 1 -s 1 -b 2 -e 2 ${strain}.${date}.tsv.gz
    """

}
*/



/* 
    ==========================================
    Generate severity tracks for cendr browser
    ==========================================
*/


process generate_severity_tracks {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    publishDir "${params.output}/tracks", mode: 'copy'

    tag { severity }

    input:
        tuple path("in.vcf.gz"), path("in.vcf.gz.tbi"), val(severity)
    output:
        set file("${date}.${severity}.bed.gz"), file("${date}.${severity}.bed.gz.tbi")

    """
        bcftools view in.vcf.gz | \
        grep -x -f ${workflow.projectDir}/bin/${severity}.txt | \
        awk '\$0 !~ "^#" { print \$1 "\\t" (\$2 - 1) "\\t" (\$2)  "\\t" \$1 ":" \$2 "\\t0\\t+"  "\\t" \$2 - 1 "\\t" \$2 "\\t0\\t1\\t1\\t0" }' | \\
        bgzip  > ${date}.${severity}.bed.gz
        tabix -p bed ${date}.${severity}.bed.gz
        fsize=\$(zcat ${date}.${severity}.bed.gz | wc -c)
        if [ \${fsize} -lt 2000 ]; then
            exit 1
        fi;
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
