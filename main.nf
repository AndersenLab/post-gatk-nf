#!/usr/bin/env nextflow 
/*
    Authors:
    - Dan Lu <dan.lu@northwestern.edu>
    - Katie Evans <katiesevans9@gmail.com>
*/

nextflow.enable.dsl=2
// NXF_VER=23.0" Require later version of nextflow
//assert System.getenv("NXF_VER") >= "23.0"

date = new Date().format( 'yyyyMMdd' )
contigs = Channel.from("I","II","III","IV","V","X")

if (params.debug) {
    params.vcf_folder = "${workflow.projectDir}/test_data/"
    params.vcf = "${params.vcf_folder}/WI.20201230.hard-filter.vcf.gz"
    params.soft_vcf = "${params.vcf}"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bam_folder = "${workflow.projectDir}/test_data/bam"
    params.output = "popgen-${date}-debug"
    params.species = "c_elegans"
    //params.anc = "ECA718"
    params.pops = "${workflow.projectDir}/input_files/WI_328_isotype.tsv"
    params.eigen_ld = "0.8,0.6"
} else {
    // Read input
    // params.vcf_folder = ""
    // params.sample_sheet = ""
    params.vcf = "${params.vcf_folder}/*.hard-filter.vcf.gz" 
    params.soft_vcf = "${params.vcf_folder}/*.soft-filter.vcf.gz"

    // folder for the bam files. currently need to put all bam in the same folder
    params.bam_folder = "${dataDir}/${params.species}/WI/alignments/"
    params.output = "popgen-${date}"
    // PCA params
    params.anc = null
    params.pops = null
    params.eigen_ld = null
    params.pca_vcf = null
    params.outlier_iterations = "5"
}

// more params
params.vcfanno_config = "${workflow.projectDir}/input_files/ANNOTATION_conf.toml"
params.eigen_par_outlier_removal = "${workflow.projectDir}/bin/eigpar"
params.eigen_par_no_removal = "${workflow.projectDir}/bin/eigpar_no_removal"
params.R_libpath = "${params.softwareDir}/R_lib_3.6.0"
params.snps = '--snps-only'


// Note that params.species is set in the config to be c_elegans (default)
if ( params.species==null ) error "Parameter --species is required. Please select c_elegans, c_briggsae, or c_tropicalis."

if(params.postgatk) {
    if ( params.vcf==null ) error "Parameter --vcf is required. Specify path to the folder containing a 'hard-filter.vcf.gz' and a 'soft-filter.vcf.gz'."
    if ( params.sample_sheet==null ) error "Parameter --sample_sheet is required. It should contain a column of strain names, a column of bam file names and a column of bai file names WITH NO HEADERS. If the bam and bai column do not contain full path to the files, specify that path with --bam_folder."
}

// check pca inputs
if(params.pca && !params.postgatk) {
    // if(params.pops == null) error "Parameter --pops is required. Specify path to file"
    // if(params.anc == null) error "Parameter --anc is required. Specify ancestor strain"
    if(params.eigen_ld == null) error "Parameter --eigen_ld is required. Specify LD value(s)"
    if(params.pca_vcf == null) error "Parameter --pca_vcf is required. Specify path to SNV-filtered VCF"
}

if(params.pca && params.postgatk) {
    //if(params.anc == null) error "Parameter --anc is required. Specify ancestor strain"
    if(params.eigen_ld == null) error "Parameter --eigen_ld is required. Specify LD value(s)"
}

// 1kb bins for all chromosomes
params.bin_bed = "${workflow.projectDir}/bin/bins_1kb_${params.species}.bed"



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
    --vcf_folder         Path to folder containing hard and soft VCF              ${params.vcf_folder}
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

// import the pca module
include {extract_ancestor_bed; annotate_small_vcf; vcf_to_eigstrat_files; run_eigenstrat_no_outlier_removal; run_eigenstrat_with_outlier_removal; HTML_report_PCA; get_singletons} from './modules/pca.nf'
include {subset_iso_ref_strains; subset_iso_ref_soft; subset_snv; make_small_vcf; convert_tree; quick_tree; plot_tree; haplotype_sweep_IBD; haplotype_sweep_plot; 
    define_divergent_region; prep_variant_coverage; count_variant_coverage; get_species_sheet} from './modules/postgatk.nf'


workflow { 

    // if only pca profile, don't rerun this stuff
    if(params.postgatk) {
        // make species sheet for nemascan
        get_species_sheet()

        // Read input
        input_vcf = Channel.fromPath("${params.vcf}")
        input_vcf_index = Channel.fromPath("${params.vcf}.tbi")

        // Also use soft-filter vcf from same path location as hard filtered
        soft_vcf = Channel.fromPath("${params.soft_vcf}")
        soft_vcf_index = Channel.fromPath("${params.soft_vcf}.tbi")

        // To convert ref strain to isotype names. Note only strains that match the first col will get changed, so won't impact ct and cb.
        isotype_convert_table = Channel.fromPath("${workflow.projectDir}/bin/ref_strain_isotype.tsv") 

        // Read sample sheet
        // No header. Columns are strain name, bam and bai
        sample_sheet = Channel.fromPath(params.sample_sheet)
                               .splitCsv(sep: "\t")
                              .map { row -> [ row[0], row[1] ]}

        // subset isotype ref strains
        input_vcf.combine(input_vcf_index).combine(isotype_convert_table) | subset_iso_ref_strains
        subset_iso_ref_strains.out.vcf | make_small_vcf
        soft_vcf.combine(soft_vcf_index).combine(isotype_convert_table) | subset_iso_ref_soft

        // subset SNV vcf for imputation
        subset_iso_ref_strains.out.vcf | subset_snv

        // // build tree
        // input_vcf.combine(input_vcf_index).concat(subset_iso_ref_strains.out.vcf) | convert_tree | quick_tree | plot_tree

        // haplotype
        subset_iso_ref_strains.out.vcf.combine(contigs) | haplotype_sweep_IBD
        subset_iso_ref_strains.out.vcf.combine(sample_sheet).combine(isotype_convert_table) | prep_variant_coverage 

        prep_variant_coverage.out
            .collectFile(name: "test.tsv", keepHeader: true)
            .splitCsv(header:true)
            .map { row -> [ row.strain, file(row.bam), file(row.bai) ] }
            .combine(subset_iso_ref_strains.out.vcf) | count_variant_coverage


        // // haplotype_sweep_plot and define_divergent_region always give error during debugging run prob b/c the debug dataset is too small. so turn it off when debugging
        // if (!params.debug) {
        //   haplotype_sweep_IBD.out.collect() | haplotype_sweep_plot
        //   count_variant_coverage.out.collect() | define_divergent_region
        // }

        // // collect snv vcf for pca
        // pca_vcf = subset_snv.out

        // // generate pops file for pca
        // pop_strains = subset_iso_ref_strains.out.pop_strains
    }

    // PCA analysis
    if(params.pca) {
        // initialize population 
        // File pop_file = new File(params.pops);
        // pop_strains = Channel.from(pop_file.collect { it.tokenize( ' ' ) })
        //              .map { POP, MAF, SM -> [POP, MAF, SM] }

        // check if pca only or also postgatk
        if(!params.postgatk) {
            pca_vcf = Channel.fromPath("${params.pca_vcf}").combine(Channel.fromPath("${params.pca_vcf}.tbi"))
            pop_strains = Channel.fromPath("${params.pops}")
        }

        // extract ancestor
        // pca_vcf | extract_ancestor_bed

        // annotate small vcf
        // pca_vcf 
         // .combine(extract_ancestor_bed.out)
         // .combine(pop_strains) | annotate_small_vcf 

        ld_range = Channel.of("${params.eigen_ld}")
                      .splitCsv()
                      .flatMap { it }

        // make vcf for eigenstrat - use LD provided
        if(params.singletons){
            pca_vcf | get_singletons

            pca_vcf
            .combine(ld_range)
            .combine(get_singletons.out)| vcf_to_eigstrat_files

        }
        else{
            singleton_ids = Channel.fromPath("${workflow.projectDir}/bin/blank_snps.txt")
        
            pca_vcf
                .combine(ld_range)
                .combine(singleton_ids) | vcf_to_eigstrat_files
        }
        
        //just a path to a blank file so we don't filter anything
        vcf_to_eigstrat_files.out.eig_strat_inputs
          .combine(Channel.fromPath(params.eigen_par_no_removal)) | run_eigenstrat_no_outlier_removal

        outlier_its = Channel.of("${params.outlier_iterations}")
                        .splitCsv()
                        .flatMap{ it }
        
        vcf_to_eigstrat_files.out.eig_strat_inputs
          .combine(Channel.fromPath(params.eigen_par_outlier_removal)) 
          .combine(outlier_its) | run_eigenstrat_with_outlier_removal

        // run html report
        // not functional quite yet...

        // run_eigenstrat_no_outlier_removal.out
          //   .join(run_eigenstrat_with_outlier_removal.out)
          //   .combine(Channel.fromPath("${workflow.projectDir}/bin/pca_report.Rmd"))
          //   .combine(Channel.fromPath("${workflow.projectDir}/bin/pca_template.Rmd"))| HTML_report_PCA
    }
    
}


workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    { Parameters }
    ---------------------------
    Species: ${params.species}
    Vcf_folder: ${params.vcf_folder}
    Sample_sheet: ${params.sample_sheet}
    Output: ${params.output}
    """

    println summary

}



