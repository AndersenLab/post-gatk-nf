#!/usr/bin/env nextflow 

nextflow.preview.dsl=2
date = new Date().format( 'yyyyMMdd' )

if (params.markers == null) error  "Parameter --markers is required."

if (params.vcf == null) error  "Parameter --pca_vcf is required."

include {convert_tree; quick_tree} from "./modules/postgatk.nf"

process create_vcf {
    
    label 'postgatk'

    input: 
        tuple file(vcf), file(markers)

    output: 
        tuple file("eiganstrat_input.vcf.gz"), file ("eiganstrat_input.vcf.gz.tbi")
    """
    bcftools view -T ${markers} -Oz -o eiganstrat_input.vcf.gz ${vcf}

    tabix eiganstrat_input.vcf.gz

    """

}


workflow{

// Index the VCF 

marker_file = Channel.fromPath(params.markers)

Channel.fromPath(params.vcf)
    .combine(marker_file) | create_vcf | convert_tree | quick_tree

} 