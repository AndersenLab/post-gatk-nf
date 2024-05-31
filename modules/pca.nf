/*
==================================
~ > *                        * < ~
~ ~ > *                    * < ~ ~
~ ~ ~ > *  ANNOTATE VCF  * < ~ ~ ~
~ ~ > *                    * < ~ ~
~ > *                        * < ~
==================================
*/


/*
------------ Extract ancestor strain from the VCF and make bed file for annotations 
*/

process extract_ancestor_bed {

    label 'postgatk'
    label 'xs'

    publishDir "${params.output}/ANNOTATE_VCF", mode: 'copy'

    input:
      tuple file(vcf), file(vcfindex)

    output:
      tuple file("ANC.bed.gz"), file("ANC.bed.gz.tbi")

      """
        bcftools query --samples ${params.anc} -f '%CHROM\\t%POS\\t%END\\t[%TGT]\\n' ${vcf} |\\
        awk -F"/" '\$1=\$1' OFS="\\t" |\\
        awk '{print \$1, \$2 = \$2 - 1, \$3, \$4}' OFS="\\t" |\\
        bgzip > ANC.bed.gz

        tabix ANC.bed.gz
        echo "ANCESTOR DONE"
      """
}

/*
------------ Annotate small variant VCF  
*/

process annotate_small_vcf {

    publishDir "${params.output}/ANNOTATE_VCF", mode: 'copy'

    label 'postgatk'
    label 'xs'

    // conda '/projects/b1059/software/conda_envs/vcffixup'

    input:
      tuple file(vcf), file(vcfindex), file("ANC.bed.gz"), file("ANC.bed.gz.tbi"), file(pop) //, val(pop), val(maf), val(sm)

    output:
      tuple file("Ce330_annotated.vcf.gz"), file("Ce330_annotated.vcf.gz.tbi")


      """
        # get vcfanno files
        cp ${workflow.projectDir}/input_files/annotations/${params.species}/* .
        cat ${params.vcfanno_config} | sed 's/species/${params.species}/' > anno_config.toml

        bcftools view -S ${pop} ${vcf} -Oz -o population-filtered.vcf.gz 

        vcfanno anno_config.toml population-filtered.vcf.gz |\\
        awk '\$0 ~ "#" || \$0 !~ "Masked" {print}' |\\
        vcffixup - |\\
        bcftools filter -i N_MISSING=0 -Oz -o Ce330_annotated.vcf.gz

        tabix -p vcf Ce330_annotated.vcf.gz
      """
}



/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  Run PCA and DAPC  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Prepare files for EIGENSTRAT
*/

process get_singletons {
  tag {"PREPARE EIGENSTRAT FILES"}

  label 'vcftools'
  label 'xs'
  
  input:
      tuple file(vcf), file(vcfindex)
  output: 
      file("singleton_ids.txt")
      
  """
  vcftools --gzvcf ${vcf} --singletons
  cat out.singletons | awk 'BEGIN{OFS=":"} {if (NR!=1) print \$1,\$2}' > singleton_ids.txt
  """
}

process vcf_to_eigstrat_files {

  tag {"PREPARE EIGENSTRAT FILES"}

  label 'postgatk'
  label 'xs'

  // conda '/projects/b1059/software/conda_envs/vcffixup'

  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/INPUTFILES", mode: 'copy'

  input:
    tuple file(vcf), file(vcfindex), val("test_ld"), file(singleton_ids)

  output:
    tuple val(test_ld), file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("eiganstrat_markers.txt"),  file ("sorted_samples.txt"), emit : eig_strat_inputs
    tuple file("ce_norm.vcf.gz"), file ("eiganstrat_input.vcf.gz"), emit: vcf_file

    """

    bcftools view --regions I,II,III,IV,V,X ${vcf} |\\
    bcftools norm -m + -Oz -o ce_norm.vcf.gz

    tabix -p vcf ce_norm.vcf.gz

    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 ${test_ld} --allow-extra-chr 

    plink --vcf ce_norm.vcf.gz --biallelic-only --set-missing-var-ids @:# --extract plink.prune.in --exclude ${singleton_ids} --geno 0 --recode 12 --out eigenstrat_input --allow-extra-chr --make-bed   
   
    awk '{print \$2}' OFS="\t" eigenstrat_input.bim | 
    awk -F":" '\$1=\$1' OFS="\t"| sort -k1,1d -k2,2n > eiganstrat_markers.txt

    bcftools view -T eiganstrat_markers.txt -Oz -o eiganstrat_input.vcf.gz ce_norm.vcf.gz

    bcftools query -l ce_norm.vcf.gz |\\
    sort > sorted_samples.txt

    cat eigenstrat_input.bim |\\
    sed 's/^III/3/g' |\\
    sed 's/^II/2/g' |\\
    sed 's/^IV/4/g' |\\
    sed 's/^I/1/g' |\\
    sed 's/^V/5/g' |\\
    cut -f-6 > eigenstrat_input.pedsnp      

    cut -f-6 -d' ' eigenstrat_input.ped |\\
    awk '{print 1, \$2, \$3, \$3, \$5, 1}'  > eigenstrat_input.pedind

    echo "rerun"
    """

}


/*
------------ input eigan data to tree (normal tree has different filtering steps)
*/

process pca_tree {

  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/Tree/", mode: 'copy'

  label 'pca'
  label 'xs'

  // conda '/projects/b1059/software/conda_envs/vcffixup'

  input:
    tuple val("test_ld"), file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file (".txt"), file ("sorted_samples.txt"), file(eigenparameters)

  output:
    tuple val(test_ld), file("eigenstrat_no_removal.evac"), file("eigenstrat_no_removal.eval"), file("logfile_no_removal.txt"), \
    file("eigenstrat_no_removal_relatedness"), file("eigenstrat_no_removal_relatedness.id"), file("TracyWidom_statistics_no_removal.tsv")


    """

    smartpca -p ${eigenparameters} > logfile_no_removal.txt

    sed -n -e '/Tracy/,\$p' logfile_no_removal.txt |\
    sed -e '/kurt/,\$d' |\
    awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
    sed -e "s/[[:space:]]\\+/ /g" |\
    sed 's/^ //g' |\
    awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
    awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_no_removal.tsv
    """

}



/*
------------ Run EIGENSTRAT without removing outlier strains
*/

process run_eigenstrat_no_outlier_removal {

  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/NO_REMOVAL/", mode: 'copy'

  label 'pca'
  label 'xs'

  // conda '/projects/b1059/software/conda_envs/vcffixup'

  input:
    tuple val("test_ld"), file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("eiganstrat_markers.txt"), file ("sorted_samples.txt"), file(eigenparameters)

  output:
    tuple val(test_ld), file("eigenstrat_no_removal.evac"), file("eigenstrat_no_removal.eval"), file("logfile_no_removal.txt"), \
    file("eigenstrat_no_removal_relatedness"), file("eigenstrat_no_removal_relatedness.id"), file("TracyWidom_statistics_no_removal.tsv")


    """

    smartpca -p ${eigenparameters} > logfile_no_removal.txt

    sed -n -e '/Tracy/,\$p' logfile_no_removal.txt |\
    sed -e '/kurt/,\$d' |\
    awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
    sed -e "s/[[:space:]]\\+/ /g" |\
    sed 's/^ //g' |\
    awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
    awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_no_removal.tsv
    """

}

/*
------------ Run EIGENSTRAT with removing outlier strains
*/

process run_eigenstrat_with_outlier_removal {

  // conda '/projects/b1059/software/conda_envs/vcffixup'

  label 'pca'
  label 'xs'

  publishDir "${params.output}/EIGESTRAT/LD_${test_ld}/OUTLIER_REMOVAL/${num_outliers}/", mode: 'copy'

  input:
    tuple val("test_ld"), file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("eiganstrat_markers.txt"), file ("sorted_samples.txt"), file(eigenparameters), val("num_outliers")

  output:
    tuple val(test_ld), file("eigenstrat_outliers_removed.evac"), file("eigenstrat_outliers_removed.eval"), file("logfile_outlier.txt"), \
    file("eigenstrat_outliers_removed_relatedness"), file("eigenstrat_outliers_removed_relatedness.id"), file("TracyWidom_statistics_outlier_removal.tsv"), val("num_outliers")

   
    """
    bash ${workflow.projectDir}/bin/edit_outlier_param.sh -f ${eigenparameters} -n ${num_outliers}

    smartpca -p outlier_eigpar > logfile_outlier.txt

    sed -n -e '/Tracy/,\$p' logfile_outlier.txt |\
    sed -e '/kurt/,\$d' |\
    awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
    sed -e "s/[[:space:]]\\+/ /g" |\
    sed 's/^ //g' |\
    awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
    awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_outlier_removal.tsv
    """

}


/*
------------ Run HTML report for PCA analysis
*/

process HTML_report_PCA {

  label 'R'
  label 'xs'

  // conda '/projects/b1059/software/conda_envs/cegwas2-nf_env'

  publishDir "${params.output}/", mode: 'copy'


  input:
  tuple val("test_ld"), file("eigenstrat_no_removal.evac"), file("eigenstrat_no_removal.eval"), file("logfile_no_removal.txt"), \
    file("eigenstrat_no_removal_relatedness"), file("eigenstrat_no_removal_relatedness.id"), file("TracyWidom_statistics_no_removal.tsv"), \
    file("eigenstrat_outliers_removed.evac"), file("eigenstrat_outliers_removed.eval"), file("logfile_outlier.txt"), \
    file("eigenstrat_outliers_removed_relatedness"), file("eigenstrat_outliers_removed_relatedness.id"), \
    file("TracyWidom_statistics_outlier_removal.tsv"), file(pca_report), file(pca_template)


  output:
   tuple file("pca_report.Rmd"), file('pca_template.Rmd'), file("*.html")


  """
  # prepare for report
  cat ${pca_report} | \\
  sed "s+LD_VALUE+${test_ld}+" | \\
  sed "s+EIGESTRAT/{test_ld}/NO_REMOVAL/++g" | \\
  sed "s+EIGESTRAT/{test_ld}/OUTLIER_REMOVAL/++g" > pca_report_LD_${test_ld}.Rmd

  Rscript -e "rmarkdown::render('pca_report_LD_${test_ld}.Rmd')"
  """


}

//  cat pca_report.Rmd | sed 's+pca_template.Rmd+${workflow.projectDir}/bin/pca_template.Rmd+' > new_pca_report.Rmd



