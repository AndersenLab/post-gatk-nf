# Post-gatk-nf

Pipeline for simple popgen analysis.

## Typical use for debugging:

```
nextflow main.nf --debug
```

### Typical use for new vcf:

```
nextflow main.nf --vcf path_to_vcf.vcf.gz --sample_sheet path_to_sample_sheet.tsv --species <species>
```

### Parameters
    parameters              description                                            Set/Default
    ==========              ===========                                            ========================
    --debug                 Use --debug to indicate debug mode                     (optional)
    --vcf                   Hard filtered vcf to calculate variant density         (required)
    --sample_sheet          TSV with column iso-ref strain, bam, bai (no header)   (required)
    --species               Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'   c_elegans
    --output                Output folder name.                                    popgen-date (in current folder)


### Overview

![Overview of post-gatk-nf](https://github.com/AndersenLab/post-gatk-nf/blob/main/img/post-gatk-nf-flow.png?raw=true)


