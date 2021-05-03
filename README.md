# Post-gatk-nf

Pipeline for vcf annotation and popgen analysis.

## Typical use for debugging:

```
nextflow main.nf --debug
```

### Typical use for new vcf:

```
nextflow main.nf --vcf path_to_vcf.vcf.gz --sample_sheet path_to_sample_sheet.tsv
```

### Parameters
    parameters              description                                            Set/Default
    ==========              ===========                                            ========================
    --debug                 Use --debug to indicate debug mode                     (optional)
    --vcf                   Hard filtered vcf to calculate variant density         (required)
    --sample_sheet          TSV with column iso-ref strain, bam, bai (no header)   (required)
    --species               Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'   c_elegans
    --project               Project name for species reference                     PRJNA13758
    --ws_build              WormBase version for species reference                 WS276
    --output                Output folder name.                                    popgen-date (in current folder)


### Overview

![Overview of post-gatk-nf](http://github.com/andersenlab/post-gatk-nf/img/post-gatk-nf-flow.png)


