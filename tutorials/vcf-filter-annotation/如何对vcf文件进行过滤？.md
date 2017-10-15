# 如何对vcf文件进行过滤？

原文链接：[ (howto) Apply hard filters to a call set](https://software.broadinstitute.org/gatk/documentation/article?id=2806)

这个内容并不多，只是看起来步骤有点多。

#### Objective

Apply hard filters to a variant callset that is too small for VQSR or for which truth/training sets are not available.

#### Caveat

This document is intended to illustrate how to compose and run the commands involved in applying the hard filtering method. The annotations and values used may not reflect the most recent recommendations. Be sure to read the documentation about [why you would use hard filters](https://www.broadinstitute.org/gatk/guide/article?id=3225) and [how to understand and improve upon the generic hard filtering recommendations](https://www.broadinstitute.org/gatk/guide/article?id=6925) that we provide.

#### Steps

1. Extract the SNPs from the call set
2. Determine parameters for filtering SNPs
3. Apply the filter to the SNP call set
4. Extract the Indels from the call set
5. Determine parameters for filtering indels
6. Apply the filter to the Indel call set

------

### 1. Extract the SNPs from the call set

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T SelectVariants \ 
    -R reference.fa \ 
    -V raw_variants.vcf \ 
    -selectType SNP \ 
    -o raw_snps.vcf 
```

#### Expected Result

This creates a VCF file called `raw_snps.vcf`, containing just the SNPs from the original file of raw variants.

------

### 2. Determine parameters for filtering SNPs

SNPs matching any of these conditions will be considered bad and filtered out, *i.e.* marked `FILTER` in the output VCF file. The program will specify which parameter was chiefly responsible for the exclusion of the SNP using the culprit annotation. SNPs that do not match any of these conditions will be considered good and marked `PASS` in the output VCF file.

- QualByDepth (QD) 2.0

This is the variant confidence (from the `QUAL` field) divided by the unfiltered depth of non-reference samples.

- FisherStrand (FS) 60.0

Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.

- RMSMappingQuality (MQ) 40.0

This is the Root Mean Square of the mapping quality of the reads across all samples.

- MappingQualityRankSumTest (MQRankSum) -12.5

This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele). Note that the mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, *i.e.* this will only be applied to heterozygous calls.

- ReadPosRankSumTest (ReadPosRankSum) -8.0

This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, *i.e.* this will only be applied to heterozygous calls.

- StrandOddsRatio (SOR) 3.0

The StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. Higher values indicate more strand bias.

------

### 3. Apply the filter to the SNP call set

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T VariantFiltration \ 
    -R reference.fa \ 
    -V raw_snps.vcf \ 
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \ 
    --filterName "my_snp_filter" \ 
    -o filtered_snps.vcf 
```

#### Expected Result

This creates a VCF file called `filtered_snps.vcf`, containing all the original SNPs from the `raw_snps.vcf` file, but now the SNPs are annotated with either `PASS` or `FILTER` depending on whether or not they passed the filters.

For SNPs that failed the filter, the variant annotation also includes the name of the filter. That way, if you apply several different filters (simultaneously or sequentially), you can keep track of which filter(s) each SNP failed, and later you can retrieve specific subsets of your calls using the SelectVariants tool. To learn more about composing different types of filtering expressions and retrieving subsets of variants using SelectVariants, please see the online GATK documentation.

------

### 4. Extract the Indels from the call set

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T SelectVariants \ 
    -R reference.fa \ 
    -V raw_HC_variants.vcf \ 
    -selectType INDEL \ 
    -o raw_indels.vcf 
```

#### Expected Result

This creates a VCF file called `raw_indels.vcf`, containing just the Indels from the original file of raw variants.

------

### 5. Determine parameters for filtering Indels.

Indels matching any of these conditions will be considered bad and filtered out, *i.e.* marked `FILTER` in the output VCF file. The program will specify which parameter was chiefly responsible for the exclusion of the indel using the culprit annotation. Indels that do not match any of these conditions will be considered good and marked `PASS` in the output VCF file.

- QualByDepth (QD) 2.0

This is the variant confidence (from the `QUAL` field) divided by the unfiltered depth of non-reference samples.

- FisherStrand (FS) 200.0

Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.

- ReadPosRankSumTest (ReadPosRankSum) 20.0

This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, *i.e.* this will only be applied to heterozygous calls.

- StrandOddsRatio (SOR) 10.0

The StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. Higher values indicate more strand bias.

------

### 6. Apply the filter to the Indel call set

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T VariantFiltration \ 
    -R reference.fa \ 
    -V raw_indels.vcf \ 
    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \ 
    --filterName "my_indel_filter" \ 
    -o filtered_indels.vcf 
```

#### Expected Result

This creates a VCF file called `filtered_indels.vcf`, containing all the original Indels from the `raw_indels.vcf` file, but now the Indels are annotated with either `PASS` or `FILTER` depending on whether or not they passed the filters.

For Indels that failed the filter, the variant annotation also includes the name of the filter. That way, if you apply several different filters (simultaneously or sequentially), you can keep track of which filter(s) each Indel failed, and later you can retrieve specific subsets of your calls using the SelectVariants tool. To learn more about composing different types of filtering expressions and retrieving subsets of variants using SelectVariants, please see the online GATK documentation.