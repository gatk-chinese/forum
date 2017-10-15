# 如何用HC找变异

原文链接：[ (howto) Call variants with HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/article?id=2803)

这个内容比较单一，就是翻译一个命令的使用

#### Objective

Call variants on a single genome with the HaplotypeCaller, producing a raw (unfiltered) VCF.

#### Caveat

This is meant only for single-sample analysis. To analyze multiple samples, see the Best Practices documentation on joint analysis.

#### Prerequisites

- TBD

#### Steps

1. Determine the basic parameters of the analysis
2. Call variants in your sequence data

------

### 1. Determine the basic parameters of the analysis

If you do not specify these parameters yourself, the program will use default values. However we recommend that you set them explicitly because it will help you understand how the results are bounded and how you can modify the program's behavior.

- Genotyping mode (`--genotyping_mode`)

This specifies how we want the program to determine the alternate alleles to use for genotyping. In the default `DISCOVERY` mode, the program will choose the most likely alleles out of those it sees in the data. In `GENOTYPE_GIVEN_ALLELES` mode, the program will only use the alleles passed in from a VCF file (using the `-alleles` argument). This is useful if you just want to determine if a sample has a specific genotype of interest and you are not interested in other alleles.

- Emission confidence threshold (`-stand_emit_conf`)

This is the minimum confidence threshold (phred-scaled) at which the program should emit sites that appear to be possibly variant.

- Calling confidence threshold (`-stand_call_conf`)

This is the minimum confidence threshold (phred-scaled) at which the program should emit variant sites as called. If a site's associated genotype has a confidence score lower than the calling threshold, the program will emit the site as filtered and will annotate it as LowQual. This threshold separates high confidence calls from low confidence calls.

*The terms "called" and "filtered" are tricky because they can mean different things depending on context. In ordinary language, people often say a site was called if it was emitted as variant. But in the GATK's technical language, saying a site was called means that that site passed the confidence threshold test. For filtered, it's even more confusing, because in ordinary language, when people say that sites were filtered, they usually mean that those sites successfully passed a filtering test. However, in the GATK's technical language, the same phrase (saying that sites were filtered) means that those sites failed the filtering test. In effect, it means that those would be filtered out if the filter was used to actually remove low-confidence calls from the callset, instead of just tagging them. In both cases, both usages are valid depending on the point of view of the person who is reporting the results. So it's always important to check what is the context when interpreting results that include these terms.*

------

### 2. Call variants in your sequence data

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T HaplotypeCaller \ 
    -R reference.fa \ 
    -I preprocessed_reads.bam \  
    -L 20 \ 
    --genotyping_mode DISCOVERY \ 
    -stand_emit_conf 10 \ 
    -stand_call_conf 30 \ 
    -o raw_variants.vcf 
```

*Note that -L specifies that we only want to run the command on a subset of the data (here, chromosome 20). This is useful for testing as well as other purposes, as documented here. For example, when running on exome data, we use -L to specify a file containing the list of exome targets corresponding to the capture kit that was used to generate the exome libraries.*

#### Expected Result

This creates a VCF file called `raw_variants.vcf`, containing all the sites that the HaplotypeCaller evaluated to be potentially variant. Note that this file contains both SNPs and Indels.

Although you now have a nice fresh set of variant calls, the variant discovery stage is not over. The distinctions made by the caller itself between low-confidence calls and the rest is very primitive, and should not be taken as a definitive guide for filtering. The GATK callers are designed to be very lenient in calling variants, so it is extremely important to apply one of the recommended filtering methods (variant recalibration or hard-filtering), in order to move on to downstream analyses with the highest-quality call set possible.



