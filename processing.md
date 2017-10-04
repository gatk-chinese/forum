# GATK翻译任务安排

参考 https://github.com/broadinstitute/picard 的翻译策略，最终形成 http://broadinstitute.github.io/picard/ 的可视化网页界面！协作的翻译稿件以markdown形式存放在统一的github里面 https://github.com/gatk-chinese/forum

## [名词解释](https://software.broadinstitute.org/gatk/documentation/topic?name=dictionary)

一个人即可

### [官方教程](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials)

初识GATK系列：

- GATK best practices 相关软件安装以及测试
- GATK workshops 相关软件安装及测试
- Queue 以及 cloud 的安装及测试 使用
- GATK4 beta的使用
- IGV 使用

**如何使用系列：**

- mapping
- Markduplicates
- local realignment 
- BQSR(recallibrate base quality scores)
- VQSR(recallibrate variant quality scores)
- haplotypeCaller
- UnifiedGenotyper

文件格式转换

- 生成uBAM
- bam -> fastq
- Genomic interval 
- GRCh38 以及 ALT contig

**评价，注释，过滤**

- bam和vcf格式的检查
- vcf的过滤(SelectVariants)
- vcf的评价
- vcf的注释

## FAQs

FAQs是对上面教程的补充，主要是有很多人没有仔细阅读文档，或者没能充分理解，提出的常见的问题。

比如VCA,BAM,interval list等文件格式的问题。

运行GATK的准备文件，包括参考基因组，bundle 数据集。

一些参数的不理解。

一些命令的不理解，是否可以省略某些步骤

是否适用于非二倍体，是否适用于转录组数据，多样本分析注意事项

##  [用法及算法](https://software.broadinstitute.org/gatk/documentation/topic?name=methods)

首先是DNA和RNA测序的数据分析最佳实践

其次是多样本的GVCF模式

还有 HaplotypeCaller流程的详细解读

 HC overview: How the HaplotypeCaller works
 HC step 1: Defining ActiveRegions by measuring data entropy
 HC step 2: Local re-assembly and haplotype determination
 HC step 3 : Evaluating the evidence for haplotypes and variant alleles
 HC step 4: Assigning per-sample genotypes

还有变异位点的深度分析，包括过滤，质量值，评价体系，统计检验模型的讲解。
