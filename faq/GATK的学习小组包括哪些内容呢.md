# GATK的研讨会包括哪些内容呢？

原文链接：[ What do GATK workshops cover?](https://software.broadinstitute.org/gatk/documentation/article?id=9218)

## 为期三天的GATK研讨会

这次研讨会将会集中在使用broad研究所开发的GATK软件来寻找变异位点的核心步骤，即GATK小组开发的最佳实践流程。你将会了解为什么最佳实践的每个步骤都是寻找变异位点这一过程必不可少的，每一步骤分别对数据做了什么，如何使用GATK工具来处理你自己的数据并得到最可信最准确的结果。同时，在研讨会的课程中，我们会突出讲解部分GATK的新功能，比如用GVCF模式来对多样本联合寻找变异，还有如何对RNA-seq数据特异的寻找变异流程，和如何使用mutech2来选择somatic变异。我们也会介绍一下即将到来的GATK4的功能，比如它的CNV流程。最后我们会演练整个GATK流程的各个步骤的衔接和实战。

## 课程设置

本次研讨会包括一天的讲座，中间会提供大量的问答环节，力图把GATK的基础知识介绍给大家。并且伴随着2天的上机实战，具体安排如下：

- 第一天全天：演讲，关于GATK最佳实践应用于高通量测序数据的前因后果，理论基础以及应用形式。
- 第二天上午：上机操作， 使用GATK寻找变异位点，包括(SNPs + Indels)
- 第二天下午：上机操作，使用GATK对找到变异位点进行过滤。
- 第三天上午：上机操作，使用GATK寻找somatic的变异位点，包括 (SNPs + Indels + CNV)
- 第三天下午：介绍通过WDL(Workflow Description Language)整合的云端流程

上机实践侧重于数据分析步骤，希望提供练习的方式帮助用户了解寻找变异位点这一过程中涉及到的各种文件格式的转换，以及如何更好的利用GATK工具来完成这一流程。在这些实践课程中，我们会有机的结合GATK和PICARDS工具，还有其它一些第三方工具，比如Samtools, IGV, RStudio and RTG Tools ，介绍它们的一些操作技巧，比较有趣的分析方法，还有中间可能会遇到的报错信息的处理。

还会讲解如何使用WDL来写GATK分析流程的脚本，这个是broad最新开发的一个流程描述语言。并且在本地或者云端执行自己写好的流程。

请注意，本次研讨会侧重于人类研究数据的处理，本次研究会讲解的大部分操作步骤都适合其它非人类研究数据的处理，但是应用起来会有一些需要解决的移植问题，不过在这一点上面我们不会讲太多。

## 目标受众

第一天的基础知识讲解的受众比较广泛，既可以是第一次接触变异位点寻找这个分析方法，或者第一次接触GATK，想从零开始了解；也可以是已经使用过GATK的用户，想提升自己对GATK的理解。但是参会者需要熟练一些基础名词，并且有基因组或者遗传学的基本背景。

上机操作是为那些想了解GATK具体使用过程的细节新手或者中级用户准备的，需要对LINUX命令行有着基础的理解能力。希望参会者可以带自己的笔记本电脑，并且预先安装好一些必备的软件工具，具体细节将会在研究会正式开始前两周公布，除非研究会的举办方提供了练习用的电脑或者云服务器。而且必须是mac或者linux系统，Windows系统是不能参加上机实践的。限30人！

## GATK最佳实践的PPT演讲内容如下：

- Introduction to variant discovery analysis and GATK Best Practices
- (coffee break)
- Marking Duplicates
- Indel Realignment (to be phased out)
- Base Recalibration
- (lunch break)
- Germline Variant Calling and Joint Genotyping
- Filtering germline variants with VQSR
- Genotype Refinement Workflow
- (coffee break)
- Callset Evaluation
- Somatic variant discovery with MuTect2
- Preview of CNV discovery with GATK4