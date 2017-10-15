# 如何从fasta序列里面模拟测序reads？

原文链接 [ (How to) Simulate reads using a reference genome ALT contig](https://software.broadinstitute.org/gatk/documentation/article?id=7859)

这个功能比较好理解，都可以自己写脚本来实现了，这里主要介绍两个外部工具的使用，讲清楚即可。



This tutorial shows how to generate simulated reads against a specific target sequence. This can be useful, e.g. if you want to simulate reads for an alternate contig in GRCh38/hg38 to see how they end up mapping between the primary assembly versus the alternate contig.

We use external tools to accomplish this. In **Section 1**, we use [Samtools](http://www.htslib.org/) to subset the target contig sequence from a reference FASTA file. In **Section 2**, we use [wgsim](https://github.com/lh3/wgsim) to generate FASTQ format paired reads against the target contig. The resulting read data is ready for alignment.

This tutorial provides example data for you to follow along and includes a mini-reference FASTA. If you are unfamiliar with terms that describe reference genome components, take a few minutes to study the *Dictionary* entry [Reference Genome Components](http://gatkforums.broadinstitute.org/dsde/discussion/7857).

