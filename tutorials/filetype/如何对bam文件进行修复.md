# 如何对bam文件进行修复

原文链接： [ (How to) Fix a badly formatted BAM](https://software.broadinstitute.org/gatk/documentation/article?id=2909)

这一个内容比较重要

Fix a BAM that is not indexed or not sorted, has not had duplicates marked, or is lacking read group information. The options on this page are listed in order of decreasing complexity.

You may ask, is all of this really necessary? The GATK imposes strict formatting guidelines, including requiring certain [read group information](http://gatkforums.broadinstitute.org/discussion/6472/), that other software packages do not require. Although this represents a small additional processing burden upfront, the downstream benefits are numerous, including the ability to process library data individually, and significant gains in speed and parallelization options.

#### Prerequisites

- Installed Picard tools
- If indexing or marking duplicates, then coordinate sorted reads
- If coordinate sorting, then reference aligned reads
- For each read group ID, a single BAM file. If you have a multiplexed file, separate to individual files per sample.

#### Jump to a section on this page

1. [Add read groups, coordinate sort and index](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#addRG) using AddOrReplaceReadGroups
2. [Coordinate sort and index](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#sort) using SortSam
3. [Index an already coordinate-sorted BAM](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#index) using BuildBamIndex
4. [Mark duplicates](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#markduplicates) using MarkDuplicates

#### Tools involved

- [AddOrReplaceReadGroups](http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups)
- [SortSam](https://broadinstitute.github.io/picard/command-line-overview.html#SortSam)
- [BuildBamIndex](https://software.broadinstitute.org/gatk/documentation/broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex)
- [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

#### Related resources

- Our [dictionary entry on read groups](http://gatkforums.broadinstitute.org/discussion/6472/) discusses the importance of assigning appropriate read group fields that differentiate samples and factors that contribute to batch effects, e.g. flow cell lane. Be sure that your read group fields conform to the recommendations.
- [Picard's standard options](http://broadinstitute.github.io/picard/command-line-overview.html#Overview) includes adding `CREATE_INDEX` to the commands of any of its tools that produce coordinate sorted BAMs.

------

### 1. Add read groups, coordinate sort and index using AddOrReplaceReadGroups

Use Picard's [AddOrReplaceReadGroups](http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) to appropriately label read group (`@RG`) fields, coordinate sort and index a BAM file. Only the five required `@RG` fields are included in the command shown. Consider the other optional `@RG` fields for better record keeping.

```
java -jar picard.jar AddOrReplaceReadGroups \ 
    INPUT=reads.bam \ 
    OUTPUT=reads_addRG.bam \ 
    RGID=H0164.2 \ #be sure to change from default of 1
    RGLB= library1 \ 
    RGPL=illumina \ 
    RGPU=H0164ALXX140820.2 \ 
    RGSM=sample1 \ 
```

This creates a file called `reads_addRG.bam` with the same content and sorting as the input file, except the SAM record header's `@RG` line will be updated with the new information for the specified fields and each read will now have an RG tag filled with the `@RG` ID field information. Because of this repetition, the length of the `@RG` ID field contributes to file size.

To additionally coordinate sort by genomic location and create a `.bai` index, add the following options to the command.

```
    SORT_ORDER=coordinate \ 
    CREATE_INDEX=true
```

To process large files, also designate a temporary directory.

```
    TMP_DIR=/path/shlee #sets environmental variable for temporary directory
```

------

### 2. Coordinate sort and index using SortSam

Picard's [SortSam](https://broadinstitute.github.io/picard/command-line-overview.html#SortSam) both sorts and indexes and converts between SAM and BAM formats. For coordinate sorting, reads must be aligned to a reference genome.

```
java -jar picard.jar SortSam \ 
    INPUT=reads.bam \ 
    OUTPUT=reads_sorted.bam \ 
    SORT_ORDER=coordinate \
```

Concurrently index by tacking on the following option.

```
    CREATE_INDEX=true
```

This creates a file called `reads_sorted.bam` containing reads sorted by genomic location, aka coordinate, and a `.bai` index file with the same prefix as the output, e.g. `reads_sorted.bai`, within the same directory.

To process large files, also designate a temporary directory.

```
    TMP_DIR=/path/shlee #sets environmental variable for temporary directory
```

------

### 3. Index an already coordinate-sorted BAM using BuildBamIndex

Picard's [BuildBamIndex](https://software.broadinstitute.org/gatk/documentation/broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex) allows you to index a BAM that is already coordinate sorted.

```
java -jar picard.jar BuildBamIndex \ 
    INPUT=reads_sorted.bam 
```

This creates a `.bai` index file with the same prefix as the input file, e.g. `reads_sorted.bai`, within the same directory as the input file. You want to keep this default behavior as many tools require the same prefix and directory location for the pair of files. Note that Picard tools do not systematically create an index file when they output a new BAM file, whereas GATK tools will always output indexed files.

To process large files, also designate a temporary directory.

```
    TMP_DIR=/path/shlee #sets environmental variable for temporary directory
```

------

### 4. Mark duplicates using MarkDuplicates

Picard's [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) flags both PCR and optical duplicate reads with a 1024 (0x400) [SAM flag](https://broadinstitute.github.io/picard/explain-flags.html). The input BAM must be coordinate sorted.

```
java -jar picard.jar MarkDuplicates \ 
    INPUT=reads_sorted.bam \ 
    OUTPUT=reads_markdup.bam \
    METRICS_FILE=metrics.txt \
    CREATE_INDEX=true
```

This creates a file called `reads_markdup.bam` with duplicate reads marked. It also creates a file called `metrics.txt` containing duplicate read data metrics and a `.bai` index file.

To process large files, also designate a temporary directory.

```
    TMP_DIR=/path/shlee #sets environmental variable for temporary directory
```

- During sequencing, which involves PCR amplification within the sequencing machine, by a stochastic process we end up sequencing a proportion of DNA molecules that arose from the same parent insert. To be stringent in our variant discovery, GATK tools discount the duplicate reads as evidence for or against a putative variant.

- Marking duplicates is less relevant to targeted amplicon sequencing and RNA-Seq analysis.

- Optical duplicates arise from a read being sequenced twice as neighboring clusters.

  ​