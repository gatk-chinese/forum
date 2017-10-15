# 如何运行BQSR

原文链接: [ (howto) Recalibrate base quality scores = run BQSR](https://software.broadinstitute.org/gatk/documentation/article?id=2801)

内容不多，需要翻译

#### Objective

Recalibrate base quality scores in order to correct sequencing errors and other experimental artifacts.

#### Prerequisites

- TBD

#### Steps

1. Analyze patterns of covariation in the sequence dataset
2. Do a second pass to analyze covariation remaining after recalibration
3. Generate before/after plots
4. Apply the recalibration to your sequence data

------

### 1. Analyze patterns of covariation in the sequence dataset

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T BaseRecalibrator \ 
    -R reference.fa \ 
    -I input_reads.bam \ 
    -L 20 \ 
    -knownSites dbsnp.vcf \ 
    -knownSites gold_indels.vcf \ 
    -o recal_data.table 
```

#### Expected Result

This creates a GATKReport file called `recal_data.table` containing several tables. These tables contain the covariation data that will be used in a later step to recalibrate the base qualities of your sequence data.

It is imperative that you provide the program with a set of known sites, otherwise it will refuse to run. The known sites are used to build the covariation model and estimate empirical base qualities. For details on what to do if there are no known sites available for your organism of study, please see the online GATK documentation.

Note that `-L 20` is used here and in the next steps to restrict analysis to only chromosome 20 in the b37 human genome reference build. To run against a different reference, you may need to change the name of the contig according to the nomenclature used in your reference.

------

### 2. Do a second pass to analyze covariation remaining after recalibration

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T BaseRecalibrator \ 
    -R reference.fa \ 
    -I input_reads.bam \ 
    -L 20 \ 
    -knownSites dbsnp.vcf \ 
    -knownSites gold_indels.vcf \ 
    -BQSR recal_data.table \ 
    -o post_recal_data.table 
```

#### Expected Result

This creates another GATKReport file, which we will use in the next step to generate plots. Note the use of the `-BQSR` flag, which tells the GATK engine to perform on-the-fly recalibration based on the first recalibration data table.

------

### 3. Generate before/after plots

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T AnalyzeCovariates \ 
    -R reference.fa \ 
    -L 20 \ 
    -before recal_data.table \
    -after post_recal_data.table \
    -plots recalibration_plots.pdf
```

#### Expected Result

This generates a document called `recalibration_plots.pdf` containing plots that show how the reported base qualities match up to the empirical qualities calculated by the BaseRecalibrator. Comparing the **before** and **after**plots allows you to check the effect of the base recalibration process before you actually apply the recalibration to your sequence data. For details on how to interpret the base recalibration plots, please see the online GATK documentation.

------

### 4. Apply the recalibration to your sequence data

#### Action

Run the following GATK command:

```
java -jar GenomeAnalysisTK.jar \ 
    -T PrintReads \ 
    -R reference.fa \ 
    -I input_reads.bam \ 
    -L 20 \ 
    -BQSR recal_data.table \ 
    -o recal_reads.bam 
```

#### Expected Result

This creates a file called `recal_reads.bam` containing all the original reads, but now with exquisitely accurate base substitution, insertion and deletion quality scores. By default, the original quality scores are discarded in order to keep the file size down. However, you have the option to retain them by adding the flag `–emit_original_quals` to the PrintReads command, in which case the original qualities will also be written in the file, tagged `OQ`.

Notice how this step uses a very simple tool, PrintReads, to apply the recalibration. What’s happening here is that we are loading in the original sequence data, having the GATK engine recalibrate the base qualities on-the-fly thanks to the `-BQSR` flag (as explained earlier), and just using PrintReads to write out the resulting data to the new file.