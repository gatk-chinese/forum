# 如何不比对就把fastq文件转为bam格式呢？

原文链接： [(How to) Generate an unmapped BAM (uBAM)](http://gatkforums.broadinstitute.org/discussion/6484/).

不需要全文翻译，只需要说清楚这个需求是为什么，然后讲一下两种方法即可。

[![image](https://us.v-cdn.net/5019796/uploads/FileUpload/31/992f952c8e9819d57bf74b0a4ac308.png)](https://us.v-cdn.net/5019796/uploads/FileUpload/31/992f952c8e9819d57bf74b0a4ac308.png)Here we outline how to generate an unmapped BAM (uBAM) from either a FASTQ or aligned BAM file. We use Picard's FastqToSam to convert a FASTQ (**Option A**) or Picard's RevertSam to convert an aligned BAM (**Option B**).

#### Jump to a section on this page

(A) [Convert FASTQ to uBAM](https://gatkforums.broadinstitute.org/gatk/discussion/6484#optionA) and add read group information using FastqToSam
(B) [Convert aligned BAM to uBAM](https://gatkforums.broadinstitute.org/gatk/discussion/6484#optionB) and discard problematic records using RevertSam

#### Tools involved

- [FastqToSam](https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam)
- [RevertSam](https://broadinstitute.github.io/picard/command-line-overview.html#RevertSam)

#### Prerequisites

- Installed Picard tools

#### Download example data

- [tutorial_6484_FastqToSam.tar.gz](https://drive.google.com/open?id=0BzI1CyccGsZiUXVNT3hsNldvUFk)
- [tutorial_6484_RevertSam.tar.gz](https://drive.google.com/open?id=0BzI1CyccGsZiMWZacmVWWnV2VFE)

Tutorial data reads were originally aligned to the [advanced tutorial bundle](http://gatkforums.broadinstitute.org/discussion/4610/)'s human_g1k_v37_decoy.fasta reference and to 10:91,000,000-92,000,000.

#### Related resources

- Our [dictionary entry on read groups](http://gatkforums.broadinstitute.org/discussion/6472/read-groups#latest) discusses the importance of assigning appropriate read group fields that differentiate samples and factors that contribute to batch effects, e.g. flow cell lane. Be sure your read group fields conform to the recommendations.
- [This post](http://gatkforums.broadinstitute.org/discussion/2909#addRG) provides an example command for AddOrReplaceReadGroups.
- This *How to* is part of a larger workflow and tutorial on [(How to) Efficiently map and clean up short read sequence data](http://gatkforums.broadinstitute.org/discussion/6483).
- To extract reads in a genomic interval from the aligned BAM, [use GATK's PrintReads](http://gatkforums.broadinstitute.org/discussion/6517/).

------

### (A) Convert FASTQ to uBAM and add read group information using FastqToSam

Picard's [FastqToSam](https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam) transforms a FASTQ file to an unmapped BAM, requires two read group fields and makes optional specification of other read group fields. In the command below we note which fields are required for GATK Best Practices Workflows. All other read group fields are optional.

```
java -Xmx8G -jar picard.jar FastqToSam \
    FASTQ=6484_snippet_1.fastq \ #first read file of pair
    FASTQ2=6484_snippet_2.fastq \ #second read file of pair
    OUTPUT=6484_snippet_fastqtosam.bam \
    READ_GROUP_NAME=H0164.2 \ #required; changed from default of A
    SAMPLE_NAME=NA12878 \ #required
    LIBRARY_NAME=Solexa-272222 \ #required 
    PLATFORM_UNIT=H0164ALXX140820.2 \ 
    PLATFORM=illumina \ #recommended
    SEQUENCING_CENTER=BI \ 
    RUN_DATE=2014-08-20T00:00:00-0400

```

Some details on select parameters:

- For paired reads, specify each FASTQ file with `FASTQ` and `FASTQ2` for the first read file and the second read file, respectively. Records in each file must be queryname sorted as the tool assumes identical ordering for pairs. The tool automatically strips the `/1` and `/2` read name suffixes and adds SAM flag values to indicate reads are paired. Do not provide a single interleaved fastq file, as the tool will assume reads are unpaired and the SAM flag values will reflect single ended reads.
- For single ended reads, specify the input file with `FASTQ`.
- `QUALITY_FORMAT` is detected automatically if unspecified.
- `SORT_ORDER` by default is queryname.
- `PLATFORM_UNIT` is often in run_barcode.lane format. Include if sample is multiplexed.
- `RUN_DATE` is in [Iso8601 date format](https://en.wikipedia.org/wiki/ISO_8601).

Paired reads will have [SAM flag](https://broadinstitute.github.io/picard/explain-flags.html) values that reflect pairing and the fact that the reads are unmapped as shown in the example read pair below.

**Original first read**

```
@H0164ALXX140820:2:1101:10003:49022/1
ACTTTAGAAATTTACTTTTAAGGACTTTTGGTTATGCTGCAGATAAGAAATATTCTTTTTTTCTCCTATGTCAGTATCCCCCATTGAAATGACAATAACCTAATTATAAATAAGAATTAGGCTTTTTTTTGAACAGTTACTAGCCTATAGA
+
-FFFFFJJJJFFAFFJFJJFJJJFJFJFJJJ<<FJJJJFJFJFJJJJ<JAJFJJFJJJJJFJJJAJJJJJJFFJFJFJJFJJFFJJJFJJJFJJFJJFJAJJJJAJFJJJJJFFJJ<<<JFJJAFJAAJJJFFFFFJJJAJJJF<AJFFFJ

```

**Original second read**

```
@H0164ALXX140820:2:1101:10003:49022/2
TGAGGATCACTAGATGGGGGAGGGAGAGAAGAGATGTGGGCTGAAGAACCATCTGTTGGGTAATATGTTTACTGTCAGTGTGATGGAATAGCTGGGACCCCAAGCGTCAGTGTTACACAACTTACATCTGTTGATCGACTGTCTATGACAG
+
AA<FFJJJAJFJFAFJJJJFAJJJJJ7FFJJ<F-FJFJJJFJJFJJFJJF<FJJA<JF-AFJFAJFJJJJJAAAFJJJJJFJJF-FF<7FJJJJJJ-JA<<J<F7-<FJFJJ7AJAF-AFFFJA--J-F######################

```

**After FastqToSam**

```
H0164ALXX140820:2:1101:10003:49022      77      *       0       0       *       *       0       0       ACTTTAGAAATTTACTTTTAAGGACTTTTGGTTATGCTGCAGATAAGAAATATTCTTTTTTTCTCCTATGTCAGTATCCCCCATTGAAATGACAATAACCTAATTATAAATAAGAATTAGGCTTTTTTTTGAACAGTTACTAGCCTATAGA -FFFFFJJJJFFAFFJFJJFJJJFJFJFJJJ<<FJJJJFJFJFJJJJ<JAJFJJFJJJJJFJJJAJJJJJJFFJFJFJJFJJFFJJJFJJJFJJFJJFJAJJJJAJFJJJJJFFJJ<<<JFJJAFJAAJJJFFFFFJJJAJJJF<AJFFFJ RG:Z:H0164.2
H0164ALXX140820:2:1101:10003:49022      141     *       0       0       *       *       0       0       TGAGGATCACTAGATGGGGGAGGGAGAGAAGAGATGTGGGCTGAAGAACCATCTGTTGGGTAATATGTTTACTGTCAGTGTGATGGAATAGCTGGGACCCCAAGCGTCAGTGTTACACAACTTACATCTGTTGATCGACTGTCTATGACAG AA<FFJJJAJFJFAFJJJJFAJJJJJ7FFJJ<F-FJFJJJFJJFJJFJJF<FJJA<JF-AFJFAJFJJJJJAAAFJJJJJFJJF-FF<7FJJJJJJ-JA<<J<F7-<FJFJJ7AJAF-AFFFJA--J-F###################### RG:Z:H0164.2

```

[back to top](https://gatkforums.broadinstitute.org/gatk/discussion/6484#top)

------

### (B) Convert aligned BAM to uBAM and discard problematic records using RevertSam

We use Picard's [RevertSam](https://broadinstitute.github.io/picard/command-line-overview.html#RevertSam) to remove alignment information and generate an unmapped BAM (uBAM). For our tutorial file we have to call on some additional parameters that we explain below. This illustrates the need to cater the tool's parameters to each dataset. As such, it is a good idea to test the reversion process on a subset of reads before committing to reverting the entirety of a large BAM. Follow the directions in [this *How to*](http://gatkforums.broadinstitute.org/discussion/6517/) to create a snippet of aligned reads corresponding to a genomic interval.

We use the following parameters.

```
java -Xmx8G -jar /path/picard.jar RevertSam \
    I=6484_snippet.bam \
    O=6484_snippet_revertsam.bam \
    SANITIZE=true \ 
    MAX_DISCARD_FRACTION=0.005 \ #informational; does not affect processing
    ATTRIBUTE_TO_CLEAR=XT \
    ATTRIBUTE_TO_CLEAR=XN \
    ATTRIBUTE_TO_CLEAR=AS \ #Picard release of 9/2015 clears AS by default
    ATTRIBUTE_TO_CLEAR=OC \
    ATTRIBUTE_TO_CLEAR=OP \
    SORT_ORDER=queryname \ #default
    RESTORE_ORIGINAL_QUALITIES=true \ #default
    REMOVE_DUPLICATE_INFORMATION=true \ #default
    REMOVE_ALIGNMENT_INFORMATION=true #default

```

To process large files, also designate a temporary directory.

```
    TMP_DIR=/path/shlee #sets environmental variable for temporary directory

```

**We invoke or change multiple RevertSam parameters to generate an unmapped BAM**

- We remove nonstandard alignment tags with the `ATTRIBUTE_TO_CLEAR` option. Standard tags cleared by default are NM, UQ, PG, MD, MQ, SA, MC, and AS tags (AS for Picard releases starting 9/2015). Additionally, the OQ tag is removed by the default `RESTORE_ORIGINAL_QUALITIES` parameter. Remove all other nonstandard tags by specifying each with the `ATTRIBUTE_TO_CLEAR` option. For example, we clear the `XT` tag using this option for our tutorial file so that it is free for use by other tools, e.g. MarkIlluminaAdapters. To list all tags within a BAM, use the command below.

  ```
  samtools view input.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '{ if(!x[$1]++) { print }}' 

  ```

  For the tutorial file, this gives RG, OC, XN, OP and XT tags as well as those removed by default. We remove all of these except the RG tag. See your aligner's documentation and the [Sequence Alignment/Map Format Specification](http://samtools.sourceforge.net/SAM1.pdf) for descriptions of tags.

- Additionally, we invoke the `SANITIZE` option to remove reads that cause problems for certain tools, e.g. MarkIlluminaAdapters. Downstream tools will have problems with paired reads with missing mates, duplicated records, and records with mismatches in length of bases and qualities. Any paired reads file subset for a genomic interval requires sanitizing to remove reads with lost mates that align outside of the interval.

- In this command, we've set `MAX_DISCARD_FRACTION` to a more strict threshold of 0.005 instead of the default 0.01. Whether or not this fraction is reached, the tool informs you of the number and fraction of reads it discards. This parameter asks the tool to additionally inform you of the discarded fraction via an exception as it finishes processing.

  ```
  Exception in thread "main" picard.PicardException: Discarded 0.787% which is above MAX_DISCARD_FRACTION of 0.500%  

  ```

**Some comments on options kept at default:**

- `SORT_ORDER`=queryname
  For paired read files, because each read in a pair has the same query name, sorting results in interleaved reads. This means that reads in a pair are listed consecutively within the same file. We make sure to alter the previous sort order. Coordinate sorted reads result in the aligner incorrectly estimating insert size from blocks of paired reads as they are not randomly distributed.
- `RESTORE_ORIGINAL_QUALITIES`=true
  Restoring original base qualities to the QUAL field requires OQ tags listing original qualities. The OQ tag uses the same encoding as the QUAL field, e.g. ASCII Phred-scaled base quality+33 for tutorial data. After restoring the QUAL field, RevertSam removes the tag.
- `REMOVE_ALIGNMENT_INFORMATION`=true will remove program group records and alignment flag and tag information. For example, [flags](https://broadinstitute.github.io/picard/explain-flags.html) reset to unmapped values, e.g. 77 and 141 for paired reads. The parameter also invokes the default `ATTRIBUTE_TO_CLEAR` parameter which removes standard alignment tags. RevertSam ignores `ATTRIBUTE_TO_CLEAR` when `REMOVE_ALIGNMENT_INFORMATION`=false.

Below we show below a read pair before and after RevertSam from the tutorial data. Notice the first listed read in the pair becomes **reverse-complemented** after RevertSam. This restores how reads are represented when they come off the sequencer--5' to 3' of the read being sequenced.

For 6484_snippet.bam, `SANITIZE` removes 2,202 out of 279,796 (0.787%) reads, leaving us with 277,594 reads.

**Original BAM**

```
H0164ALXX140820:2:1101:10003:23460  83  10  91515318    60  151M    =   91515130    -339    CCCATCCCCTTCCCCTTCCCTTTCCCTTTCCCTTTTCTTTCCTCTTTTAAAGAGACAAGGTCTTGTTCTGTCACCCAGGCTGGAATGCAGTGGTGCAGTCATGGCTCACTGCCGCCTCAGACTTCAGGGCAAAAGCAATCTTTCCAGCTCA :<<=>@AAB@AA@AA>6@@A:>,*@A@<@??@8?9>@==8?:?@?;?:><??@>==9?>8>@:?>>=>;<==>>;>?=?>>=<==>>=>9<=>??>?>;8>?><?<=:>>>;4>=>7=6>=>>=><;=;>===?=>=>>?9>>>>??==== MC:Z:60M91S MD:Z:151    PG:Z:MarkDuplicates RG:Z:H0164.2    NM:i:0  MQ:i:0  OQ:Z:<FJFFJJJJFJJJJJF7JJJ<F--JJJFJJJJ<J<FJFF<JAJJJAJAJFFJJJFJAFJAJJAJJJJJFJJJJJFJJFJJJJFJFJJJJFFJJJJJJJFAJJJFJFJFJJJFFJJJ<J7JJJJFJ<AFAJJJJJFJJJJJAJFJJAFFFFA    UQ:i:0  AS:i:151

H0164ALXX140820:2:1101:10003:23460  163 10  91515130    0   60M91S  =   91515318    339 TCTTTCCTTCCTTCCTTCCTTGCTCCCTCCCTCCCTCCTTTCCTTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTCCCCTCTCCCACCCCTCTCTCCCCCCCTCCCACCC :0;.=;8?7==?794<<;:>769=,<;0:=<0=:9===/,:-==29>;,5,98=599;<=########################################################################################### SA:Z:2,33141573,-,37S69M45S,0,1;    MC:Z:151M   MD:Z:48T4T6 PG:Z:MarkDuplicates RG:Z:H0164.2    NM:i:2  MQ:i:60 OQ:Z:<-<-FA<F<FJF<A7AFAAJ<<AA-FF-AJF-FA<AFF--A-FA7AJA-7-A<F7<<AFF###########################################################################################    UQ:i:49 AS:i:50

```

**After RevertSam**

```
H0164ALXX140820:2:1101:10003:23460  77  *   0   0   *   *   0   0   TGAGCTGGAAAGATTGCTTTTGCCCTGAAGTCTGAGGCGGCAGTGAGCCATGACTGCACCACTGCATTCCAGCCTGGGTGACAGAACAAGACCTTGTCTCTTTAAAAGAGGAAAGAAAAGGGAAAGGGAAAGGGAAGGGGAAGGGGATGGG AFFFFAJJFJAJJJJJFJJJJJAFA<JFJJJJ7J<JJJFFJJJFJFJFJJJAFJJJJJJJFFJJJJFJFJJJJFJJFJJJJJFJJJJJAJJAJFAJFJJJFFJAJAJJJAJ<FFJF<J<JJJJFJJJ--F<JJJ7FJJJJJFJJJJFFJF< RG:Z:H0164.2

H0164ALXX140820:2:1101:10003:23460  141 *   0   0   *   *   0   0   TCTTTCCTTCCTTCCTTCCTTGCTCCCTCCCTCCCTCCTTTCCTTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTCCCCTCTCCCACCCCTCTCTCCCCCCCTCCCACCC <-<-FA<F<FJF<A7AFAAJ<<AA-FF-AJF-FA<AFF--A-FA7AJA-7-A<F7<<AFF########################################################################################### RG:Z:H0164.2

```

[](https://gatkforums.broadinstitute.org/gatk/discussion/6484#top)



 