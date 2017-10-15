# GATK软件需要的输入文件是哪些？

原文链接：[ What input files does the GATK accept / require? ](https://software.broadinstitute.org/gatk/documentation/article?id=1204)

All analyses done with the GATK typically involve several (though not necessarily all) of the following inputs:

- Reference genome sequence
- Sequencing reads
- Intervals of interest
- Reference-ordered data

This article describes the corresponding file formats that are acceptable for use with the GATK.

------

### 1. Reference Genome Sequence

The GATK requires the reference sequence in a single reference sequence in FASTA format, with all contigs in the same file. The GATK requires strict adherence to the FASTA standard. All the standard IUPAC bases are accepted, but keep in mind that non-standard bases (i.e. other than ACGT, such as W for example) will be ignored (i.e. those positions in the genome will be skipped).

**Some users have reported having issues with reference files that have been stored or modified on Windows filesystems. The issues manifest as "10" characters (corresponding to encoded newlines) inserted in the sequence, which cause the GATK to quit with an error. If you encounter this issue, you will need to re-download a valid master copy of the reference file, or clean it up yourself.**

Gzipped fasta files will not work with the GATK, so please make sure to unzip them first. Please see [this article](http://www.broadinstitute.org/gatk/guide/article?id=1601) for more information on preparing FASTA reference sequences for use with the GATK.

#### Important note about human genome reference versions

If you are using human data, your reads must be aligned to one of the official b3x (e.g. b36, b37) or hg1x (e.g. hg18, hg19) references. The names and order of the contigs in the reference you used must exactly match that of one of the official references canonical orderings. These are defined by historical karotyping of largest to smallest chromosomes, followed by the X, Y, and MT for the b3x references; the order is thus 1, 2, 3, ..., 10, 11, 12, ... 20, 21, 22, X, Y, MT. The hg1x references differ in that the chromosome names are prefixed with "chr" and chrM appears first instead of last. The GATK will detect misordered contigs (for example, lexicographically sorted) and throw an error. This draconian approach, though unnecessary technically, ensures that all supplementary data provided with the GATK works correctly. You can use ReorderSam to fix a BAM file aligned to a missorted reference sequence.

**Our Best Practice recommendation is that you use a standard GATK reference from the GATK resource bundle.**

------

### 2. Sequencing Reads

The only input format for sequence reads that the GATK itself supports is the [Sequence Alignment/Map (SAM)] format. See [SAM/BAM] for more details on the SAM/BAM format as well as [Samtools](http://samtools.sourceforge.net/) and [Picard](http://picard.sourceforge.net/), two complementary sets of utilities for working with SAM/BAM files.

If you don't find the information you need in this section, please see our [FAQs on BAM files](http://www.broadinstitute.org/gatk/guide/article?id=1317).

If you are starting out your pipeline with raw reads (typically in FASTQ format) you'll need to make sure that when you map those reads to the reference and produce a BAM file, the resulting BAM file is fully compliant with the GATK requirements. See the Best Practices documentation for detailed instructions on how to do this.

In addition to being in SAM format, we require the following additional constraints in order to use your file with the GATK:

- The file must be binary (with `.bam` file extension).
- The file must be indexed.
- The file must be sorted in coordinate order with respect to the reference (i.e. the contig ordering in your bam must exactly match that of the reference you are using).
- The file must have a proper bam header with read groups. Each read group must contain the platform (PL) and sample (SM) tags. For the platform value, we currently support 454, LS454, Illumina, Solid, ABI_Solid, and CG (all case-insensitive).
- Each read in the file must be associated with exactly one read group.

Below is an example well-formed SAM field header and fields (with @SQ dictionary truncated to show only the first two chromosomes for brevity):

```
@HD     VN:1.0  GO:none SO:coordinate
@SQ     SN:1    LN:249250621    AS:NCBI37       UR:file:/lustre/scratch102/projects/g1k/ref/main_project/human_g1k_v37.fasta    M5:1b22b98cdeb4a9304cb5d48026a85128
@SQ     SN:2    LN:243199373    AS:NCBI37       UR:file:/lustre/scratch102/projects/g1k/ref/main_project/human_g1k_v37.fasta    M5:a0d9851da00400dec1098a9255ac712e
@RG     ID:ERR000162    PL:ILLUMINA     LB:g1k-sc-NA12776-CEU-1 PI:200  DS:SRP000031    SM:NA12776      CN:SC
@RG     ID:ERR000252    PL:ILLUMINA     LB:g1k-sc-NA12776-CEU-1 PI:200  DS:SRP000031    SM:NA12776      CN:SC
@RG     ID:ERR001684    PL:ILLUMINA     LB:g1k-sc-NA12776-CEU-1 PI:200  DS:SRP000031    SM:NA12776      CN:SC
@RG     ID:ERR001685    PL:ILLUMINA     LB:g1k-sc-NA12776-CEU-1 PI:200  DS:SRP000031    SM:NA12776      CN:SC
@PG     ID:GATK TableRecalibration      VN:v2.2.16      CL:Covariates=[ReadGroupCovariate, QualityScoreCovariate, DinucCovariate, CycleCovariate], use_original_quals=true, defau 
t_read_group=DefaultReadGroup, default_platform=Illumina, force_read_group=null, force_platform=null, solid_recal_mode=SET_Q_ZERO, window_size_nqs=5, homopolymer_nback=7, except on_if_no_tile=false, pQ=5, maxQ=40, smoothing=137       UR:file:/lustre/scratch102/projects/g1k/ref/main_project/human_g1k_v37.fasta    M5:b4eb71ee878d3706246b7c1dbef69299
@PG     ID:bwa  VN:0.5.5
ERR001685.4315085       16      1       9997    25      35M     *       0       0       CCGATCTCCCTAACCCTAACCCTAACCCTAACCCT     ?8:C7ACAABBCBAAB?CCAABBEBA@ACEBBB@?     XT:A:U  XN:i:4    X0:i:1  X1:i:0  XM:i:2  XO:i:0  XG:i:0  RG:Z:ERR001685  NM:i:6  MD:Z:0N0N0N0N1A0A28     OQ:Z:>>:>2>>>>>>>>>>>>>>>>>>?>>>>??>???>
ERR001689.1165834       117     1       9997    0       *       =       9997    0       CCGATCTAGGGTTAGGGTTAGGGTTAGGGTTAGGG     >7AA<@@C?@?B?B??>9?B??>A?B???BAB??@     RG:Z:ERR001689    OQ:Z:>:<<8<<<><<><><<>7<>>>?>>??>???????
ERR001689.1165834       185     1       9997    25      35M     =       9997    0       CCGATCTCCCTAACCCTAACCCTAACCCTAACCCT     758A:?>>8?=@@>>?;4<>=??@@==??@?==?8     XT:A:U  XN:i:4    SM:i:25 AM:i:0  X0:i:1  X1:i:0  XM:i:2  XO:i:0  XG:i:0  RG:Z:ERR001689  NM:i:6  MD:Z:0N0N0N0N1A0A28     OQ:Z:;74>7><><><>>>>><:<>>>>>>>>>>>>>>>>
ERR001688.2681347       117     1       9998    0       *       =       9998    0       CGATCTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG     5@BA@A6B???A?B??>B@B??>B@B??>BAB???     RG:Z:ERR001688    OQ:Z:=>>>><4><<?><??????????????????????       
```

#### Note about fixing BAM files with alternative sortings

The GATK requires that the BAM file be sorted in the same order as the reference. Unfortunately, many BAM files have headers that are sorted in some other order -- lexicographical order is a common alternative. To resort the BAM file please use [ReorderSam](http://picard.sourceforge.net/command-line-overview.shtml#ReorderSam).

------

### 3. Intervals of interest

The GATK accept interval files for processing subsets of the genome in several different formats. Please see the [FAQs on interval lists](http://www.broadinstitute.org/gatk/guide/article?id=1319) for details.

------

### 4. Reference Ordered Data (ROD) file formats

The GATK can associate arbitrary reference ordered data (ROD) files with named tracks for all tools. Some tools require specific ROD data files for processing, and developers are free to write tools that access arbitrary data sets using the ROD interface. The general ROD system has the following syntax:

```
-argumentName:name,type file
```

Where `name` is the name in the GATK tool (like "eval" in VariantEval), `type` is the type of the file, such as VCF or dbSNP, and `file` is the path to the file containing the ROD data.

The GATK supports several common file formats for reading ROD data:

- [VCF](http://www.1000genomes.org/wiki/analysis/variant-call-format/) : VCF type, the recommended format for representing variant loci and genotype calls. The GATK will only process valid VCF files; [VCFTools](http://vcftools.sourceforge.net/) provides the official VCF validator. See [here](http://vcftools.sourceforge.net/VCF-poster.pdf) for a useful poster detailing the VCF specification.
- UCSC formated dbSNP : dbSNP type, UCSC dbSNP database output
- BED : BED type, a general purpose format for representing genomic interval data, useful for masks and other interval outputs. **Please note that the bed format is 0-based while most other formats are 1-based.**

**Note that we no longer support the PED format. See here for converting .ped files to VCF.**

If you need additional information on VCF files, please see our FAQs on VCF files [here](http://www.broadinstitute.org/gatk/guide/article?id=1318) and [here](http://www.broadinstitute.org/gatk/guide/article?id=1268).