# GATK workshop的变异寻找流程

原文链接： [(howto) Discover variants with GATK - A GATK Workshop Tutorial](https://software.broadinstitute.org/gatk/documentation/article?id=7869)

跟其它最佳实践是大部分重复的，应该是不需要翻译。

## GATK TUTORIAL :: Variant Discovery :: Worksheet

**June 2016 - GATK 3.6**

This tutorial covers material taught at GATK workshops, and focuses on key steps of the GATK Best Practices for Germline SNP and Indel Discovery in Whole Genomes and Exomes. If you aren't already, please set up your computer using the [workshop-specific installation instructions](https://www.broadinstitute.org/gatk/guide/article?id=7098). You can find additional background information relevant to this tutorial in the [Variant Discovery Appendix](https://www.broadinstitute.org/gatk/guide/article?id=7870).

Our main purpose is to **demonstrate an effective workflow for calling germline SNPs and indels** in cohorts of multiple samples. This workflow can be applied to **whole genomes** as well as **exomes** and other targeted sequencing datasets.

We’ll start by examining the **differences between data types** (whole genomes, exomes and RNAseq) to highlight the properties of the data that influence what we need to do to analyze it as well as what we can expect to get out of it.

Once we understand our data, we will demonstrate how **key features of the HaplotypeCaller** enable it to produce better results than position-based callers like UnifiedGenotyper. In particular, we’ll show how **local assembly of haplotypes and realignment of reads** are crucial to producing superior indel calls. Along the way we’ll show you useful tips and tricks for **troubleshooting variant calls** with HaplotypeCaller and the IGV genome browser.

All this will build up to demonstrating the **GVCF workflow for joint variant analysis**, as applied to a trio of whole-genome samples. We hope to convince you that this workflow has substantial practical advantages over a joint analysis that is achieved by calling variants simultaneously on all samples, while producing **results that are just as good** or even better.

The tutorial dataset is available for public download [here](https://drive.google.com/folderview?id=0BwTg3aXzGxEDNTF3M2hhSnBPU2s&usp=sharing).

------

### Table of Contents

1. WORKING WITH DATASETS FROM DIFFERENT EXPERIMENTAL DESIGNS 1.1 [The genome reference: b37](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#1.1)1.2 [The test sample: NA12878 Whole-Genome Sequence (WGS)](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#1.2) 1.3 [For comparison: NA12878 Exome Sequence](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#1.3) 1.4 [Another comparison: NA12878 RNAseq](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#1.4)
2. DIAGNOSING UNKNOWN BAMS 2.1 [View header and check read groups](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#2.1) 2.2 [Validate the file](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#2.2)
3. VARIANT DISCOVERY 3.1 [Call variants with a position-based caller: UnifiedGenotyper](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.1) 3.2 [Call variants with HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.2)  3.2.1 [View realigned reads and assembled haplotypes](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.2.1)  3.2.2 [Run more samples](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.2.2) 3.3 [Run HaplotypeCaller on a single bam file in GVCF mode](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.3)  3.3.1 [View resulting GVCF file in the terminal](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.3.1)  3.3.2 [View variants in IGV](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.3.2)  3.3.3 [Run joint genotyping on the CEU Trio GVCFs to generate the final VCF](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.3.3)  3.3.4 [View variants in IGV and compare callsets](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials#3.3.4)

------

### 1 WORKING WITH DATASETS FROM DIFFERENT EXPERIMENTAL DESIGNS

### 1.1 The genome reference: b37

We are using a version of the b37 human genome reference containing only a subset of chromosome 20, which we prepared specially for this tutorial in order to provide a reasonable bundle size for download. It is accompanied by its index and sequence dictionary.

|                            |      |                     |
| -------------------------- | ---- | ------------------- |
| ref/                       |      |                     |
| human_g1k_b37_20.fasta     |      | genome reference    |
| human_g1k_b37_20.fasta.fai |      | fasta index         |
| human_g1k_b37_20.dict      |      | sequence dictionary |

Open up IGV, and load the **Human (1kg, b37+decoy)** reference available on the IGV server (Genomes>Load Genome from Server). We use this reference in IGV because it has a pre-loaded gene track, whereas our custom chromosome-20-only reference does not.

![img](https://us.v-cdn.net/5019796/uploads/FileUpload/50/ea77178d2c4407230d0a0282777951.png)

### 1.2 The test sample: NA12878 Whole-Genome Sequence (WGS)

The biological sample from which the example sequence data was obtained comes from individual NA12878, a member of a 17 sample collection known as CEPH Pedigree 1463, taken from a family in Utah, USA. A trio of two parents and one child from this data set is often referred to as the CEU Trio and is widely used as an evaluation standard (e.g. in the Illumina Platinum Genomes dataset). Note that an alternative trio constituted of the mother (NA12878) and her parents is often also referred to as a CEU Trio. Our trio corresponds to the 2nd generation and one of the 11 grandchildren.

We will begin with a bit of data exploration by looking at the following BAM files derived from NA12878:

1. `NA12878_wgs_20.bam`

   Whole genome sequence (WGS) dataset, paired-end 151 bp reads sequenced on Illumina HiSeqX and fully pre-processed according to the GATK Best Practices for germline DNA.

2. `NA12878_rnaseq_20.bam`

   RNAseq dataset, paired-end 75 bp reads sequenced on Illumina HiSeqX and aligned using STAR 2-pass according to the GATK Best Practices for RNAseq.

3. `NA12878_ICE_20.bam`

   Exome dataset, Illumina Capture Exome (ICE) library, paired-end 76 bp reads sequenced on Illumina HiSeqX, fully pre-processed according to the GATK Best Practices for germline DNA.

4. `NA12878_NEX_20.bam`

   Exome dataset, Illumina Nextera Rapid Capture Exome (NEX) library, paired-end 76 bp reads sequenced on Illumina HiSeqX, fully pre-processed according to the GATK Best Practices for germline DNA.

The sequence data files have been specially prepared as well to match our custom chromosome 20-only reference. They only contain data on chromosome 20, in two pre-determined intervals of interest ranging from positions 20:10,000,000-10,200,000 and 20:15,800,000-16,100,00 to keep file sizes down.

Let’s start by loading the DNA WGS sample of NA12878 (`bams/exp_design/NA12878_wgs_20.bam`), as shown in the screenshots below.

Initially you will not see any data displayed. You need to zoom in to a smaller region for IGV to start displaying reads. You can do that by using the -/+ zoom controls, or by typing in some genome regions coordinates. Here, we’ll zoom into a predetermined interval of interest, so type **20:16,029,744-16,030,079** into the coordinates box. Once you hit the `[Go]` button, you should see something like this:

The top track shows depth of coverage, i.e. the amount of sequence reads present at each position. The mostly grey horizontal bars filling the viewport are the reads. Grey means that those bases match the reference, while colored stripes or base letters (depending on your zoom level) indicate mismatches. You will also see some reads with mapping insertions and deletions, indicated by purple `I` symbols and crossed-out gaps, respectively.

> **TOOL TIP** ![img](https://us.v-cdn.net/5019796/uploads/FileUpload/19/e73d751f1284ba7c8ea6434f80b12c.png)Read details are shown when you hover over them with your mouse--which can be convenient when troubleshooting, but gets annoying quickly. To turn it off, Click the yellow speech bubble in the toolbar and select “Show details on click”.

### 1.3 For comparison: NA12878 Exome Sequence

Next, let’s load our two Exome data sets (File>Load from File), `NA12878_ICE_20.bam` and `NA12878_NEX_20.bam`, and go to position **20:15,873,697-15,875,416**.

You can see from the coverage graph that the ICE sample has more breadth and depth of coverage at this target site, in comparison to the NEX sample. This directly affects our ability to call variants in the leftmost peak, since ICE provides much more depth and NEX has a particularly lopsided distribution of coverage at that site. That’s not to say that ICE is better in general--just that for this target site, in this sequencing run, it provided more even coverage. The overarching point here is that exome kits are not all equivalent and you should evaluate which kit provides the results you need in the regions you care about, before committing to a particular kit for a whole project. As a corollary, comparing exome datasets generated with different kits can be complicated and requires careful evaluation.

### 1.4 Another comparison: NA12878 RNAseq

Lastly, let’s load (File>Load from File) the aligned RNAseq dataset that we have for NA12878 (`NA12878_rnaseq_20.bam`).

 

You’ll notice pale blue lines to the right of center instead of reads. This is because it’s an intronic region! The blue lines connect to reads that are located in the exon. Click on one to see the N operator in the CIGAR string: in the example here, 32M91225N43M indicates that the read covers a 91225 bp intron.

------

### 2 DIAGNOSING UNKNOWN BAMS

### 2.1 View header and check read groups

Now let’s say that you have been brought on to a new project: you will be analyzing sequenced genomes for particular variants in chromosome 20--since you are the chromosome 20 specialist. Your coworker has given you some files that they sequenced a while back. Unfortunately, their lab notebook mostly illegible and lacking in detail where you can read it. So how do you know what’s been done to these files already? Or even if they are good to use still?

Enter Samtools. You can use this tool to open up the bam file your coworker gave you, and check the bam’s record log. Open up your terminal and execute the following:

```
samtools view -H bams/exp_design/NA12878_wgs_20.bam | grep ‘@RG’
```

The bam records log information in the header, so we use `view -H` to ask it to just show us the header. Since we want to see what this sample is, we will also add `| grep ‘@RG’`, which will only grab the line of the header that starts with @RG.

> @RG ID:H0164.2 PL:illumina PU:H0164ALXX140820.2 LB:Solexa-272222 PI:0 DT:2014-08-20T00:00:00-0400 SM:NA12878 CN:BI

You can use the read group information to confirm that this file is what your coworker’s notebook scribbles say it is. You can see that it is indeed the NA12878 sample (SM), and the read group ID H0164.2 (ID) matches, etc. After checking that these identifiers match what you can decipher from your coworker’s writing, call Samtools again. This time we will look at `@PG` to see what tools have been used on this bam file.

```
samtools view -H bams/exp_design/NA12878_wgs_20.bam | grep ‘@PG’
```

Again, this only grabs `@PG` lines from the header, but you will still get a rather long print out in the terminal; we show a single `@PG` entry below.

> @PG ID:bwamem PN:bwamem VN:0.7.7-r441 CL:/seq/software/picard/1.750/3rd_party/bwa_mem/bwa mem -M -t 10 -p /ref/b37.fasta /dev/stdin > /dev/stdout

At the very beginning of each `@PG` entry, there will be a program ID. From this entry, you can see that **BWA MEM**was run on the bam file your coworker gave you--the rest of the entry describes the specific parameters that the tool was run with. Scanning through all the entries, you should see that your coworker ran GATK IndelRealigner, GATK PrintReads, MarkDuplicates, and BWA MEM. These tools correlate with the pre-processing steps that your coworker told you they took: mapping with BWA MEM, duplicate marking with MarkDuplicates, indel realignment with IndelRealigner, and lastly, BQSR with PrintReads*. *How does BQSR correspond to PrintReads? Well, PrintReads is the tool used after running BQSR to apply the recalibration to the bam file itself. Since running BaseRecalibrator didn’t modify the bam file, it isn’t recorded in the bam header, but you can infer that it was run because PrintReads shows up in the header.

### 2.2 Validate the file

Now satisfied that the file your coworker gave you is properly pre-processed from looking at its header, you want to make sure that the body of the bam file wasn’t broken at some point. We will try [diagnosing possible problems](https://www.broadinstitute.org/gatk/blog?id=7567)in the bam using ValidateSamFile.

```
java -jar picard.jar ValidateSamFile \
    I=input.bam \
    MODE=SUMMARY
```

Since we don’t know what kind of errors or warnings we will find, we first run the tool in `SUMMARY` mode. This will output a histogram listing all the errors and warnings in our file.

> \## HISTOGRAM java.lang.String Error Type Count ERROR:MATE_NOT_FOUND 77

That many errors? The file could be badly damaged, but let’s take a closer look. The error here is a MATE_NOT_FOUND, indicating that a read was marked as paired, but that its mate is not found in the file. Now, usually this would be a point of concern, but your coworker told you that this file was subset to a small part of chromosome 20, so it would make sense that some reads mapped within this region and their mates mapped outside the region.

We can safely ignore this warning. For more details on errors and warnings that ValidateSamFile can produce (since you won’t just be running your coworker’s samples forever), check out [this article](https://www.broadinstitute.org/gatk/guide/article?id=7571). For your coworker’s file, though, you are finally ready to move on to…

------

### 3 VARIANT DISCOVERY

### 3.1 Call variants with a position-based caller: UnifiedGenotyper

You found a (typed!) copy of your coworker's variant discovery protocol, so you want to run their bam file following it. It tells you to run the following command:

```
java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -o sandbox/NA12878_wgs_20_UG_calls.vcf \
    -glm BOTH \
    -L 20:10,000,000-10,200,000
```

Reading from the protocol, you see that `-glm BOTH` tells the tool to call both indels and SNPs, while `-L` gives the interval that the bam was subset to--no use wasting time trying to run on the whole genome when you only have data for a small amount.

When the results return, load the original bam file (`bams/exp_design/NA12878_wgs_20.bam`) and the output VCF (`sandbox/NA12878_wgs_20_UG_calls.vcf`) in IGV. Zooming to the coordinates **20:10,002,371-10,002,546**, you will see something like the screenshot below.

 

The variant track shows only variant calls--so at this particular site, there is a homozygous SNP call. (You can click on the variant call for more information on it, too.) The bam track below shows the supporting read data that led to a variant call at that site.

Since this laptop screen is so tiny (our budget went to reagents rather than monitors…) and we can’t zoom out any more vertically, right-click on the bam track and select “Collapsed” view.

This gives us a better overview of what the data looks like in this region: good even coverage, not too much noise in the region, and reasonable allele balance (mostly variant supports the homozygous variant call). Based on the information we see here, this should be a clear variant site.

### 3.2 Call variants with HaplotypeCaller

While preparing for this project, though, you recall hearing about another variant caller: HaplotypeCaller. And, looking on GATK’s website, you see that it recommends calling your variants using HaplotypeCaller over the old UnifiedGenotyper. The new algorithm calls both SNP and indel variants simultaneously via local de-novo assembly of haplotypes in an active region. Essentially, when this variant caller finds a region with signs of variation, it tosses out the old alignment information (from BWA MEM) and performs a local realignment of reads in that region. This makes HaplotypeCaller more accurate in regions that are traditionally difficult to call--such as areas that contain different types of variants close together. Position-based callers like UnifiedGenotyper simply can’t compete.

You decide to re-run your sample with the new variant caller to see if it makes a difference. Tool documentation on the website gives you a basic command to run, and you add your coworker’s interval trick (`-L`) in as well.

```
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -o sandbox/NA12878_wgs_20_HC_calls.vcf \
    -L 20:10,000,000-10,200,000
```

Load the output VCF (`sandbox/NA12878_wgs_20_HC_calls.vcf`) in IGV to compare the HC calls to the previously-loaded UG calls.

We see that HC called the same C/T SNP as UG, but it also called another variant, a homozygous variant insertion of three T bases. How is this possible when so few reads seem to support an insertion at this position?

> **TOOL TIP** ![img](https://us.v-cdn.net/5019796/uploads/FileUpload/29/654cffb26222c4bbacf32e415fed54.png)When you encounter indel-related weirdness, turn on the display of soft-clips, which IGV turns off by default. Go to View > Preferences > Alignments and select “Show soft-clipped bases”

With soft clip display turned on, the region lights up with variants. This tells us that the aligner (here, BWA MEM) had a lot of trouble mapping reads in the region. It suggests that HaplotypeCaller may have found a different alignment after performing its local graph assembly step. This reassembled region provided HaplotypeCaller with enough support to call the indel that UnifiedGenotyper missed.

### *3.2.1 View realigned reads and assembled haplotypes*

But we’re not satisfied with “probably” here. Let’s take a peek under the hood of HaplotypeCaller. You find that HaplotypeCaller has a parameter called `-bamout`, which allows you to ask for the realigned version of the bam. That realigned version is what HaplotypeCaller uses to make its variant calls, so you will be able to see if a realignment fixed the messy region in the original bam.

You decide to run the following command:

```
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -o sandbox/NA12878_wgs_20_HC_calls_debug.vcf \
    -bamout sandbox/NA12878_wgs_20.HC_out.bam \
    -forceActive -disableOptimizations \
    -L 20:10,002,371-10,002,546 -ip 100
```

Since you are only interested in looking at that messy region, you decide to give the tool a narrowed interval with `-L 20:10,002,371-10,002,546`, with a 100 bp padding on either side using `-ip 100`. To make sure the tool does perform the reassembly in that region, you add in the `-forceActive` and `-disableOptimizations`arguments.

Load the output BAM (`sandbox/NA12878_wgs_20.HC_out.bam`) in IGV, and switch to Collapsed view once again. You should still be zoomed in on coordinates **20:10,002,371-10,002,546**, and have the original bam track loaded for comparison.

After realignment by HaplotypeCaller (the bottom track), almost all the reads show the insertion, and the messy soft clips from the original bam are gone. Expand the reads in the output BAM (right click>Expanded view), and you can see that all the insertions are in phase with the C/T SNP.

There is more to a BAM than meets the eye--or at least, what you can see in this view of IGV. Right-click on the reads to bring up the view options menu. Select **Color alignments by**, and choose **read group**. Your gray reads should now be colored similar to the screenshot below.

Some of the first reads, shown in red at the top of the pile, are not real reads. These represent artificial haplotypes that were constructed by HaplotypeCaller, and are tagged with a special read group identifier, “ArtificialHaplotype,” so they can be visualized in IGV. You can click on an artificial read to see this tag under **RG**.

We see that HaplotypeCaller considered six possible haplotypes, because there is more than one variant in the same ActiveRegion. Zoom out further , and we can see that two ActiveRegions were examined within the scope of the interval we provided (with padding).

![img](https://us.v-cdn.net/5019796/uploads/FileUpload/ee/19be1d025100c211f410161ea1917f.png)

### *3.2.2 Run more samples*

You’ve decided that perhaps HaplotypeCaller will work better for your project. However, since you have been working on this protocol update, your coworker found two more samples--they were in a different folder on their computer for reasons you can’t figure out. Regardless, you now need to joint call all the samples together. So, using the same command as before, you’ve tacked on the two additional bam files.

```
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -I bams/trio-calling/NA12877_wgs_20.bam \
    -I bams/trio-calling/NA12882_wgs_20.bam \
    -o sandbox/NA12878_wgs_20_HC_jointcalls.vcf \
    -L 20:10,000,000-10,200,000
```

You notice after entering that, that HaplotypeCaller takes a much longer time to return than other tasks we have run so far. You decide to check the results of this command later, and do some digging on how to make things go faster.

### 3.3 Run HaplotypeCaller on a single bam file in GVCF mode

Every time your coworker finds a new folder of samples, you’ll have to re-run all the samples using this increasingly slower HaplotypeCaller command. You’ve also been approved for a grant and intend to send your own samples out for sequencing, so there are those to add in as well. You could just wait until you have all the samples gathered, but that could be a while and your PI wants to see some preliminary results soon. You read about a new GATK workflow that lets you make everyone happy: the GVCF workflow.

The first step in variant discovery is to run HaplotypeCaller in GVCF mode on each individual bam file. This is basically running HaplotypeCaller as you did before, but with `-ERC GVCF` added to the command. You first want to run HaplotypeCaller in GVCF mode on the NA12878 bam. (In the interest of time, we have supplied the other sample GVCFs in the bundle, but normally you would run them individually in the same way as the first.) This will produce a GVCF file that contains genotype likelihoods for each variant position as well as blocks for each interval where no variant is likely. You’ll see what this looks like more in a minute.

```
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -o sandbox/NA12878_wgs_20.g.vcf \
    -ERC GVCF \
    -L 20:10,000,000-10,200,000
```

### *3.3.1 View resulting GVCF file in the terminal*

Since a GVCF is a new file type for your workflow, let’s take a look at the actual content first. You can do this in the terminal by typing this command:

more sandbox/NA12878_wgs_20.g.vcf

As you scroll through the file (hit `[ENTER]` to scroll, `[CTRL]`+`[C]` to exit), note the NON_REF allele defined in the header.

> \##ALT=<ID=NON_REF,Description=”Represents any possible alternative allele at this location”>

Also note the GVCF blocks defined later in the header. The reference (non-variant) blocks are recorded in the GVCF file, in blocks separated by genotype quality.

> \##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive) ##GVCFBlock1-2=minGQ=1(inclusive),maxGQ=2(exclusive) ##GVCFBlock10-11=minGQ=10(inclusive),maxGQ=11(exclusive) ##GVCFBlock11-12=minGQ=11(inclusive),maxGQ=12(exclusive)

Finally, while scrolling through the records, we can see the reference blocks and **variant sites**.

> 20 10000115 . G . . END=10000116 GT:DP:GQ:MIN_DP:PL 0/0:25:69:25:0,69,1035 **20 10000117 . C T, 262.77 . BaseQRankSum=-0.831;ClippingRankSum=-0.092;DP=23;MLEAC=1,0;MLEAF=0.500,0.00;MQ=60.47;MQRankSum=1.446;ReadPosRankSum=0.462 GT:AD:DP:GQ:PL:SB 0/1:11,12,0:23:99:291,0,292,324,327,652:9,2,9,3** 20 10000118 . T . . END=10000123 GT:DP:GQ:MIN_DP:PL 0/0:25:63:24:0,63,945

Every site in the interval we analyzed is represented here--whether it be by a variant call, a reference call, or a reference block. This helps to distinguish between a “no call” (we don’t have enough data to make a call) and a “reference call” (we have evidence that the site matches the reference).

### *3.3.2 View variants in IGV*

Now, text in a terminal window can be rather hard to read, so let’s take a look at the GVCFs in IGV. Start a new session to clear your IGV screen, then load the three GVCFs (`sandbox/NA12878_wgs_20.g.vcf`, `gvcfs/NA12877_wgs_20.g.vcf`, `gvcfs/NA12882_wgs_20.g.vcf`). You should already be zoomed in on **20:10,002,371-10,002,546** from our previous section, and see this:

Notice anything different from the VCF? Along with the colorful variant sites, you see many gray blocks in the GVCF representing the non-variant intervals. Most of the gray blocks are next to each other, but are not grouped together, because they belong to different GQ blocks. The chief difference between the GVCF here and the next step’s VCF is the lack of reference blocks (the gray bits). Only very low-confidence variant sites will be removed in the VCF, based on the QUAL score.

### *3.3.3 Run joint genotyping on the CEU Trio GVCFs to generate the final VCF*

The last step is to joint call all your GVCF files using the GATK tool GenotypeGVCFs. After looking in the tool documentation, you run this command:

```
java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs \
    -R ref/human_g1k_b37_20.fasta \
    -V sandbox/NA12878_wgs_20.g.vcf \
    -V gvcfs/NA12877_wgs_20.g.vcf \
    -V gvcfs/NA12882_wgs_20.g.vcf \
    -o sandbox/CEUTrio_wgs_20_GGVCFs_jointcalls.vcf \
    -L 20:10,000,000-10,200,000
```

That returned much faster than the HaplotypeCaller step--and a good thing, too, since this step is the one you’ll need to re-run every time your coworker finds a “new” sample buried in their messy file structure. But does calling this way really give you good results? Let’s take a look.

### *3.3.4 View variants in IGV and compare callsets*

Load the joint called VCF from normal HaplotypeCaller, section 3.2.1 (`sandbox/NA12878_wgs_20_HC_jointcalls.vcf`), and GenotypeGVCFs, section 3.3.3 (`sandbox/CEUTrio_wgs_20_GGVCFs_jointcalls.vcf`). Change your view to look at **20:10,002,584-10,002,665**, and you will see:

At this site, the father NA12877 is heterozygous for a G/T SNP, and the mother, NA12878, and son, NA12882, are homozygous variant for the same SNP. These calls match up, and you figure that the calls between GenotypeGVCFs and HaplotypeCaller, when run in multisample mode, are essentially equivalent. (And if you did some digging, you would find some marginal differences in borderline calls.) However, the GVCF workflow allows you to be more flexible. Every time your PI wants an update on the project, you can simply re-run the quick GenotypeGVCFs step on all the samples you have gathered so far. The expensive and time-consuming part of calculating genotype likelihoods only needs to be done once on each sample, so you won’t have to spend all your grant money on compute to rerun the whole cohort every time you have a new sample.

You have successfully run your coworker’s samples, and you’ve found that the most effective workflow for you is the most recent GVCF workflow. Your next step takes you to filtering the callset with either VQSR or hard filters--but you decide to take a break before tackling the next part of the workflow.