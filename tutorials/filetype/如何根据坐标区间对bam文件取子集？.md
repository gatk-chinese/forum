# 如何根据坐标区间对bam文件取子集？

原文链接：[ (How to) Create a snippet of reads corresponding to a genomic interval](https://software.broadinstitute.org/gatk/documentation/article?id=6517)

 

#### Tools involved

- [PrintReads](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php)

#### Prerequisites

- Installed GATK tools
- Reference genome
- Coordinate-sorted, aligned and indexed BAM

#### Download example data

- Use the [advanced tutorial bundle](http://gatkforums.broadinstitute.org/discussion/4610/)'s human_g1k_v37_decoy.fasta as reference
- [tutorial_6517.tar.gz](https://drive.google.com/open?id=0BzI1CyccGsZiTmlDLW13MXdTSG8) contains four files: 6517_2Mbp_input.bam and .bai covering reads aligning to 10:90,000,000-92,000,000 and 6517_1Mbp_output.bam and .bai covering 10:91,000,000-92,000,000

#### Related resources

- This *How to* is referenced in a tutorial on [(How to) Generate an unmapped BAM (uBAM)](http://gatkforums.broadinstitute.org/discussion/6484/).
- See [this tutorial](http://gatkforums.broadinstitute.org/discussion/2909/) to coordinate-sort and index a BAM.
- See [here](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--unsafe) for command line parameters accepted by all GATK tools.
- For applying interval lists, e.g. to whole exome data, see [When should I use L to pass in a list of intervals](http://gatkforums.broadinstitute.org/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals).

------

### Create a snippet of reads corresponding to a genomic interval using PrintReads

PrintReads merges or subsets sequence data. The tool automatically applies MalformedReadFilter and BadCigarFilter to filter out certain types of reads that cause problems for downstream GATK tools, e.g. reads with mismatching numbers of bases and base qualities or reads with CIGAR strings containing the N operator.

- To create a test snippet of RNA-Seq data that retains reads with Ns in CIGAR strings, use `-U ALLOW_N_CIGAR_READS`.

Subsetting reads corresponding to a genomic interval using PrintReads requires reads that are aligned to a reference genome, coordinate-sorted and indexed. Place the `.bai` index in the same directory as the `.bam` file.

```
java -Xmx8G -jar /path/GenomeAnalysisTK.jar \
    -T PrintReads \ 
    -R /path/human_g1k_v37_decoy.fasta \ #reference fasta
    -L 10:91000000-92000000 \ #desired genomic interval chr:start-end
    -I 6517_2Mbp_input.bam \ #input
    -o 6517_1Mbp_output.bam 
```

This creates a subset of reads from the input file, `6517_2Mbp_input.bam`, that align to the interval defined by the `-L` option, here a 1 Mbp region on chromosome 10. The tool creates two new files, `6517_1Mbp_output.bam`and corresponding index `6517_1Mbp_output.bai`.

- For paired reads, the tool does not account for reads whose mate aligns outside of the defined interval. To filter these lost mate reads, use RevertSam's `SANITIZE` option.

To process large files, also designate a temporary directory.

```
    TMP_DIR=/path/shlee #sets environmental variable for temporary directory
```