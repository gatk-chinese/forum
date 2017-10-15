# 如何安装GATK最佳实践的必备软件

原文链接：[ (howto) Install all software packages required to follow the GATK Best Practices](https://software.broadinstitute.org/gatk/documentation/article?id=2899)

这个软件安装有点多，跟GATK最佳实践大部分内容是重合的

#### Objective

Install all software packages required to follow the GATK Best Practices.

#### Prerequisites

To follow these instructions, you will need to have a basic understanding of the meaning of the following words and command-line operations. If you are unfamiliar with any of the following, you should consult a more experienced colleague or your systems administrator if you have one. There are also many good online tutorials you can use to learn the necessary notions.

- Basic Unix environment commands
- Binary / Executable
- Compiling a binary
- Adding a binary to your path
- Command-line shell, terminal or console
- Software library

You will also need to have access to an ANSI compliant C++ compiler and the tools needed for normal compilations (make, shell, the standard library, tar, gunzip). These tools are usually pre-installed on Linux/Unix systems. **On MacOS X, you may need to install the MacOS Xcode tools.** See <https://developer.apple.com/xcode/> for relevant information and software downloads. The XCode tools are free but an AppleID may be required to download them.

Starting with version 3.6, the GATK requires Java Runtime Environment version 1.8 (Java 8). Previous versions down to 2.6 required JRE 1.7, and earlier versions required 1.6. All Linux/Unix and MacOS X systems should have a JRE pre-installed, but the version may vary. To test your Java version, run the following command in the shell:

```
java -version 
```

This should return a message along the lines of ”java version 1.8.0_25” as well as some details on the Runtime Environment (JRE) and Virtual Machine (VM). If you have a version that does not match the requirements stated above for the version of GATK you are running, the GATK may not run correctly or at all. The simplest solution is to install an additional JRE and specify which you want to use at the command-line. To find out how to do so, you should seek help from your systems administrator.

#### Software packages

1. BWA
2. SAMtools
3. Picard
4. Genome Analysis Toolkit (GATK)
5. IGV
6. RStudio IDE and R libraries ggplot2 and gsalib

*Note that the version numbers of packages you download may be different than shown in the instructions below. If so, please adapt the number accordingly in the commands.*

------

### 1. BWA

Read the overview of the BWA software on the [BWA project homepage](http://bio-bwa.sourceforge.net/), then download the [latest version of the software package](http://sourceforge.net/projects/bio-bwa/files/).

- Installation

Unpack the tar file using:

```
tar xvzf bwa-0.7.12.tar.bz2 
```

This will produce a directory called `bwa-0.7.12` containing the files necessary to compile the BWA binary. Move to this directory and compile using:

```
cd bwa-0.7.12
make
```

The compiled binary is called `bwa`. You should find it within the same folder (`bwa-0.7.12` in this example). You may also find other compiled binaries; at time of writing, a second binary called `bwamem-lite` is also included. You can disregard this file for now. Finally, just add the BWA binary to your path to make it available on the command line. This completes the installation process.

- Testing

Open a shell and run:

```
bwa 
```

This should print out some version and author information as well as a list of commands. As the **Usage** line states, to use BWA you will always build your command lines like this:

```
bwa <command> [options] 
```

This means you first make the call to the binary (`bwa`), then you specify which command (method) you wish to use (e.g. `index`) then any options (*i.e.* arguments such as input files or parameters) used by the program to perform that command.

------

### 2. SAMtools

Read the overview of the SAMtools software on the [SAMtools project homepage](http://samtools.sourceforge.net/), then download the [latest version of the software package](http://sourceforge.net/projects/samtools/files/).

- Installation

Unpack the tar file using:

```
tar xvjf samtools-0.1.2.tar.bz2 
```

This will produce a directory called `samtools-0.1.2` containing the files necessary to compile the SAMtools binary. Move to this directory and compile using:

```
cd samtools-0.1.2 
make 
```

The compiled binary is called `samtools`. You should find it within the same folder (`samtools-0.1.2` in this example). Finally, add the SAMtools binary to your path to make it available on the command line. This completes the installation process.

- Testing

Open a shell and run:

```
samtools 
```

This should print out some version information as well as a list of commands. As the **Usage** line states, to use SAMtools you will always build your command lines like this:

```
samtools <command> [options] 
```

This means you first make the call to the binary (`samtools`), then you specify which command (method) you wish to use (e.g. `index`) then any options (*i.e.* arguments such as input files or parameters) used by the program to perform that command. This is a similar convention as used by BWA.

------

### 3. Picard

Read the overview of the Picard software on the [Picard project homepage](http://broadinstitute.github.io/picard/), then download the [latest version](https://github.com/broadinstitute/picard/releases/)(currently 2.4.1) of the package containing the pre-compiled program file (the picard-tools-2.x.y.zip file).

- Installation

Unpack the zip file using:

```
tar xjf picard-tools-2.4.1.zip 
```

This will produce a directory called `picard-tools-2.4.1` containing the Picard jar files. Picard tools are distributed as a pre-compiled Java executable (jar file) so there is no need to compile them.

Note that it is not possible to add jar files to your path to make the tools available on the command line; you have to specify the full path to the jar file in your java command, which would look like this:

```
java -jar ~/my_tools/jars/picard.jar <Toolname> [options]
```

*This syntax will be explained in a little more detail further below.*

However, you can set up a shortcut called an "environment variable" in your shell profile configuration to make this easier. The idea is that you create a variable that tells your system where to find a given jar, like this:

```
PICARD = "~/my_tools/jars/picard.jar"
```

So then when you want to run a Picard tool, you just need to call the jar by its shortcut, like this:

```
java -jar $PICARD <Toolname> [options]
```

The exact way to set this up depends on what shell you're using and how your environment is configured. We like [this overview and tutorial](https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-a-linux-vps) which explains how it all works; but if you are new to the command line environment and you find this too much too deal with, we recommend asking for help from your institution's IT support group.

This completes the installation process.

- Testing

Open a shell and run:

```
java -jar picard.jar -h 
```

This should print out some version and usage information about the `AddOrReplaceReadGroups.jar` tool. At this point you will have noticed an important difference between BWA and Picard tools. To use BWA, we called on the BWA program and specified which of its internal tools we wanted to apply. To use Picard, we called on Java itself as the main program, then specified which jar file to use, knowing that one jar file = one tool. This applies to all Picard tools; to use them you will always build your command lines like this:

```
java -jar picard.jar <ToolName> [options] 
```

This means you first make the call to Java itself as the main program, then specify the `picard.jar` file, then specify which tool you want, and finally you pass whatever other arguments (input files, parameters etc.) are needed for the analysis.

Note that the command-line syntax of Picard tools has recently changed from `java -jar <ToolName>.jar` to `java -jar picard.jar <ToolName>`. We are using the newer syntax in this document, but some of our other documents may not have been updated yet. If you encounter any documents using the old syntax, let us know and we'll update them accordingly. If you are already using an older version of Picard, either adapt the commands or better, upgrade your version!

Next we will see that GATK tools are called in essentially the same way, although the way the options are specified is a little different. The reasons for how tools in a given software package are organized and invoked are largely due to the preferences of the software developers. They generally do not reflect strict technical requirements, although they can have an effect on speed and efficiency.

### 4. Genome Analysis Toolkit (GATK)

Hopefully if you're reading this, you're already acquainted with the [purpose of the GATK](http://www.broadinstitute.org/gatk/about), so go ahead and download the [latest version of the software package](http://www.broadinstitute.org/gatk/download).

In order to access the downloads, you need to register for a free account on the [GATK support forum](http://gatkforums.broadinstitute.org/). You will also need to read and accept the license agreement before downloading the GATK software package. Note that if you intend to use the GATK for commercial purposes, you will need to purchase a license. See the [licensing page](https://www.broadinstitute.org/gatk/about/#licensing) for an overview of the commercial licensing conditions.

- Installation

Unpack the tar file using:

```
tar xjf GenomeAnalysisTK-3.3-0.tar.bz2 
```

This will produce a directory called `GenomeAnalysisTK-3.3-0` containing the GATK jar file, which is called `GenomeAnalysisTK.jar`, as well as a directory of example files called `resources`. GATK tools are distributed as a single pre-compiled Java executable so there is no need to compile them. Just like we discussed for Picard, it's not possible to add the GATK to your path, but you can set up a shortcut to the jar file using environment variables as described above.

This completes the installation process.

- Testing

Open a shell and run:

```
java -jar GenomeAnalysisTK.jar -h 
```

This should print out some version and usage information, as well as a list of the tools included in the GATK. As the **Usage** line states, to use GATK you will always build your command lines like this:

```
java -jar GenomeAnalysisTK.jar -T <ToolName> [arguments] 
```

This means that just like for Picard, you first make the call to Java itself as the main program, then specify the `GenomeAnalysisTK.jar` file, then specify which tool you want, and finally you pass whatever other arguments (input files, parameters etc.) are needed for the analysis.

------

### 5. IGV

The Integrated Genomics Viewer is a genome browser that allows you to view BAM, VCF and other genomic file information in context. It has a graphical user interface that is very easy to use, and can be downloaded for free (though registration is required) from [this website](https://www.broadinstitute.org/igv/home). We encourage you to read through IGV's very helpful [user guide](https://www.broadinstitute.org/software/igv/UserGuide), which includes many detailed tutorials that will help you use the program most effectively.

------

### 6. RStudio IDE and R libraries ggplot2 and gsalib

Download the [latest version of RStudio IDE](http://www.rstudio.com/). The webpage should automatically detect what platform you are running on and recommend the version most suitable for your system.

- Installation

Follow the installation instructions provided. Binaries are provided for all major platforms; typically they just need to be placed in your Applications (or Programs) directory. Open RStudio and type the following command in the console window:

```
install.packages("ggplot2") 
```

This will download and install the ggplot2 library as well as any other library packages that ggplot2 depends on for its operation. Note that some users have reported having to install two additional package themselves, called `reshape` and `gplots`, which you can do as follows:

```
install.packages("reshape")
install.packages("gplots")
```

Finally, do the same thing to install the gsalib library:

```
install.packages("gsalib")
```

This will download and install the gsalib library.

**Important note**

If you are using a recent version of `ggplot2` and a version of GATK older than 3.2, you may encounter an error when trying to generate the BQSR or VQSR recalibration plots. This is because until recently our scripts were still using an older version of certain `ggplot2` functions. This has been fixed in GATK 3.2, so you should either upgrade your version of GATK (recommended) or downgrade your version of ggplot2. If you experience further issues generating the BQSR recalibration plots, please see [this tutorial](http://www.broadinstitute.org/gatk/guide/article?id=4294).