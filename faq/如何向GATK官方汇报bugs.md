# 如何向GATK官方汇报bugs?

原文链接：[ What is "Phone Home" and how does it affect me? ](https://software.broadinstitute.org/gatk/documentation/article?id=1250)

In GATK versions produced between September 2010 and May 2016, the GATK had a "Phone Home" usage reporting feature that sent us information about each GATK run via the Broad filesystem (within the Broad) and Amazon's S3 cloud storage service (outside the Broad). This feature was enabled by default and required a key to be disabled (for running offline or for regulatory reasons).

**The Phone Home feature was removed in version 3.6.** Keys are no longer necessary, so if you had one, you can stop using it. We do not expect that including Phone Home arguments in GATK command lines would cause any errors (so this should not break any scripts), but let us know if you run into any trouble.

Note that keys remain necessary for disabling Phone Home in older versions of GATK. See further below for details on how to obtain a key.

------

### How Phone Home helped development

At the time, the information provided by the Phone Home feature was critical in driving improvements to the GATK:

- By recording detailed information about each error that occurs, it enabled GATK developers to **identify and fix previously-unknown bugs**in the GATK.
- It allowed us to better understand how the GATK is used in practice and **adjust our documentation and development goals** for common use cases.
- It gave us a picture of **which versions** of the GATK are in use over time, and how successful we've been at encouraging users to migrate from obsolete or broken versions of the GATK to newer, improved versions.
- It told us **which tools** were most commonly used, allowing us to monitor the adoption of newly-released tools and abandonment of outdated tools.
- It provided us with a sense of the **overall size of our user base** and the major organizations/institutions using the GATK.

------

### What information was sent to us

Below are two example GATK Run Reports showing exactly what information is sent to us each time the GATK phones home.