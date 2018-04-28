#! /usr/bin/env python
"""Module to do the Dstatistics, with the error correction. Works on a joint
vcf file as input, and the sample configurations in a separate file, one line
per configuration. One important note is that the vcf file used must be an
all sites - all sample vcf, not a polymorphisms only or a single sample vcf.
The statistics computed in this paper are based on the math presented in
the Neanderthal paper, Green et al.

There is an accompanying module, TTStatistic which computes the 2x2
statistic for computing split times.
"""
__author__ = "Shyam Gopalakrishnan <shyam@snm.ku.dk>"
__date__ = "15th March 2018"

__version__ = "$Revision: 01 $"

__credits__ = """Marc de Manuel Montero, Jazmin Madrigal Ramos for discussions
and motivating me to finally write this code."""

import sys
import vcf
import argparse
import numpy as np
import gzip as gz

baseToAllele = {'A': 0,'C': 1,'G': 2,'T': 3}

def countGenotypeCombinationsSampled(vcfrecord, sampleConfigs):
    """Counts genotypes from list of genoytpes.

    Given genotypes from a line in the vcf, and an array of possible sample configs,
    count the number of genotype combinations, viz., AAAA, AAAB,...ABCD, by sampling
    a single base for each sample.

    Args:
        vcfrecord: a single record from a vcf.
        sampleConfig: Array with each element being a list of 4 indices with
            samples (H1, H2, H3, O).

    Returns:
        Array of list of counts for each genotype class.

    Raises:
        LookupError: Sample in sample config file is not present in vcf.
    """

    curGenoCount = np.zeros((15, len(sampleConfigs)))
    alleles = [vcfrecord.REF]
    alleles.extend(vcfrecord.ALT)
    indices = [ baseToAllele[x] for x in alleles ]
    for (cnt, (H1samps, H2samps, H3samps, Osamps)) in enumerate(sampleConfigs):
        H1Alleles = np.zeros(4)
        H2Alleles = np.zeros(4)
        H3Alleles = np.zeros(4)
        OAlleles = np.zeros(4)
        for samp in Osamps:
            depths = np.array(vcfrecord.genotype(samp)['AD'])
            depths /= np.sum(depths)
            index = np.random.choice(indices, size=1, p=depths)
            OAlleles[index] += 1
        for samp in H1samps:
            depths = np.array(vcfrecord.genotype(samp)['AD'])
            depths /= np.sum(depths)
            index = np.random.choice(indices, size=1, p=depths)
            H1Alleles[index] += 1
        for samp in H2samps:
            depths = np.array(vcfrecord.genotype(samp)['AD'])
            depths /= np.sum(depths)
            index = np.random.choice(indices, size=1, p=depths)
            H2Alleles[index] += 1
        for samp in H3samps:
            depths = np.array(vcfrecord.genotype(samp)['AD'])
            depths /= np.sum(depths)
            index = np.random.choice(indices, size=1, p=depths)
            H3Alleles[index] += 1
        curGenoCount[:,cnt] = collapseAlleleCounts(H1Alleles, H2Alleles, H3Alleles, OAlleles, indices)


def collapseAlleleCounts(H1Alleles, H2Alleles, H3Alleles, OAlleles, indices):
    """Convert the allele counts in each set to fractional counts for this siteself.

    Convert allele counts to fractional count for each of the 15 classes. In the
    case of a single sample with randomly sampled allele, this is going to be
    only 1 class with a count of 1 and all others zero. Order of classes is the
    same as found in supplementary material of Green et al. 2010 paper A draft
    sequence of the Neanderthal genome.

    Args:
        H1Alleles: Counts of the 4 alleles in H1 samples.
        H2Alleles: Counts of the 4 alleles in H2 samples.
        H3Alleles: Counts of the 4 alleles in H3 samples.
        OAlleles: Counts of the 4 alleles in O samples.
        indices: The indices of the bases of ref and alt alleles.

    Returns:
        A np 1-D array, summing to 1, with fractions of sites in each class,
        going from AAAA, AABA, ... to BCDA.

    Raises:
        Nothing.
    """
    curgenocnt = np.zeros(15)
    ## at least 1 tree branch has no info.
    if np.sum(OAlleles) == 0 or np.sum(H1Alleles) = 0 or
       np.sum(H2Alleles) == 0 or np.sum(H3Alleles) == 0:
       return curgenocnt
    ## outgroup is heterozygous :(
    if np.sum(OAlleles == 0) != 3:
        return curgenocnt
    Aindex = np.argwhere(OAlleles > 0)[0,0]
    if len(indices) == 1:
        #all AAAA case so easy to fix,
        curgenocnt[0] += 1
        return curgenocnt
    Bindex = 
    if len(indices) == 2:
        ## fix outgroup allele to A

def readSampleConfigFile(configFilename):
    """Read the sample configuration file, and return 4-tuples of sample names.

    Read a sample configuration file, with each line representing a sample config,
    and convert this to a tuple of sample names (H1,H2,H3,O). The input file is a
    4 column file, with each column containing a comma separated list of samples
    that belong to H1, H2, H3 and O.

    Args:
        configFilename: Name of configuration file, with 4 names per line, and
            each line representing a configuration (H1, H2, H3, O)

    Returns:
        Array of list of sample names for each configuration in the file.

    Raises:
        IOError: When the format of the input config file name is incorrect.
    """
    configFile = open(configFilename)
    sampleConfigs = []
    for line in configFile:
        toks = line.strip().split()
        if len(toks) != 4:
            raise IOError("Incorrect number of columns in config file.")
        sampleConfigs.append((toks[0].split(","), toks[1].split(","),
                              toks[2].split(","), toks[3].split(",")))
        if len(toks[0]) > 1 or len(toks[1]) > 1 or len(toks[2]) > 1 or len(toks[3]) > 1:
            raise IOError("Currently only implements one sample per tree branch.")
    configFile.close()
    return(sampleConfigs)

def computeDStats(genoCounts, blockSize):
    """Compute Dstats and error give genoypte counts and block size.

    Compute the D statistics from the ABBA and BABA counts, given the genotype
    counts for each configuration. Also compute the D statistic obtained from the
    error alone. Further, compute the variance of the statistic using the block
    jackknife, using the block size specified by user.

    Args:
        genoCount: The dictionary containing the list of genoypte counts for each samples
            configuration, with key as chrom and posself.
        blockSize: The size of the block used for block jackknife.

    Returns:
        Array of D-statistics and the variances from block jackknife

    Raises:
        ValueError: When computations fail due to some numerical error.
    """
    pass

def main(args):
    """The boss function. Calls all the subfunctions to get the job done. Basically
    pointy haired boss from Dilbert. At least it is not a Wally.

    Args:
        args: arguments from the argparser.

    Returns:
        Nothing, what did you expect from the boss.

    Raises:
        IOError: When file is not found or raised by subfunctions.
        LookupError: Raised by subfunctions
    """
    if args.vcf[-3:] == ".gz":
        vcffile = gz.open(args.vcf, 'r')
    else:
        vcffile = open(args.vcf, 'r')
    sampleConfigs = readSampleConfigFile(args.sampleConfigs)
    genoCounts = np.zeros()

    ## Dictionary to store the array counts.
    ## key is tuple (chr, pos) and value is array
    ## of array of genotype counts - 1 per sample
    ## configuration.
    genoCounts = {}
    reader = vcf.Reader(vcffile)
    for record in reader:
        if record.QUAL
        genoCounts[(record.CHROM, int(record.POS))] = countGenotypeCombinations(record, sampleConfigs)


if __name__ == "__main__":
    ## Write the argument parser.
    parser = argparse.ArgumentParser(
        description="Compute D-statistics using 4 samples in the list of samples in the vcf file.\n")

    parser.add_argument("-v", "--vcf",
        help="path of vcf file to be processed.",
        required=True,
        type=str)

    parser.add_argument("-s", "--sampleConfigs",
        help="file with sample configurations (H1, H2, H3, O) - one per line.",
        required=True,
        type=str)

    parser.add_argument("-o", "--output",
        help="output file name.",
        required=False,
        default="",
        type=str)

    parser.add_argument("-d", "--min-depth",
        help="minimum per individual depth at site.",
        required=False,
        default=0,
        type=int)

    parser.add_argument("-D", "--max-depth",
        help="maximum per individual depth at site.",
        required=False,
        default=10000,
        type=int)

    parser.add_argument("-g", "--min-genoqual",
        help="minimum genotype quality.",
        required=False,
        default=0,
        type=int)

    parser.add_argument("-m", "--mode",
        help="1: sample one read (default)\n2: use genotype",
        required=False,
        default=1,
        type=int)

    args = parser.parse_args()

    main(args)
