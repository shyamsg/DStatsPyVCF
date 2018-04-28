#! /usr/bin/env python
"""Module to do the Two-by-Two statistic from vcf files. Works on a joint
vcf file as input, and pair of samples in a separate file, one line per pair.
One important note is that the vcf file used must be an all sites - all sample
vcf, not a polymorphisms only or a single sample vcf. The statistics computed
in this script are computed using math presented in the African paper, and in
Skoglund et al. The vcf file has to be sorted, will not be checked but ... if it
is not.. then trouble!!!
"""
__author__ = "Shyam Gopalakrishnan <shyam@snm.ku.dk>"
__date__ = "15th March 2018"

__version__ = "$Revision: 01 $"

__credits__ = """Marc de Manuel Montero, Mikkel Sinding and Jazmin Madrigal
Ramos for discussions and motivating me to finally write this code. And python!"""

import sys
import vcf
import argparse
import numpy as np
import gzip as gz

def countGenotypeCombinationsSampled(vcfrecord, sampleConfigs, mind, maxd):
    """Counts genotypes from list of genoytpes.

    Given a vcf record, count the genoypte combinations for each of the pairs of
    samples, using a given ancestral sample. The combination of genotypes is
    indexed like as follows:
    0: both ancestral homozygotes
    1: samp1 het, samp2 ancestral homozygote
    2: samp1 ancestral homozygote, samp2 het
    3: samp1 alternate homozygote, samp1 ancestral homozygote
    4: samp1 ancestral homozygote, samp1 alternate homozygote
    5: both heterozygote
    6: samp1 alternate homozygote, samp2 heterozygote
    7: samp1 heterozygote, samp2 alternate homozygote
    8: both alternate homozygote

    Args:
        vcfrecord: a single record from a vcf.
        sampleConfig: Array with each element being a list of 3 samples
        names (S1, S2, O).

    Returns:
        return array of list of indices, from list above.

    Raises:
        LookupError: Sample in sample config file is not present in vcf.
    """
    genoClass = []
    for (s1, s2, o) in sampleConfigs:
        if o == "REF":
            oallele = vcfrecord.REF
        else:
            ocall = vcfrecord.genotype(o)
            odepth = ocall['DP']
            if odepth < mind or odepth > maxd:
                genoClass.append("-1")
                continue
            ogeno = ocall['GT']
            if ogeno == "./." or ogeno == "0/1":
                genoClass.append("-1")
                continue
            else:
                oallele = ogeno.split("/")[0]
        s1call = vcfrecord.genotype(s1)
        s1depth = s1call['DP']
        s2call = vcfrecord.genotype(s2)
        s2depth = s2call['DP']
        if s1depth < mind or s1depth > maxd or s2depth < mind or s2depth > maxd:
            genoClass.append('-1')
            continue
        else:
            if s1call['GT'] == "./." or s2call['GT'] == "./.":
                genoClass.append("-1")
            s1cnt = s1call['GT'].count(oallele)
            s2cnt = s2call['GT'].count(oallele)
            if s1cnt == 0 and s2cnt == 0: genoClass.append(0)
            elif s1cnt == 1 and s2cnt == 0: genoClass.append(1)
            elif s1cnt == 0 and s2cnt == 1: genoClass.append(2)
            elif s1cnt == 2 and s2cnt == 0: genoClass.append(3)
            elif s1cnt == 0 and s2cnt == 2: genoClass.append(4)
            elif s1cnt == 1 and s2cnt == 1: genoClass.append(5)
            elif s1cnt == 2 and s2cnt == 1: genoClass.append(6)
            elif s1cnt == 1 and s2cnt == 2: genoClass.append(7)
            elif s1cnt == 2 and s2cnt == 2: genoClass.append(8)
    return genoClass

def readSampleConfigFile(configFilename):
    """Read the sample configuration file, and return list of pairs of sample names.

    Read a sample configuration file, with each line representing a sample config,
    and convert this to a tuple of sample names (S1, S2, O). The input file is a
    3 column file, with each column containing a one sample name.

    Args:
        configFilename: Name of configuration file, with 2 names per line, and
            each line representing a configuration (S1, S2, O)

    Returns:
        Array of list of sample names for each configuration in the file.

    Raises:
        IOError: When the format of the input config file name is incorrect.
    """
    configFile = open(configFilename)
    sampleConfigs = []
    for line in configFile:
        toks = line.strip().split()
        if len(toks) != 3:
            raise IOError("Incorrect number of columns in config file.")
        sampleConfigs.append((toks[0], toks[1], toks[2]))
    configFile.close()
    return(sampleConfigs)

def computeSplitTime(counts, generationTime, mutationRate):
    """Compute split time given the counts of configs, generation time and mu.

    Given a count of the different configurations, where sample 1 and sample 2
    have a defined set of genotype combinations, this function computes the two-
    by-two split time using the generation time and mutaion rate.

    Args:
        configFilename: counts of genotype classess.
        generationTime: generation time in number of years.
        mutationRate: mutaion rate in /bp/gen

    Returns:
        List of 2 numbers, t1 and t2.
    """
    f1 = generationTime / float(mutationRate*np.sum(counts))
    f21 = (counts[2]/2.0) + counts[4]
    f22 = ((6*counts[6] + counts[5]) * (2*counts[7] + counts[5])) / (8*counts[5])
    t1 = f1 * (f21 - f22)
    f31 = (counts[1]/2.0) + counts[3]
    f32 = ((2*counts[6] + counts[5]) * (6*counts[7] + counts[5])) / (8*counts[5])
    t2 = f1 * (f31 - f32)
    return t1, t2

def computeBlockJackknifeEstimates(genoCounts):
    """Compute the jackknife estimates of t1 and t2 given the counts of geno classes
    in each block.

    Takes in the lists of genoypte class counts many blocks and returns the total
    estimate, the block jackknife estimate, and the block jackknife estimate of
    the variance.

    Args:
        genoCounts: counts of the genotype classes, one for each block.

    Returns:
        A list of 3 items (each one a pair): the overall estimate of t1 and t2,
        the block jackknife estimate of t1 and t2, and the variance of t1 and t2.

    Raises:
        ValueError: When the genoCounts is empty.
    """
    npairs = np.shape(genoCounts)[1]
    nblocks = np.shape(genoCounts)[0]
    nsnps = np.zeros((nblocks, npairs))
    # count number of snps for each pair in each block
    for b in xrange(nblocks):
        for pair in xrange(npairs):
            nsnps[b,pair] = np.sum(genoCounts[b,pair,:])
    # compute block jack knife for each pair:
    estimatesJackknife = np.zeros((nblocks, npairs, 2))
    brange = np.arange(nblocks)
    for b in xrange(nblocks):
        brangeCur = np.delete(brange, b)
        genoCountsCur = np.sum(genoCounts[brangeCur, :, :], axis=0)
        for pair in xrange(npairs):
            ## assume both gen and mut are 1, so no scaling
            temp = computeSplitTime(genoCountsCur[pair,:], 1, 1)
            estimatesJackknife[b, pair, 0] = temp[0]
            estimatesJackknife[b, pair, 0] = temp[1]
    overallEstimate = np.zeros((npairs, 2))
    for pair in xrange(npairs):
        genoCountsCur = np.sum(genoCounts, axis=0)
        temp = computeSplitTime(genoCountsCur, 1, 1)
        overallEstimate[pair, 0] = temp[0]
        overallEstimate[pair, 1] = temp[1]
    blockEstimate = np.zeros((npairs,2))
    for pair in xrange(npairs):
        totsnps = np.sum(nsnps[:,pair])*1.0
        weights = (1-nsnps[:,pair]/totsnps)
        for index in xrange(2):
            blockEstimate[pair, index] = nblocks*overallEstimate[pair, index] \
                                         - np.sum(weights*estimatesJackknife[:, pair, index])
            pseudoValues = (overallEstimate[pair, index] - (1.0-weights)*estimatesJackknife[:, pair, index])
            pseudoValues = pseudoValues / weights
            blockVariance[pair, index] = np.sum((weights/(1-weights))*((blockEstimate[pair, index] - pseudoValues)**2))/nblocks
    return((overallEstimate, blockEstimate, blockVariance))

def jackknife(fileOfVals, blockSize, outfile):
    """Compute jackknife estimates of the split time for each pair of samples.

    Reads the output file from the split time estimate, compute jack knife estimates
    and overall estimate. Do this for each pair.

    Args:
        fileOfVals: file with value of genoclass at each site, for each pair.
        blockSize: Size of block for the jackknife.
        outfile: file to write jackknife estimates to.

    Returns:
        Jackknife estimate, overall estimate, sdev of the split times.

    Raises:
        None
    """
    valfile = open(fileOfVals)
    line = valfile.readline()
    samps = line.strip().split()[2:]
    npairs = len(samps)
    ## np array with each pair as rows, and each geno class as column
    genoCountsCurBlock = np.zeros((npairs, 9))
    ## array with one element for each pair of samples.
    genoCountsBlock = []
    prevChrom = ""
    prevBlock = -1
    for line in valfile:
        line = line.strip().split()
        genoclasses = [int(x) for x in line[2:]]
        chrom = line[0]
        block = int(int(line[1])/blockSize)
        ## This is new block, so add the
        if chrom != curchrom or block != curblock:
            genoCountBlocks.append(genoCountsCurBlock)
            nsnpsBlock.append(nsnpsCurBlock)
            genoCountsCurBlock = np.zeros((npairs, 9))
            nsnpsCurBlock = np.zeros((npairs))
        else:
            for x in xrange(npairs):
                genoclass = genoclasses[x+2]
                if genoclass == -1: continue
                genoCountsCurBlock[x,genoclass] += 1
                nsnpsCurBlock[x] += 1
        prevChrom = chrom
        prevBlock = block
    genoCountsBlock = np.array(genoCountsBlock)
    # Done with counting genotype classes in blocks.
    # So now call the workhorse function.
    estimatesAndVars = computeBlockJackknifeEstimates(genoCountsBlock)
    outf = open(outfile, "w")
    outf.write("S1\tS2\tO\tT1\tT1.BJK\tSE(T2).BJK\tT2\tT2.BJK\tSE(T2).BJK\n")
    for pair in xrange(npairs):
        snames = samps[pair].split(",")
        outf.write("\t".join(snames)+"\t")
        outf.write(str(np.round(estimatesAndVars[0][pair,0],2))+"\t")
        outf.write(str(np.round(estimatesAndVars[1][pair,0],2))+"\t")
        outf.write(str(np.round(np.sqrt(estimatesAndVars[2][pair,0]),4))+"\t")
        outf.write(str(np.round(estimatesAndVars[0][pair,1],2))+"\t")
        outf.write(str(np.round(estimatesAndVars[1][pair,1],2))+"\t")
        outf.write(str(np.round(np.sqrt(estimatesAndVars[2][pair,1]),4))+"\n")
    outf.close()


def main(args):
    """The boss function. Calls all the subfunctions - engineers, pointy haired
    boss of Dilbert.

    Args:
        args: arguments from the argparser.

    Returns:
        Nothing, expected output from the boss.

    Raises:
        IOError: When file is not found or raised by subfunctions.
    """
    if args.vcf[-3:] == ".gz":
        vcffile = gz.open(args.vcf, 'r')
    else:
        vcffile = open(args.vcf, 'r')
    sampleConfigs = readSampleConfigFile(args.sampleConfigs)

    ## Dictionary to store the array counts.
    ## key is tuple (chr, pos) and value is array
    ## of array of genotype counts - 1 per sample
    ## configuration.
    outfile = open(args.outroot + ".genoCounts", "w")
    outfile.write("Chrom\tPosition\t")
    outfile.write("\t".join([x+","+y+","+z for x,y,z in sampleConfigs]))
    outfile.write("\n")
    reader = vcf.Reader(vcffile)
    for record in reader:
        genoClass = countGenotypeCombinations(record, sampleConfigs,
                                              args.mind, args.maxd)
        outfile.write(record.CHROM+"\t"+record.POS+"\t")
        outfile.write("\t".join(genoClass))
        outfile.write("\n")
    outfile.close()
    jackknife(args.outroot + ".genoCounts", args.blocksize, args.outroot + ".ttout")

if __name__ == "__main__":
    ## Write the argument parser.
    parser = argparse.ArgumentParser(
        description="Compute the Two-by-Two split time statistics using 2 samples in the list of samples in the vcf file.\n")

    parser.add_argument("-v", "--vcf",
        help="path of vcf file to be processed.",
        required=True,
        type=str)

    parser.add_argument("-s", "--sampleConfigs",
        help="file with sample configurations (S1 S2 anc) - one per line.",
        required=True,
        type=str)

    parser.add_argument("-o", "--outroot",
        help="output files root name.",
        required=False,
        default="test",
        type=str)

    parser.add_argument("-d", "--min-depth",
        help="minimum per individual depth at site.",
        dest="mind",
        required=False,
        default=0,
        type=int)

    parser.add_argument("-D", "--max-depth",
        help="maximum per individual depth at site.",
        dest="maxd",
        required=False,
        default=10000,
        type=int)

    parser.add_argument("-g", "--min-genoqual",
        help="minimum genotype quality.",
        required=False,
        default=0,
        type=int)

    parser.add_argument("-b", "--blocksize",
        help="block size for jackknife",
        required=False,
        default=2000000,
        type=int)

    args = parser.parse_args()

    main(args)
