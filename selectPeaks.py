## Title: From the gff file with the merged peaks, select the peaks that fullfil certain criteria. Output a bed file with peak locations and fasta file with peak seq
## Author: I. Moustakas, i.m*@lumc.nl

import re
import argparse
import pandas as pd
import numpy
import pybedtools as bt


def CountGroups(row):
    """
    Get the peak name string, separate the sample names, count the occurrence of groups in the sample names
    :param row: An object of type pandas.Series. It is a single row from the gff file. 
    :return: An object of type pandas.Series. A row from gff file plus 4 extra fields with the counts of the number of replicates a peak has been found in for each of the 4 groups
    """

    ## In the gene ID the info of the samples (group, number of replicates) is saved
    ## Extract this info with regex
    sampleNameStr = row["geneID"]
    p = re.compile('gene_id "(.+)"')
    m = p.match(sampleNameStr)
    allPeaks = m.group(1)
    sampleNameList = allPeaks.split(",")
    ## Pattern to capture sample group
    regexGroup = 'S_\d+_[GR,pCREB]+_(\w+)_\d+'
    p2 = re.compile(regexGroup)
    group2ReplList = {}
    for sampleName in sampleNameList:
        m = p2.match(sampleName)
        try:
            groupReplicate = m.group(1)
            group = groupReplicate[0]
            replicate = groupReplicate[1]
            if group not in group2ReplList.keys():
                group2ReplList[group] = [replicate]
            else:
                group2ReplList[group].append(replicate)
        except AttributeError:
            print("Sample Name: {0} does not fit in the format described in regex: {1}".format(sampleName, regexGroup))
            exit("Exiting script, see above message")

    groupCounts = {'1':0, '2':0, '3':0, '4':0}
    for key in group2ReplList.keys():
        counts = len(numpy.unique(group2ReplList[key]))
        groupCounts[key] = counts
    groupCountsSeries = pd.Series(groupCounts)
    ## Rename group "#" to "G#" for improved clarity
    gCS = groupCountsSeries.rename(lambda x: "G" + x)
    rowGroupCounts = row.append(gCS)
    return rowGroupCounts



if __name__ == '__main__':
    ### Parse input args
    parser = argparse.ArgumentParser(description='Read the .gff file with the peaks. Filter on group and numbr of replicates')
    parser.add_argument('--gff', type=str, help='gff file', required=True)
    parser.add_argument('--group', type = str, help = 'Group name to filter on', required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--equal_or_greater_than', type=int, help="Number of replicates is equal to or greater than this integer")
    group.add_argument('--equal_to', type=int, help="Number of replicates is exactly equal to this integer")
    args = parser.parse_args()

    ## read gff, name fields accordingly
    headerNames = [ "chr", "source", "method", "start", "end", "score", "strand", "phase", "geneID"]
    gff = pd.read_csv(args.gff, sep = "\t", skiprows=2, header=None, names=headerNames)

    ## Get the trancription factor name from the gff file
    fileName = args.gff.split("/")[-1]
    tfName = fileName.split("_")[0]

    ## Apply counting function
    gffGroupCounts = gff.apply(CountGroups, axis=1, reduce=False)

    group = args.group
    if args.equal_or_greater_than is not None:
        numberRepl = args.equal_or_greater_than
        slice = gffGroupCounts[ gffGroupCounts[group] >= numberRepl ]
        fileName = "_{0}_replicatesEqualOrGreaterThan_{1}".format(group, numberRepl)
    else:
        numberRepl = args.equal_to
        slice = gffGroupCounts[gffGroupCounts[group] == numberRepl]
        fileName = "_{0}_replicatesEqualTo_{1}".format(group, numberRepl)

    bed = slice[["chr", "start", "end"]]
    bed.to_csv( tfName+fileName+".bed", sep="\t", header=False, index=False)
    bedToString = bed.to_string(header = False, index = False)
    bedString = bt.BedTool(bedToString, from_string=True)
    fasta = bt.example_filename("/exports/genomes/species/R.norvegicus/Rnor_6.0/reference.fa")
    bedString.sequence(fi=fasta, fo= tfName+fileName+".fasta")

