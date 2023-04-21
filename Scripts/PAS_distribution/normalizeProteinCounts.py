#!/usr/bin/python
#Author: Vadim M. Gumerov
import sys
import fileinput

USAGE = sys.argv[0] + " inputFile > outputFile"

INPUT_FILE = str(sys.argv[1])
if INPUT_FILE == "-h" or INPUT_FILE == "--help":
    print(USAGE)
    exit(0)

def getNormalizedCounts():
    for line in fileinput.input(INPUT_FILE):
        if "Uniprot Id" not in line:
            lineSplitted = line.strip().split("\t")
            if len(lineSplitted) == 13:
                #fullProteinCount = int(lineSplitted[9])
                fullGeneCount = int(lineSplitted[10])
                pasProteinWithIsoformCount = int(lineSplitted[5])
                pasProteinCount = int(lineSplitted[11])
                normalizedPasProteinCount = pasProteinCount/float(fullGeneCount)
                normalizedPasProteinWithIsoformCount = pasProteinWithIsoformCount/float(fullGeneCount)
                pasDomainCount = int(lineSplitted[12])
                normalizedPasDomainCount= pasDomainCount/float(fullGeneCount)
                print ("\t".join([line.strip(), str(normalizedPasProteinCount), str(normalizedPasDomainCount), str(normalizedPasProteinWithIsoformCount)]))

def main():
    getNormalizedCounts()

main()