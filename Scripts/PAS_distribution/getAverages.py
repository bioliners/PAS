#!/usr/bin/python
#Author: Vadim M. Gumerov
import sys
import fileinput
import collections

USAGE = sys.argv[0] + " inputFile taxSource('NCBI' or 'GTDB') taxLevel(one of these: 'genus', 'family', 'order', 'class', 'phylum')  mode(1 or 2, see the script) taxRepres('named' or 'notNamed') > outputFile"

if len(sys.argv) < 6:
    print ("The number of input parameters should be equal to 5.\n")
    print ("Usage: " + USAGE + "\n")
    print ("The script was not executed. Exiting!")
    exit(1)

INPUT_FILE = str(sys.argv[1])
TAXONOMY_SOURCE = str(sys.argv[2])
TAX_LEVEL = str(sys.argv[3])
#1 or 2. If 1 - input is original file, if 2 - input is a file processed by this script
MODE = int(sys.argv[4])

#"named" or "notNamed" taxonomy:
#Exmaple of the 'named' taxonomy: no rank__cellular organisms;superkingdom__Eukaryota;clade__Sar;clade__Alveolata;class__Dinophyceae;order__Suessiales;family__Symbiodiniaceae;genus__Symbiodinium
#Exampe of the 'notNamed' taxonomy: cellular organisms;Eukaryota;Sar;Alveolata;Dinophyceae;Suessiales;Symbiodiniaceae;Symbiodinium
#Another example of the 'notNamed' taxonomy: d__Bacteria;p__Bdellovibrionota;c__Bdellovibrionia;o__Bdellovibrionales;f__Bdellovibrionaceae;g__Bdellovibrio
TAXONOMY_REPRESENTATION = str(sys.argv[5])


"{'Tax1': [12, 19, 225], ...}"
TAXONOMY_TO_PROTEIN_WITH_ISOFORM_COUNTS = collections.defaultdict(list)
"{'Tax1': [12, 14, 145], ...}"
TAXONOMY_TO_PROTEIN_COUNTS = collections.defaultdict(list)
"{'Tax1': [4, 8, 11], ...}"
TAXONOMY_TO_DOMAIN_COUNTS = collections.defaultdict(list)

"{'Tax1': [4465, 2158, 1111], ...}"
#Importatnly: I subtract from the total gene count PAS containing protein counts to reduce possible selfe-correlation in the downstream analysis
TAXONOMY_TO_TOTAL_GENE_COUNT = collections.defaultdict(list)

"====== Normalized values====="
"{'Tax1': [0.000012, 0.002, 0.0008], ...}"
TAXONOMY_TO_PROTEIN_WITH_ISOFORM_COUNTS_NORMALIZED = collections.defaultdict(list) 
"{'Tax1': [0.000011, 0.0006, 0.00036], ...}"
TAXONOMY_TO_PROTEIN_COUNTS_NORMALIZED = collections.defaultdict(list)
"{'Tax1': [0.0009, 0.003, 0.0025], ...}"
TAXONOMY_TO_DOMAIN_COUNTS_NORMALIZED = collections.defaultdict(list)

"{'Tax1': [12.4<protein count>, 13.5<PAS domain count>, 0.005<normalized protein count>, 0.002<normalized PAS domain count>, 1350<total protein count per genome minus pas proteins>], ...}"
TAXONOMY_TO_AVERAGES = collections.defaultdict(list)

PROTEIN_WITH_ISOFORM_COUNT_COLUMN = 5
PROTEIN_COUNT_COLUMN = 11
DOMAIN_COUNT_COLUMN = 12
PROTEIN_COUNT_NORMALIZED_COLUMN = 13
DOMAIN_COUNT_NORMALIZED_COLUMN = 14
PROTEIN_WITH_ISOFORM_COUNT_NORMALIZED_COLUMN = 15

TOTAL_GENE_COUNT_COLUMN = 10
TAX_LEVEL_FIELD = 2
TAX_COLUMN = 7

def putDataNotNamed():
    for line in fileinput.input(INPUT_FILE):
        if "Uniprot Id" not in line:
            lineSplitted = line.strip().split("\t")
            taxonomyLine = lineSplitted[TAX_COLUMN].strip().split(";")
            if taxonomyLine != "Not found":
                taxString = ";".join(taxonomyLine[:TAX_LEVEL_FIELD])
                protein_with_isoform_count = float(lineSplitted[PROTEIN_WITH_ISOFORM_COUNT_COLUMN])
                protein_count = float(lineSplitted[PROTEIN_COUNT_COLUMN])
                domain_count = float(lineSplitted[DOMAIN_COUNT_COLUMN])

                protein_with_isoform_count_normalized = float(lineSplitted[PROTEIN_WITH_ISOFORM_COUNT_NORMALIZED_COLUMN])
                protein_count_normalized = float(lineSplitted[PROTEIN_COUNT_NORMALIZED_COLUMN])
                domain_count_normalized = float(lineSplitted[DOMAIN_COUNT_NORMALIZED_COLUMN])
                total_gene_count = float(lineSplitted[TOTAL_GENE_COUNT_COLUMN])
                if MODE == 1:
                    #Subtracting pas domain protein counts from the total number of genes to reduce selfe-correlation when correaltion test will be done
                    total_gene_count = total_gene_count - protein_count

                TAXONOMY_TO_PROTEIN_WITH_ISOFORM_COUNTS[taxString].append(protein_with_isoform_count)
                TAXONOMY_TO_PROTEIN_COUNTS[taxString].append(protein_count)
                TAXONOMY_TO_DOMAIN_COUNTS[taxString].append(domain_count)

                TAXONOMY_TO_PROTEIN_WITH_ISOFORM_COUNTS_NORMALIZED[taxString].append(protein_with_isoform_count_normalized)
                TAXONOMY_TO_PROTEIN_COUNTS_NORMALIZED[taxString].append(protein_count_normalized)
                TAXONOMY_TO_DOMAIN_COUNTS_NORMALIZED[taxString].append(domain_count_normalized)

                TAXONOMY_TO_TOTAL_GENE_COUNT[taxString].append(total_gene_count)

def putDataNamed():
    for line in fileinput.input(INPUT_FILE):
        if "Uniprot Id" not in line:
            lineSplitted = line.strip().split("\t")
            taxonomyLine = lineSplitted[TAX_COLUMN].strip()
            if taxonomyLine and taxonomyLine != "Not found":
                #I put ";" before TAX_LEVEL to calculate properly, since there are levels called ex., super-class
                #I put + 1 in order to shoft the beginning from the found semicolon, so that the next "find" will find the correct ending semicolon
                levelBegins = taxonomyLine.find(";"+TAX_LEVEL) + 1 
                levelEnds = taxonomyLine[levelBegins:].find(";")
                if levelBegins != 0 and levelEnds != -1:
                    taxString = taxonomyLine[:levelBegins + levelEnds]
                else:
                    taxString = taxonomyLine

                protein_with_isoform_count = float(lineSplitted[PROTEIN_WITH_ISOFORM_COUNT_COLUMN])
                protein_count = float(lineSplitted[PROTEIN_COUNT_COLUMN])
                domain_count = float(lineSplitted[DOMAIN_COUNT_COLUMN])

                protein_with_isoform_count_normalized = float(lineSplitted[PROTEIN_WITH_ISOFORM_COUNT_NORMALIZED_COLUMN])
                protein_count_normalized = float(lineSplitted[PROTEIN_COUNT_NORMALIZED_COLUMN])
                domain_count_normalized = float(lineSplitted[DOMAIN_COUNT_NORMALIZED_COLUMN])
                total_gene_count = float(lineSplitted[TOTAL_GENE_COUNT_COLUMN])
                if MODE == 1:
                    #Subtracting pas domain protein counts from the total number of genes to reduce selfe-correlation when correaltion test will be done
                    total_gene_count = total_gene_count - protein_count

                TAXONOMY_TO_PROTEIN_WITH_ISOFORM_COUNTS[taxString].append(protein_with_isoform_count)
                TAXONOMY_TO_PROTEIN_COUNTS[taxString].append(protein_count)
                TAXONOMY_TO_DOMAIN_COUNTS[taxString].append(domain_count)

                TAXONOMY_TO_PROTEIN_WITH_ISOFORM_COUNTS_NORMALIZED[taxString].append(protein_with_isoform_count_normalized)
                TAXONOMY_TO_PROTEIN_COUNTS_NORMALIZED[taxString].append(protein_count_normalized)
                TAXONOMY_TO_DOMAIN_COUNTS_NORMALIZED[taxString].append(domain_count_normalized)

                TAXONOMY_TO_TOTAL_GENE_COUNT[taxString].append(total_gene_count)

def getAverages():
    calculateAverages(TAXONOMY_TO_PROTEIN_WITH_ISOFORM_COUNTS)
    calculateAverages(TAXONOMY_TO_PROTEIN_COUNTS)
    calculateAverages(TAXONOMY_TO_DOMAIN_COUNTS)
    calculateAverages(TAXONOMY_TO_PROTEIN_WITH_ISOFORM_COUNTS_NORMALIZED)
    calculateAverages(TAXONOMY_TO_PROTEIN_COUNTS_NORMALIZED)
    calculateAverages(TAXONOMY_TO_DOMAIN_COUNTS_NORMALIZED)
    calculateAverages(TAXONOMY_TO_TOTAL_GENE_COUNT)
    for taxonomy, averageList in TAXONOMY_TO_AVERAGES.items():
       print(taxonomy.replace("no rank__cellular organisms;superkingdom__Eukaryota;", "") + "\t" + "\t".join(averageList))

def calculateAverages(taxonomyToListOfValues, roundValue=False):
    for tax, valueList in taxonomyToListOfValues.items():
        average = sum(valueList)/float(len(valueList))
        TAXONOMY_TO_AVERAGES[tax].append(str(average))
        if roundValue:
            averageRounded = int(round(average))
            TAXONOMY_TO_AVERAGES[tax].append(str(averageRounded))

def main():
    global TAX_COLUMN, PROTEIN_COUNT_COLUMN, PROTEIN_WITH_ISOFORM_COUNT_COLUMN, PROTEIN_WITH_ISOFORM_COUNT_NORMALIZED_COLUMN, DOMAIN_COUNT_COLUMN, PROTEIN_COUNT_NORMALIZED_COLUMN, DOMAIN_COUNT_NORMALIZED_COLUMN, TOTAL_GENE_COUNT_COLUMN, TAX_LEVEL_FIELD
    if MODE == 1:
        if TAXONOMY_SOURCE == "NCBI":
            TAX_COLUMN = 7
        elif TAXONOMY_SOURCE == "GTDB":
            TAX_COLUMN = 8
    elif MODE == 2:
        PROTEIN_WITH_ISOFORM_COUNT_COLUMN = 1
        PROTEIN_COUNT_COLUMN = 2
        DOMAIN_COUNT_COLUMN = 3
        PROTEIN_WITH_ISOFORM_COUNT_NORMALIZED_COLUMN = 4
        PROTEIN_COUNT_NORMALIZED_COLUMN = 5
        DOMAIN_COUNT_NORMALIZED_COLUMN = 6
        TOTAL_GENE_COUNT_COLUMN = 7
        TAX_COLUMN = 0
    else:
        print("MODE parameter should either be 1 or 2")
        exit(1)

    if TAX_LEVEL == "phylum":
        TAX_LEVEL_FIELD = 2
    elif TAX_LEVEL == "class":
        TAX_LEVEL_FIELD = 3
    elif TAX_LEVEL == "order":
        TAX_LEVEL_FIELD = 4
    elif TAX_LEVEL == "family":
        TAX_LEVEL_FIELD = 5
    elif TAX_LEVEL == "genus":
        TAX_LEVEL_FIELD = 6
    
    if TAXONOMY_REPRESENTATION == "named":
        putDataNamed()
    elif TAXONOMY_REPRESENTATION == "notNamed":
        putDataNotNamed()

    getAverages()

main()