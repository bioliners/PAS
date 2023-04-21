#!/usr/bin/python3
#Author: Vadim M. Gumerov
import re
import sys
import fileinput
from Bio import Entrez
import time
Entrez.email = 'gumerov.1@osu.edu'

try:
	inputFile = str(sys.argv[1])
	# Two options: prok [menaing prokaryote], euk [meaning eukaryote]
	prokOrEuk = "prok"
	if len(sys.argv) == 3:
		prokOrEuk = str(sys.argv[2])
	for line in fileinput.input(inputFile):
		taxId = line.strip()
		if len(taxId):
			try: 
				handle = Entrez.efetch(db='taxonomy', id=taxId, retmode='xml')
				record = Entrez.read(handle, validate=False)
				handle.close()

				tax_list = record[0]['LineageEx']
				taxStr = ""
				if prokOrEuk == "prok":
					for tax_element in tax_list:
						taxStr+=tax_element['ScientificName'] + ";"
				elif prokOrEuk == "euk":
					for tax_element in tax_list:
						taxStr+=tax_element['Rank'] + "__" + tax_element['ScientificName'] + ";"
				print (taxId + "\t" + taxStr.strip(";"))
			except:
				print (sys.exc_info()[0])
				print (taxId + " raised Exception")
		time.sleep (0.4)
except:
	print (sys.exc_info()[0])
	print (taxId + " aised Exception. Outer Try.")

finally:
	fileinput.close()




