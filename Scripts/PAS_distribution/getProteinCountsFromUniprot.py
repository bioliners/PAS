#!/usr/bin/python3
#Author: Vadim M. Gumerov
import sys
import urllib.request, urllib.error
import json
import fileinput

USAGE = "The script retrieves protein counts of proteomes from the Uniprot Database by Uniprot proteome ids"

PROTEOMES_URL = "https://rest.uniprot.org/proteomes/"
INPUT_FILE = str(sys.argv[1])
TIMEOUT_FILE = str(sys.argv[2])

def getProteinCounts():
	noDataAnymore = False
	for line in fileinput.input(INPUT_FILE):
		for iteration in range(0, 10):
			proteomeId = line.strip()
			if len(proteomeId):
				try:
					url_ready = PROTEOMES_URL + proteomeId
					result = urllib.request.urlopen(url_ready)
					resultAsJson = json.loads(result.read().decode("utf-8"))
					if len(resultAsJson) == 0:   #No data anymore from this page on
						noDataAnymore = True
						break
					if "name" in resultAsJson:	#404 NotFoundError
						break
				#except ValueError: #504 Gateway timeouts
				except urllib.error.HTTPError:
				#except json.decoder.JSONDecodeError:   #504 Gateway timeouts  From Python 3.5+
					if iteration == 9:
						with open (TIMEOUT_FILE, "a") as timeoutFile:
							timeoutFile.write(proteomeId + " Caused the error\n")
					continue
				extractProteinCounts(resultAsJson, proteomeId)
				break
		if noDataAnymore:
			break
	
def extractProteinCounts(jsonObject, proteomeId):
	proteinCount = jsonObject["proteinCount"]
	geneCount = jsonObject["geneCount"]
	print (proteomeId + "\t" + str(proteinCount) + "\t" + str(geneCount))


getProteinCounts()
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
