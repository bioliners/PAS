#!/usr/bin/python3
#Author: Vadim M. Gumerov
import sys, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep

#disable SSL verification to avoid config issues
context = ssl._create_unverified_context()

INPUT_FILE = str(sys.argv[1])
OUTPUT_FILE = str(sys.argv[2])
PAS_DOMAINS = set(["PAS", "PAS_5", "PAS_6", "PAS_2", "PAS_3", "PAS_4", "MEKHLA", "PAS_7", "PAS_8", "PAS_9", "PAS_10", "PAS_11", "CpxA_peri", "MLTR_LBD", "DUF5593", "PAS_12", "AbfS_sensor"])

def get_proteins_of_proteomes():
  counter = 1
  with open(INPUT_FILE) as inputFile, open(OUTPUT_FILE, "w") as outputFile:
    for line in inputFile:
      stripLine = line.strip()
      splitLine = stripLine.split("\t")
      if len(splitLine) >= 2:
        proteinAccession = splitLine[1]
        if len(proteinAccession):
          print (str(counter) +  "\t" + stripLine)
          counter+=1
          url = "https://rest.uniprot.org/uniprotkb/" + proteinAccession + "?format=json"
          get_domains(url, stripLine, outputFile)

def get_domains(url, line, outputFile):
  attempts = 0
  resultReceived = False
  while not resultReceived:
    try:
      sleep(0.5)
      req = request.Request(url, headers={"Accept": "application/json"})
      res = request.urlopen(req, context=context)
      # If the API times out due a long running query
      if res.status == 408:
        # wait just over a minute
        sleep(61)
        # then continue this loop with the same URL
        continue
      elif res.status == 204:
        #no data so leave loop
        break
      payload = json.loads(res.read().decode())
      resultReceived = True
      attempts = 0
    except HTTPError as e:
      if e.code == 408:
        sleep(61)
        continue
      else:
        # If there is a different HTTP error, it will re-try 3 times before failing
        if attempts < 3:
          attempts += 1
          sleep(61)
          continue
        else:
          sys.stderr.write("Could not receive response for protein: " + line + "\n")
          raise e
    
    numberOfPasDomains = 0
    pasDomainVariety = set()
    domainArchitecture = []
    if "uniProtKBCrossReferences" in payload:
      for ref in payload["uniProtKBCrossReferences"]:
        if ref["database"] == "Pfam":
          domain = ref["properties"][0]["value"]
          numOfDuplicates = int(ref["properties"][1]["value"])
          for elem in range(numOfDuplicates):
            domainArchitecture.append(domain)
      for domain in domainArchitecture:
        if domain in PAS_DOMAINS:
          pasDomainVariety.add(domain)
          numberOfPasDomains+=1
      #Sorting domains, so that it would be easier to count domain combinations later
      domainArchitecture.sort()
      outputFile.write(line + "\t" + ",".join(domainArchitecture) + "\t" + str(numberOfPasDomains) + "\t" + str(len(pasDomainVariety)) + "\n")
    else:
      with open ("not_in_Uniprot_anymore.txt", "a") as outputFile2:
        outputFile2.write(line + "\n")

if __name__ == "__main__":
  get_proteins_of_proteomes()
