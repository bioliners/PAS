#!/usr/bin/python3
#Author: Vadim M. Gumerov
import sys, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
import collections


#disable SSL verification to avoid config issues
context = ssl._create_unverified_context()

INPUT_FILE = str(sys.argv[1])
OUTPUT_FILE = str(sys.argv[2])

# entryType is one of "UniProtKB reviewed (Swiss-Prot)" or "UniProtKB unreviewed (TrEMBL)"
#{"Proteome1":{"p1isoform1": [("p1isoform1", len<152>, entryType), ("p1isoform2", 565), ...], "p2isoform1":[...], ...}. "Proteome2":{}, ...}
PROTEOME_TO_ISOFORMS = {}
ALL_PROTEINS_ISOFORMS = set()
PROTEOME_TO_SINGLE_PROTEIN = collections.defaultdict(list)

def processIsoforms():
  with open(OUTPUT_FILE, "w") as outputFile:
    for proteome, isoformDict in PROTEOME_TO_ISOFORMS.items():
      for isoform, isoformList in isoformDict.items():
        isoformListSortedByLen = sorted(isoformList, key=lambda x: x[2])
        longestIsoform = isoformListSortedByLen[-1]
        longestIsoformList = map(str, longestIsoform)
        outputFile.write(proteome + "\t" +  "\t".join(longestIsoformList) + "\n")
    #Save those proteins for which no results were obtained from the queried urls
    for proteome, proteinList in PROTEOME_TO_SINGLE_PROTEIN.items():
      for proteinLine in set(proteinList):
        outputFile.write(proteinLine + "\n")

def getIsoformsOfProteoms():
  counter = 1
  with open(INPUT_FILE) as inputFile:
    for line in inputFile:
      lineSplit = line.strip().split("\t")
      #if len(lineSplit) == 5 and lineSplit[4] == "Eukaryota":
      if len(lineSplit) == 4:
        proteomeId = lineSplit[0]
        proteinAccession = lineSplit[1]
        url = "https://rest.uniprot.org/genecentric/search?query=accession%3A" + proteinAccession
        print (str(counter) + "\t" + proteomeId + "\t" + proteinAccession)
        counter+=1
        getIsoforms(url, proteinAccession, proteomeId, line.strip())
        #Save datastructures to output every 500 iteration
        if (counter % 500) == 0:
          try:
            with open(OUTPUT_FILE.split(".")[0] +".json", "w") as outfile:
              json.dump(PROTEOME_TO_ISOFORMS, outfile)
            with open(OUTPUT_FILE.split(".")[0] +"_allIsoforms.json", "w") as outfile:
              json.dump(list(ALL_PROTEINS_ISOFORMS), outfile)
            with open(OUTPUT_FILE.split(".")[0] +"_norResFromUrl.json", "w") as outfile:
              json.dump(PROTEOME_TO_SINGLE_PROTEIN, outfile)
          except:
            print ("Could not save to the file")

def getIsoforms(url, proteinAccession, proteomeId, line):
  if proteinAccession not in ALL_PROTEINS_ISOFORMS:
    attempts = 0
    resultReceived = False
    while not resultReceived:
      try:
        # Don't overload the server, give it time before asking for more
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
          # If there is a different HTTP error, it wil re-try 3 times before failing
          if attempts < 3:
            attempts += 1
            sleep(61)
            continue
          else:
            sys.stderr.write("Could not receive response with url: " + url)
            raise e
      if len(payload["results"]):
        isoformList = []
        canonicalProtAccession = payload["results"][0]["canonicalProtein"]["id"]
        canonicalProtLength = payload["results"][0]["canonicalProtein"]["sequence"]["length"]
        canonicalProtEntryType = payload["results"][0]["canonicalProtein"]["entryType"]
        canonicalProtName = payload["results"][0]["canonicalProtein"]["proteinName"]
        isoformList.append((canonicalProtAccession, canonicalProtName, canonicalProtLength, canonicalProtEntryType))
        ALL_PROTEINS_ISOFORMS.add(canonicalProtAccession)
        #collect isoforms
        if "relatedProteins" in payload["results"][0]:
          for relatedProtein in payload["results"][0]["relatedProteins"]:
            #additional safety check; perhaps unnecessary
            relatedProtAccession = relatedProtein["id"]
            relatedProtLength = relatedProtein["sequence"]["length"]
            relatedProtEntryType = relatedProtein["entryType"]
            relatedProtName = relatedProtein["proteinName"]
            ALL_PROTEINS_ISOFORMS.add(relatedProtAccession)
            isoformList.append((relatedProtAccession, relatedProtName, relatedProtLength, relatedProtEntryType))
        if proteomeId not in PROTEOME_TO_ISOFORMS:
          PROTEOME_TO_ISOFORMS[proteomeId] = {}
        PROTEOME_TO_ISOFORMS[proteomeId][proteinAccession] = isoformList
      else:
        #save the only available protein
        print ("Results have not been obtained: " +  proteomeId + "\t" + proteinAccession)
        PROTEOME_TO_SINGLE_PROTEIN[proteomeId].append(line)
  else:
    print (proteinAccession + " has already been processed")

if __name__ == "__main__":
  getIsoformsOfProteoms()
  processIsoforms()

