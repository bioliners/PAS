#!/usr/bin/python3
#Author: Vadim M. Gumerov
import sys, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep

#disable SSL verification to avoid config issues
context = ssl._create_unverified_context()

INPUT_FILE = str(sys.argv[1])

def get_proteins_of_proteomes():
  with open(INPUT_FILE) as inputFile:
    for line in inputFile:
      proteomeId = line.strip().split("\t")[0]
      if len(proteomeId):
        url = "https://www.ebi.ac.uk/interpro/api/protein/UniProt/set/pfam/CL0183/proteome/uniprot/" + proteomeId + "/?extra_fields=counters&amp;format=json&amp;page_size=20"
        get_proteins(url, proteomeId)

def get_proteins(url, proteomeId):
  next = url
  attempts = 0
  while next:
    try:
      req = request.Request(next, headers={"Accept": "application/json"})
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
      next = payload["next"]
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
          sys.stderr.write("LAST URL: " + next)
          raise e
    for i, item in enumerate(payload["results"]):
      uniprotAccession = item["metadata"]["accession"]
      proteinName = item["metadata"]["name"]
      proteinLength = str(item["metadata"]["length"])
      sys.stdout.write(proteomeId + "\t" + "\t".join([uniprotAccession, proteinName, proteinLength]) + "\n")
         
    # Don't overload the server, give it time before asking for more
    if next:
      sleep(1)


if __name__ == "__main__":
  get_proteins_of_proteomes()

