#!/usr/bin/python3
#Downloader from InterPro with my changes
# standard library modules
import sys, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep

BASE_URL = "https://www.ebi.ac.uk:443/interpro/api/proteome/uniprot/set/pfam/CL0183/?extra_fields=counters,assembly&page_size=100"

def parse_items(items):
  if type(items)==list:
    return ",".join(items)
  return ""
def parse_member_databases(dbs):
  if type(dbs)==dict:
    return ";".join([f"{db}:{','.join(dbs[db])}" for db in dbs.keys()])
  return ""
def parse_go_terms(gos):
  if type(gos)==list:
    return ",".join([go["identifier"] for go in gos])
  return ""
def parse_locations(locations):
  if type(locations)==list:
    return ",".join(
      [",".join([f"{fragment['start']}..{fragment['end']}" 
                for fragment in location["fragments"]
                ])
      for location in locations
      ])
  return ""
def parse_group_column(values, selector):
  return ",".join([parse_column(value, selector) for value in values])

def parse_column(value, selector):
  if value is None:
    return ""
  elif "member_databases" in selector:
    return parse_member_databases(value)
  elif "go_terms" in selector: 
    return parse_go_terms(value)
  elif "children" in selector: 
    return parse_items(value)
  elif "locations" in selector:
    return parse_locations(value)
  return str(value)

def output_list():
  #disable SSL verification to avoid config issues
  context = ssl._create_unverified_context()

  next = BASE_URL
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
      sys.stdout.write(parse_column(item["metadata"]["accession"], 'metadata.accession') + "\t")
      sys.stdout.write(parse_column(item["metadata"]["source_database"], 'metadata.source_database') + "\t")
      sys.stdout.write(parse_column(item["metadata"]["name"], 'metadata.name') + "\t")
      sys.stdout.write(parse_column(item["metadata"]["taxonomy"], 'metadata.taxonomy') + "\t")
      sys.stdout.write(parse_column(item["extra_fields"]["counters"]["entries"], 'extra_fields.counters.entries') + "\t")
      sys.stdout.write(parse_column(item["extra_fields"]["counters"]["proteins"], 'extra_fields.counters.proteins') + "\t")
      sys.stdout.write(parse_column(item["extra_fields"]["assembly"], 'extra_fields.assembly') + "\t")     
      sys.stdout.write("\n")
      
    # Don't overload the server, give it time before asking for more
    if next:
      sleep(1)


if __name__ == "__main__":
  output_list()

