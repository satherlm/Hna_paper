'''
eFetch Taxonomy Information Script
-------------------------------------------
Code used to obtain taxonomy information from a csv file of BLAST hits.
Ideally should only be run in its entirety once (to avoid making a lot of ncbi requests).

Was going to extract taxIDs from protein XML files and then use those taxIDs to
eFetch taxonomy entries, from which I could then extract all the relevant lineage info.

However, I couldn't easily obtain the taxIDs from the XML files (in GBSeq_feature_table),
even though you'd think that would be straightforward. I can't figure out which tags/keys
to use to get what I want.

So instead I went with a "simpler" and inelegant option (which might require a bit more manual checking):

1. Read in CSV file (blast hit table--generate this online for however many proteins)
2. Extract protein accessions from CSV file into list
3. Obtain protein entries for each accession (xml format--convert to Python object)
4. Find lineage/taxonomy information associated with each protein
5. Write lineage for each protein to csv file. Can then manipulate/analyze further. 
'''
import csv
from Bio import Entrez

Entrez.email = #EMAIL FOR NCBI API
Entrez.api_key = #API KEY

#Reads in csv file
#Each row becomes a list in python_csv
with open("QC85sma2245-Alignment-HitTable.csv", newline='') as csvfile:
    python_csv = csv.reader(csvfile, delimiter = ",")
    csv_list = list(python_csv)
    
#Making list of only protein accessions from csv_list
#Had to delete empty rows at bottom of Excel file for this to work
accessions_list = [item[1] for item in csv_list]


#Efetching xml file for list of accessions
#Would like to save this to a local file to avoid always requesting from Entrez
handle = Entrez.efetch(db="protein", id=accessions_list, rettype="gp", retmode="xml")

#.parse converts handle into a Python object entry by entry
records = Entrez.parse(handle)

'''
Use .parse instead of .read for large file to prevent memory errors
e.g. records = Entrez.parse(handle)
Can then iterate over records and extract required info (in my case the taxid)
for record in records:
etc.
'''

#Make list of taxonomy information for each accession:

taxonomy = []

for record in records:
    taxonomy.append(record["GBSeq_taxonomy"])

#taxonomy is a list of strings--names in strings separated by semicolons
#Split taxonomy strings into lists themselves to make writing to csv easier:

taxonomy_lists = []

for lineage in taxonomy:
    taxonomy_lists.append(lineage.split(sep = ";"))


#Now to write list of lists to CSV file (so I can manipulate it in Excel):

with open("QC85sma2245-taxonomy.csv", "w", newline='') as csvfile2:
    writer = csv.writer(csvfile2)
    writer.writerows(taxonomy_lists)

#Can now analyse taxonomy information in csv file. 

#Parsing handle into Python and printing it to screen to check:
#print(handle.read())
#handle.close()
