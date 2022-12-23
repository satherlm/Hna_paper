'''
Genbank Neighbours Finder
-----------------------------------
get_neighbour_info extracts locus tag, product description, and protein id for the ten coding sequences upstream and downstream (10 each)
of a cds of interest (specified by locus tag) from Genbank file and saves it to a csv file (named with locus tag).

get_all_neighbourhoods can do this automatically for a given list of GenBank files and specified locus tags
of interest (tags read in from a csv file (NON UTF-8 encoded by Excel)).

Once you have protein ids for neighbouring cds, it is relatively easy to obtain amino acid sequences for further analysis (e.g. by EFetch).
'''

from Bio import SeqIO
import csv

def get_coding_seqs(file):
    '''Reads a Genbank file into a SeqRecord object and creates a list of all CDS features in the record. Assumes single record for each file.'''
    record = SeqIO.read(file, "genbank")
    coding_seqs = [feature for feature in record.features if feature.type == 'CDS']
    return coding_seqs

def get_index_of_interest(coding_seqs, tag):
    '''Returns index of cds in coding_seqs list with specified locus tag'''
    index = [coding_seqs.index(cds) for cds in coding_seqs if cds.qualifiers['locus_tag'] == [tag]]
    return index[0]

def get_neighbour_info(coding_seqs, index_of_interest):
    '''Returns locus tag, product, protein id for up- and downstream neighbours (10 +/-) of indexed cds as a list of lists'''
    if len(coding_seqs) < index_of_interest + 9: #if there are not 10 cds ahead... 
        indices_down = range(index_of_interest + 1, len(coding_seqs))
        print("Length of coding_seqs:", len(coding_seqs))
        print(indices_down)
    else:
        indices_down = range(index_of_interest + 1, index_of_interest + 11)
        print("Length of coding_seqs:", len(coding_seqs))
        print(indices_down)

    #Generate range of indices of 10 cds upstream:
    if index_of_interest > 10:
        indices_up = range(index_of_interest - 10, index_of_interest)
        print("Index of interest:", index_of_interest)
        print(indices_up)
    else:
        indices_up = range(0, index_of_interest)

    neighbours = []

    for i in indices_up:
        locus_tag = coding_seqs[i].qualifiers['locus_tag'][0] #The 0 index is to extract string instead of single-item list
        product = coding_seqs[i].qualifiers['product'][0]
        protein_id = coding_seqs[i].qualifiers.get('protein_id', 'No id given') #Trying to handle dictionary KeyErrors... 
        #protein_id = coding_seqs[i].qualifiers['protein_id'][0]
        #neighbours.append([locus_tag, product, protein_id])
        neighbours.append([locus_tag, product, protein_id])

    for i in indices_down:
        locus_tag = coding_seqs[i].qualifiers['locus_tag'][0] #The 0 index is to extract string instead of single-item list
        product = coding_seqs[i].qualifiers['product'][0]
        protein_id = coding_seqs[i].qualifiers.get('protein_id', 'No id given')
        #protein_id = coding_seqs[i].qualifiers['protein_id'][0]
        #neighbours.append([locus_tag, product, protein_id])
        neighbours.append([locus_tag, product, protein_id])

    return neighbours


def genbank_neighbourhood(file, locus_tag):
    '''General function that calls the earlier ones--fetches neighbourhood info from file for specific locus tag and writes it to a csv file.
Input path to file and locus tag of interest and this does the rest.'''
    coding_seqs = get_coding_seqs(file)
    index = get_index_of_interest(coding_seqs, locus_tag)
    neighbours_info = get_neighbour_info(coding_seqs, index)

    #Writing results to file:
    with open(locus_tag, "w", newline = '') as file_to_write:
        writer = csv.writer(file_to_write)
        writer.writerows(neighbours_info)
        
    return print("Neighbour data for locus {} from {} written to file {}".format(locus_tag, file, locus_tag))
    #Is the print necessary here?? 

def get_all_neighbourhoods(file_names, locus_tags_file):
    '''Extracts neighbourhood data for the locus tag specified for each file in a given list.
Reads in a csv file with the locus tags of interest in the same order as the list of files'''
    #Read locus tags into a list from specified locus_tags_file (a csv file)
    with open(locus_tags_file, newline = '') as csv_file:
        python_csv = csv.reader(csv_file, dialect = 'excel', delimiter = ",")
        csv_list = list(python_csv)
        locus_tags_list = [item[0] for item in csv_list]
        print(locus_tags_list)
        
    i = 0
    for file in file_names:
        genbank_neighbourhood(file, locus_tags_list[i])
        i = i + 1
    return "Neighbour data extracted from {} GenBank files".format(i)

'''
Created list of genbank file names:
>>> species_list = []
>>> for x in range(1, 22):
	species_list.append("species_{}.gb".format(x))

Called function as follows:
>>> get_all_neighbourhoods(species_list, "locus_tags.csv")
'''
