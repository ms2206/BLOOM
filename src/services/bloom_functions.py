"""This module provides functions to access NCBI database and BLAST.

"""

from Bio import Entrez, SeqIO           # For queries to NCBI web
from Bio.Blast import NCBIWWW, NCBIXML  # To use BLAST

from config.config_manager import ConfigManager

# Creates instance of ConfigManager to get ncbi parameters for queries from config file
config = ConfigManager()

EMAIL = config.get('ncbi.email')
NCBI_API_KEY = config.get('ncbi.api_key')


def get_NCBI_Sequences(barcode, taxid, email=EMAIL, api_key=NCBI_API_KEY):
    """Gets the barcode sequences uploaded to NCBI for a given taxonomic ID and barcode type.

    Args:
        barcode (str): The name of the barcode.
        organism (str): The taxonomical ID of the organism.
        email (str): Email address to use with NCBI Entrez.
        api_key (str): API key to use with NCBI Entrez (optional).

    Returns:
        A list of tuples. Each tuple contains the fasta header and sequence of a barcode sequence.
    """
    sequences = []
    # Entrez parameters
    Entrez.email = email
    query = f'{barcode}[gene] AND txid{taxid}[organism]'
    # Perform search
    handle = Entrez.esearch(db = 'nucleotide',
                            term = query,
                            api_key = api_key, 
                            usehistory = 'y') 
    search_results = Entrez.read(handle)
    handle.close()
    # Get sequences
    num_seqs = int(search_results['Count']) 
    webenv = search_results['WebEnv']
    query_key = search_results['QueryKey']
    handle = Entrez.efetch(db='nucleotide',
                           rettype='gb',
                           retmode='text',
                           retstart=0,
                           retmax=num_seqs,
                           webenv=webenv,
                           query_key=query_key,
                           idtype='acc')
    # Get barcode sequences found
    for record in SeqIO.parse(handle, 'genbank'):
            for feature in record.features:
                if feature.type == 'gene' and 'gene' in feature.qualifiers:
                    if barcode in feature.qualifiers['gene']:
                        # Get nucleotide sequence
                        gene_seq = str(feature.location.extract(record.seq))
                        # Create fasta header
                        fasta_header = f'>{record.id} | {barcode} | {record.description}'
                        # Create tuple with fasta header and nucleotide sequence
                        datapoint = (fasta_header, gene_seq)
                        # Append tuple to list
                        sequences.append(datapoint)
    handle.close()
    # Return results
    return sequences


def get_taxonomy(taxid, email=EMAIL, api_key=NCBI_API_KEY):
    """Gets the taxonomical lineage of a given taxonomical ID.

    Args:
        taxid (str): The ID of the organism.
        email (str): Email address to use with NCBI Entrez.
        api_key (str): API key to use with NCBI Entrez (optional).

    Returns:
        A list of dictionaries with the taxonomical lineage. Keys: Generic name of the rank.
            Values: [Taxonomical ID, Scientific name, Name of the rank]
    """
    Entrez.email = email
    # Perform data fetching
    handle = Entrez.efetch(db = "taxonomy", 
                                 id=taxid, 
                                 api_key=api_key,
                                 retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    # Extract lineage information
    lineage = records[0]['LineageEx']
    return lineage


def get_taxa_id(organism_name, email=EMAIL, api_key=NCBI_API_KEY):
    """Gets the taxonomy ID from NCBI for a given organism name.
    
    Args:
        organism_name (str): The scientific name of the organism.
        email (str): Email address to use with NCBI Entrez.

    Returns:
        str: Taxonomy ID of the organism, or None if not found.
    """
    Entrez.email = email
    try:
        # Search the taxonomy database
        handle = Entrez.esearch(db="taxonomy", 
                                term=organism_name,
                                api_key=api_key)
        record = Entrez.read(handle)
        handle.close()
        # Get the taxonomy ID
        if record["IdList"]:
            taxa_id = record["IdList"][0]
            return taxa_id
        else:
            print(f"No taxonomy ID found for '{organism_name}'")
            return None
    except Exception as e:
        print(f"Error fetching taxonomy ID: {e}")
        return None
 
    
import re

def blast(barcode, rank, email=EMAIL, api_key=NCBI_API_KEY):
    # Define element to return
    data = []
    # Define query
    query = f'{barcode.get_query()}[gene] AND {rank}[Organism]'
    # Blast
    print("")
    print("Blasting...")
    result_handle = NCBIWWW.qblast(
        program         = "blastn",
        database        = "nt",
        sequence        = str(barcode),
        entrez_query    = query,
        hitlist_size    = 20000,
        megablast       = True,
    )
    # Parse blast results
    blast_records = NCBIXML.parse(result_handle)
    # Extract info and print
    for blast_record in blast_records:
        if blast_record.alignments:
            for index, alignment in enumerate(blast_record.alignments):
                for subindex, hsp in enumerate(alignment.hsps):
                    # Get similiarity percentage
                    similarity_percentage = (hsp.identities / hsp.align_length) * 100
                    # Get name of sample
                    name_match = re.search("\|\s\w+\s\w+", alignment.title)
                    code_match = re.search("\w+\|\w+\|\w+\|\w+", alignment.title)
                    name = name_match.group(0)[2:] if name_match else code_match.group(0)
                    # Data point with info of the alignment
                    datapoint = (index, name, similarity_percentage, hsp.score)
                    data.append(datapoint)
    return filter_data(data)

def filter_data(data):
    unique = {}
    for element in data:
        if element[1] in unique.keys():
            if element[2] > unique[element[1]]:
                unique[element[1]] = element[2]
        else:
            unique[element[1]] = element[2]
    return unique