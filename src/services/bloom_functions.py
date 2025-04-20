"""This module provides functions to access NCBI database and BLAST.

"""

from Bio import Entrez, SeqIO           # For queries to NCBI web
from Bio.Blast import NCBIWWW, NCBIXML  # To use BLAST
from config.config_manager import ConfigManager
import re

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
 

def blast(sequence, barcode_query, rank, megablast_use, email=EMAIL, api_key=NCBI_API_KEY):
    # Define element to return
    data = []
    # Define query
    query = f'({barcode_query}[gene] OR (complete genome[all])) AND {rank}[Organism]'
    # Blast
    result_handle = NCBIWWW.qblast(
        program         = "blastn",
        database        = "nt",
        sequence        = sequence,
        entrez_query    = query,
        hitlist_size    = 20000,
        megablast       = megablast_use,
    )
    # Parse blast results
    blast_records = NCBIXML.parse(result_handle)
    # Extract info and print
    for blast_record in blast_records:
        if blast_record.alignments:
            for alignment in blast_record.alignments:
                # Get the scientific name of the organism aligned
                pattern = r"(?:\[\s*)?([A-Z][a-z]+)(?:\s*\])?\s+([a-z\-]+)"
                matches = re.findall(pattern, alignment.title)
                if not matches:
                    continue # When the name is like A.sativa skip it. Next time write the full name >:(
                name = " ".join(matches[0])
                # Check only the first hsp and get parameters
                best_hsp = alignment.hsps[0]
                identities = best_hsp.identities
                gaps = best_hsp.gaps
                diffs = len(sequence) - identities + gaps
                score = best_hsp.score
                sbjct_seq = best_hsp.sbjct
                # Get percentage of identity 
                identity_pct = 100*score/len(sequence)
                # Create datapoint to save
                datapoint = {
                    "Title": alignment.title,
                    "Scientific name": name,
                    "Hit ID": alignment.hit_id,
                    "Accession Number": alignment.accession,
                    "Identities": identities,
                    "Differences": diffs,
                    "Identity percentage": identity_pct,
                    "Score": score,
                    "Subject sequence": sbjct_seq
                }
                # Save datapoint
                data.append(datapoint)
    return data


def filter_data2(data, target_key):
    """Filters the high scoring pairs to get just the ones with the best score."""
    print("Filtering data")
    unique = {}
    for element in data:
        name = element["Scientific name"]
        num_differences = element[target_key]
        if name in unique.keys():
            if num_differences < unique[name]:
                unique[name] = num_differences
        else:
            unique[name] = num_differences
    print("done")
    return unique

def filter_data(data):
    """Filters the data to get only on entry per scientific name"""
    unique = {}
    for element in data:
        name = element["Scientific name"]
        identities = element["Identities"]
        if name in unique.keys():
            if identities > unique[name]["Identities"]:
                unique[name] = element
        else:
            unique[name] = element
    return list(unique.values())