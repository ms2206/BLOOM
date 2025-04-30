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
        taxis (str): The taxonomical ID of the organism.
        email (str): Email address to use with NCBI Entrez.
        api_key (str): API key to use with NCBI Entrez (optional).

    Returns:
        A list of tuples. Each tuple contains the fasta header and sequence of a barcode found in 
        NCBI database.
    """
    sequences = []
    # Entrez parameters
    Entrez.email = email
    query = f'{barcode}[gene] AND txid{taxid}[organism]'
    # Perform search
    try:
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
    except Exception as e:
        # Raise exception
        return (None, e)
    else:
        # Return results
        return sequences
    finally:
        handle.close()


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
    try:
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
    except:
        return None


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
        taxa_id = record["IdList"][0]
        return taxa_id
    except:
        return None
 

def blast(sequence, barcode_query, rank, megablast_use):
    """ Returns BLAST results as a list of dictionaries for the alignment of the target barcode sequence 
    with species in the specified taxonomy rank.

    Args:
        sequence (str): Nucleotide sequence to blast
        barcode_query (str): name of the barcode to filter fasta documents in organisms blasted
        rank (str): biological name of the taxonomical rank of the organism that contains the desired species
            to blast
        megablast_use (bool): True if the users want to use Megablast or False otherwise

    Returns:
        List with dictionaries with information for every hit that BLAST reports
    """
    # Define element to return
    data = []
    # Define query
    query = f'({barcode_query}[gene] OR (complete genome[all])) AND {rank}[Organism]'
    try:
        # Blast
        result_handle = NCBIWWW.qblast(
            program         = "blastn",
            database        = "nt",
            sequence        = sequence,
            entrez_query    = query,
            hitlist_size    = 20000,            # Increase hitlist to avoid falling short
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
                    # Get percentage of identity with the score
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
    except Exception as e:
        return (None, e)
    else:
        return data


def filter_data(data):
    """
    Filters the data obtained from BLAST to get only one entry per scientific name. 
    Results are filtered by identity percentage, the data corresponding to the organism kept will
    be the one with the higher identity percentage as it is the worst-case scenario.
    """
    # Dictionary to store unique species names
    unique = {}
    for element in data:
        name = element["Scientific name"]
        identity_pct = element["Identity percentage"]
        # Check if the name has already been stored
        if name in unique.keys():
            # Store the one with the greatest identity percentage
            if identity_pct > unique[name]["Identity percentage"]:
                unique[name] = element
        else:
            unique[name] = element
    return list(unique.values())