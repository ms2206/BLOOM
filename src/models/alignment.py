"""Contaings aligners configurations and functions to align nucleotide sequences.

"""

from Bio import Align
from Bio.Seq import Seq

from config.config_manager import ConfigManager

# Creates instance of ConfigManager to get aligner values from config file
config = ConfigManager()

""" Aligner to detect primers in sequences. It prefers mismatches over gaps."""
primer_aligner = Align.PairwiseAligner(mode = 'local',
                                       match_score = config.get('primer_aligner.match_score'),
                                       mismatch_score = config.get('primer_aligner.mismatch_score'),
                                       open_gap_score = config.get('primer_aligner.open_gap_score'),
                                       extend_gap_score = config.get('primer_aligner.extend_gap_score'))

""" Aligner to compare two sequences of similar length. It prefers gaps over mismatches."""
long_aligner = Align.PairwiseAligner(mode = 'global',
                                     match_score = config.get('long_aligner.match_score'),
                                     mismatch_score = config.get('long_aligner.mismatch_score'),
                                     open_gap_score = config.get('long_aligner.open_gap_score'),
                                     extend_gap_score = config.get('long_aligner.extend_gap_score'))


def find_primers(sequence, primers):
    """ Finds the starting and ending position of each primer in the sequence.

    Args:
        sequence (str): nucleotide sequence
        primers (tuple): tuple of sequences for the starting primer and ending primer (in that order)
    
    Returns:
        List of tuples with the intervals where the nucleotide sequence is the same as the primer
    """
    coordinates = []
    for i, primer in enumerate(primers):
        primer_seq = Seq(primer) if not i else Seq(primer).reverse_complement()
        best_score = 0
        best_coords = [(0, 0)] if not i else [(len(sequence), len(sequence))] # Default values if there is no alignment
        alignments = primer_aligner.align(sequence, primer_seq)
        for alignment in alignments:
            gaps = alignment.counts()[0]
            mismatches = alignment.counts()[2]
            if alignment.score > best_score and (gaps + mismatches) < 10:
                best_coords = get_correct_intervals(alignment)
                best_score = alignment.score
        coordinates.append(best_coords)
    return coordinates

def get_seqs_alignment(seq1, seq2):
    """Gets the alignment object for two long sequences."""
    alignments = long_aligner.align(seq1, seq2)
    best_score = 0
    best_alignment = None
    for i, alignment in enumerate(alignments):
        # If there are a lot of alignments to check just look at the first 1000
        if i > 100:
            break
        else:
            if alignment.score > best_score:
                best_alignment = alignment
    return best_alignment
    

def get_correct_intervals(alignment):
    """
    Given a pairwise alignment object (from Biopython’s align package)
    where the first (top) sequence is considered the reference,
    return a list of intervals (tuples of integer indexes) corresponding to
    regions where the reference sequence is "correct" (i.e. the top character matches
    the bottom character). Gaps in the top sequence are skipped since they do not
    correspond to positions in the real (ungapped) reference sequence.
    
    For example, for an alignment with:
    
    Top:    CG-AATCGAGTTA-A--TTACG
    Marker: ||-|||||-|-||-|--.||||
    Bottom: CGAAATCG-G-TAGACGCTACG
    
    the result is:
    
      [(0, 6), (8, 8), (10, 12), (14, 17)]
    """
    # Get the aligned sequences as strings.
    top_seq = str(alignment[0])
    bottom_seq = str(alignment[1])

    intervals = []  # list to hold (start, end) intervals of correct positions.
    current_start = None  # start coordinate for the current interval
    ref_index = alignment.coordinates[0][0]  # coordinate in the ungapped reference sequence

    # Iterate over alignment columns.
    for top_char, bottom_char in zip(top_seq, bottom_seq):
        # If the top sequence has a gap, it does not correspond to a real base.
        if top_char == '-':
            continue
        
        # Check if the base is correct (a match between top and bottom).
        if top_char == bottom_char:
            if current_start is None:
                # Start a new interval.
                current_start = ref_index
        else:
            # Mismatch: if an interval is open, close it.
            if current_start is not None:
                intervals.append((current_start, ref_index - 1))
                current_start = None

        # Increment the reference coordinate only when top_char is not a gap.
        ref_index += 1

    # Close any open interval at the end.
    if current_start is not None:
        intervals.append((current_start, ref_index - 1))
    
    return intervals