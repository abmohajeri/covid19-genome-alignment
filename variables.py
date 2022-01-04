import enum

genome1_path = "Coronavirus.fasta"
genome2_path = "BetaCoronavirus.fasta"


class Score(enum.Enum):
    gap = -1
    match = +1
    mismatch = -1


codes = {
    'A': '00',
    'C': '01',
    'G': '10',
    'T': '11'
}

minimum_subseq_len = 100    # aka k
best_cut_threshold = 20
base_seq_len = 6
num_seq_bits = base_seq_len * 2
middle_fault_tolerance = 50
