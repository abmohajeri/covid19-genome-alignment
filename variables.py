import enum

genome1_path = "alpha_corona.txt"
genome2_path = "delta_corona.txt"


class Score(enum.Enum):
    gap = -1
    match = +1
    mismatch = -1


minimum_subseq_len = 100
best_cut_threshold = 20
base_seq_len = 6
num_seq_bits = base_seq_len * 2
middle_fault_tolerance = 5
