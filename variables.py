import enum

genome1_path = "alpha_corona.txt"
genome2_path = "delta_corona.txt"


class Score(enum.Enum):
    gap = -1
    match = +1
    mismatch = -1


minimum_subseq = 100
base_seq_length = 6
num_seq_bits = base_seq_length * 2
