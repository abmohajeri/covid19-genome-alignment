import enum


class Score(enum.Enum):
    gap = -1
    match = +1
    mismatch = -1


genome1_path = "Coronavirus.fasta"
genome2_path = "BetaCoronavirus.fasta"

codes = {
    'A': '00',
    'C': '01',
    'G': '10',
    'T': '11'
}

maximum_subseq_len = 200    # aka k
best_cut_threshold = 50
base_seq_len = 6
num_seq_bits = base_seq_len * 2
middle_fault_tolerance = 10

# Plot Variables
colors = {
    'A': (0, 0, 255),
    'C': (0, 255, 0),
    'T': (225, 150, 7),
    'G': (225, 10, 190),
    '-': (200, 200, 200),
}
gene_parts = 2
font_scale = 0.5
thickness = 1
margin = 30
line_height = 25
character_space = 30
