import numpy as np
import variables
from genome_hash_table import HashTable


def nw(genome1, genome2, match=variables.Score.match.value, mismatch=variables.Score.mismatch.value,
       gap=variables.Score.gap.value):
    nx = len(genome1)
    ny = len(genome2)
    f = np.zeros((nx + 1, ny + 1))
    f[:, 0] = np.linspace(0, -nx, nx + 1)
    f[0, :] = np.linspace(0, -ny, ny + 1)
    # Pointers to trace through an optimal alignment.
    p = np.zeros((nx + 1, ny + 1))
    p[:, 0] = 3
    p[0, :] = 4
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if genome1[i] == genome2[j]:
                t[0] = f[i, j] + match
            else:
                t[0] = f[i, j] + mismatch
            t[1] = f[i, j + 1] + gap
            t[2] = f[i + 1, j] + gap
            tmax = np.max(t)
            f[i + 1, j + 1] = tmax
            if t[0] == tmax:
                p[i + 1, j + 1] += 2
            if t[1] == tmax:
                p[i + 1, j + 1] += 3
            if t[2] == tmax:
                p[i + 1, j + 1] += 4
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if p[i, j] in [2, 5, 6, 9]:
            rx.append(genome1[i - 1])
            ry.append(genome2[j - 1])
            i -= 1
            j -= 1
        elif p[i, j] in [3, 5, 7, 9]:
            rx.append(genome1[i - 1])
            ry.append('-')
            i -= 1
        elif p[i, j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(genome2[j - 1])
            j -= 1
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    score = f[-1, -1]
    return score


def read_genomes(genome_path):
    genome = []
    with open(genome_path) as f:
        while True:
            # Read from file
            c = f.read(1)
            if not c:
                break
            elif c == '\n':
                continue
            else:
                genome.append(c)
    return genome


def create_genome_arrays():
    genome1 = read_genomes(variables.genome1_path)
    genome2 = read_genomes(variables.genome2_path)
    return genome1, genome2


def create_hash_table(genome1):
    ht = HashTable(2 ** variables.num_seq_bits)
    for i in range(len(genome1) - 5):
        substr = genome1[i:i + 6]
        ht.add(substr, i)
    return ht


# extend seed from back and front
def seed_n_extend(genome1, genome2, start1, start2):
    forward_equality = variables.base_seq_len
    backward_equality = 0
    genome1_forward_subseq = genome1[start1:start1 + forward_equality]
    seed_forward = genome2[start2:start2 + forward_equality]
    genome1_backward_subseq = genome1[start1:start1 + forward_equality]
    seed_backward = genome2[start2:start2 + forward_equality]

    while genome1_forward_subseq == seed_forward or genome1_backward_subseq == seed_backward:
        equality = forward_equality + 1
        backward_equality = backward_equality + 1
        genome1_forward_subseq = genome1[start1:start1 + equality]
        seed_forward = genome2[start2:start2 + equality]
        genome1_backward_subseq = genome1[start1 - backward_equality:start1 + equality]
        seed_backward = genome2[start2 - backward_equality:start2 + equality]
        if forward_equality - 1 >= variables.best_cut_threshold:
            return True, start1, start2
        elif backward_equality - 1 >= variables.best_cut_threshold - variables.base_seq_len:
            return True, start1 - backward_equality, start2 - backward_equality
    return False, -1, -1


def reboot_best_cut(genome2, genome2_mid):
    seed = genome2[genome2_mid + 1:genome2_mid + variables.base_seq_len + 1]
    genome2_mid = genome2_mid + 1
    return seed, genome2_mid


def find_best_cut(genome1, genome2):
    genome2_mid = int(len(genome2) / 2)
    initial_genome2_cut = genome2_mid
    initial_genome1_cut = len(genome1) / 2
    seed = genome2[genome2_mid:genome2_mid + variables.base_seq_len]

    while genome2_mid <= len(genome2) - variables.base_seq_len:
        found, genome1_mid = ht.get(seed, len(genome1))
        if found:
            best_cunt_found, genome1_cut, genome2_cut = seed_n_extend(genome1, genome2, genome1_mid, genome2_mid)
            if best_cunt_found:
                return int(genome1_cut), int(genome2_cut)
            else:
                seed, genome2_mid = reboot_best_cut(genome2, genome2_mid)
                continue
        else:
            seed, genome2_mid = reboot_best_cut(genome2, genome2_mid)
            continue
    return int(initial_genome1_cut), int(initial_genome2_cut)


def divide_n_conquer_seq_align(genome1, genome2):
    global total_score
    if len(genome1) <= variables.minimum_subseq_len and len(genome2) <= variables.minimum_subseq_len:
        return nw(genome1, genome2)
    else:
        genome1_cut, genome2_cut = find_best_cut(genome1, genome2)
        return divide_n_conquer_seq_align(genome1[0:genome1_cut], genome2[0:genome2_cut]) + \
                      divide_n_conquer_seq_align(genome1[genome1_cut:-1], genome2[genome2_cut:-1])


genome1, genome2 = create_genome_arrays()
ht = create_hash_table(genome1)


def main():
    print(divide_n_conquer_seq_align(genome1, genome2))


if __name__ == '__main__':
    main()
