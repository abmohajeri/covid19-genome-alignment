import numpy as np
import math
import variables


def nw(genome1, genome2, match=variables.Score.match.value, mismatch=variables.Score.mismatch.value, gap=variables.Score.gap.value):
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
            t[1] = f[i, j+1] + gap
            t[2] = f[i+1, j] + gap
            tmax = np.max(t)
            f[i+1, j+1] = tmax
            if t[0] == tmax:
                p[i+1, j+1] += 2
            if t[1] == tmax:
                p[i+1, j+1] += 3
            if t[2] == tmax:
                p[i+1, j+1] += 4
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
    return '\n'.join([rx, ry])


class Node:
    def __init__(self, value, next=None):
        self.value = value
        self.next = next


class LinkedList:
    def __init__(self):
        self.head = None
        self.size = 0

    def insert_at_first(self, value):
        self.head = Node(value, self.head)
        self.size += 1

    def insert_at_last(self, value):
        current = self.head
        while current.next:
            current = current.next
        current.next = Node(value)
        self.size += 1

    def find(self, value):
        current = self.head
        while current:
            if value + variables.middle_fault_tolerance >= current.value >= value - variables.middle_fault_tolerance:
                return True, current.value
            current = current.next
        return False, -1

    def print_all(self, idx):
        current = self.head
        print("index: ", idx, "values: ", end=' ')
        while current:
            print(current.value, end=' ')
            current = current.next
        print()


class HashTable:
    def __init__(self, size):
        self.table = [None] * size

    @staticmethod
    def hash_key(sequence) -> int:
        codes = {
            'a': '00',
            'c': '01',
            'g': '10',
            't': '11'
        }
        chars = list(sequence)
        binary_str = ''
        for char in chars:
            binary_str += codes[char]
        return int(binary_str, 2)

    def add(self, sequence, value):
        idx = self.hash_key(sequence)
        if self.table[idx] is None:
            new_linked_list = LinkedList()
            new_linked_list.insert_at_first(value)
            self.table[idx] = new_linked_list
        else:
            self.table[idx].insert_at_last(value)

    def get(self, sequence, genome1_len) -> any:
        idx = self.hash_key(sequence)
        if self.table[idx] is not None:
            return self.table[idx].find(genome1_len/2)
        else:
            return False

    def print_all(self):
        for i in range(len(self.table)):
            if self.table[i]:
                self.table[i].print_all(i)

# x = "AGATTCGATTACAAGAGATTACAGATTACAAATTAGATTCGATTACAAGAGATTACAGATTACAAATTAGATTCGATTACAAGAGATTACAGATTACAAATTAGATTCGATTACAAGAGATTACAGATTACAAATT"
# y = "AGATTCGCAGGATTACAAGATGATTACAATACAA"
# ht = HashTable(math.ceil(len(x)/4))
# for i in range(len(x) - 5):
#     substr = x[i:i+6]
#     ht.add(substr, i)
# ht.print_all()


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
    ht = HashTable(2**variables.num_seq_bits)
    for i in range(len(genome1) - 5):
        substr = genome1[i:i+6]
        ht.add(substr, i)
    return ht

# extend seed from back and front
def seed_n_extend(genome1, genome2, start1, start2):
    forward_equality = variables.base_seq_len
    backward_equality = 0
    genome1_forward_subseq = genome1[start1, start1 + forward_equality]
    seed_forward = genome2[start2, start2 + forward_equality]
    genome1_backward_subseq = genome1[start1, start1 + forward_equality]
    seed_backward = genome2[start2, start2 + forward_equality]

    while genome1_forward_subseq == seed_forward or genome1_backward_subseq == seed_backward:
        equality = forward_equality + 1
        backward_equality = backward_equality + 1
        genome1_forward_subseq = genome1[start1, start1+equality]
        seed_forward = genome2[start2, start2 + equality]
        genome1_backward_subseq = genome1[start1 - backward_equality, start1+equality]
        seed_backward = genome2[start2 - backward_equality, start2 + equality]

    if forward_equality - 1 >= variables.best_cut_threshold:
        return True, start1, start2
    elif backward_equality - 1 >= variables.best_cut_threshold - variables.base_seq_len:
        return True, start1 - backward_equality, start2 - backward_equality
    else:
        return False, -1, -1

def reboot_best_cut(genome2, genome2_mid):
    seed = genome2[genome2_mid + 1, genome2_mid + variables.base_seq_len + 1]
    genome2_mid = genome2_mid + 1
    return seed, genome2_mid

def find_best_cut(ht, genome1, genome2):
    genome2_mid = len(genome2) / 2
    genome1_cut = 0
    genome2_cut = 0
    seed = genome2[genome2_mid, genome2_mid + variables.base_seq_len]
    while len(seed) <= len(genome2) - variables.base_seq_len:
        found, genome1_mid = ht.get(seed)
        if found:
            best_cunt_found, genome1_cut, genome2_cut = seed_n_extend(genome1, genome2, genome1_mid, genome2_mid)
            if best_cunt_found:
                return True, genome1_cut, genome2_cut
            else:
                seed, genome2_mid = reboot_best_cut(genome2, genome2_mid)
                continue
        else:
            seed, genome2_mid = reboot_best_cut(genome2, genome2_mid)
            continue
    return False, genome1_cut, genome2_cut

def divide_n_conquer_seq_align():
    pass

,
def main():
    genome1, genome2 = create_genome_arrays()
    ht = create_hash_table(genome1)
    find_best_cut(ht, genome1, genome2)


if __name__ == '__main__':
    main()
