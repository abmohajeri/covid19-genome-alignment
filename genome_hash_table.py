import variables
from genome_linked_list import LinkedList


class HashTable:
    def __init__(self, size):
        self.table = [None] * size

    @staticmethod
    def hash_key(sequence) -> int:
        chars = list(sequence)
        binary_str = ''
        for char in chars:
            binary_str += variables.codes[char]
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
            return self.table[idx].find(genome1_len / 2)
        else:
            return False, -1

    def print_all(self):
        for i in range(len(self.table)):
            if self.table[i]:
                self.table[i].print_all(i)
