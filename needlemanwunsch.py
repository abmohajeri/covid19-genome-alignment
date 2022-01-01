import numpy as np
import math

base_seq_length = 6
num_seq_bits = base_seq_length * 2


def nw(x, y, match=1, mismatch=-1, gap=-1):
    nx = len(x)
    ny = len(y)
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
            if x[i] == y[j]:
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
            rx.append(x[i-1])
            ry.append(y[j-1])
            i -= 1
            j -= 1
        elif p[i, j] in [3, 5, 7, 9]:
            rx.append(x[i-1])
            ry.append('-')
            i -= 1
        elif p[i, j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j-1])
            j -= 1
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    return '\n'.join([rx, ry])


x = "AGATTCGATTACAAGAGATTACAGATTACAAATT"
y = "AGATTCGCAGGATTACAAGATGATTACAATACAA"
print(nw(x, y))

class Node:
    def __init__(self, key, value, next=None):
        self.key = key
        self.value = value
        self.next = next


class LinkedList:
    def __init__(self):
        self.head = None
        self.size = 0

    def insertAtFirst(self, key, value):
        self.head = Node(key, value, self.head)
        self.size += 1

    def insertAtLast(self,key,value):
        current = self.head
        while current.next:
            current = current.next
        current.next = Node(key,value)
        self.size+=1

    def find(self, key):
        current = self.head
        while current:
            if current.key == key:
                return current.value
            current = current.next

    def printAll(self, idx):
        current = self.head
        print("index: ", idx, "string: ", current.key, "values: ",end=' ')
        while current:
            print(current.value, end=' ')
            current = current.next
        print()


class HashTable:
    def __init__(self, size):
        self.table = [None] * (2**num_seq_bits)
        self.size = size

    def hashKey(self, key) -> int:
        CODES = {
            'A': '00',
            'C': '01',
            'G': '10',
            'T': '11'
        }
        chars = list(key)
        binary_str = ''
        for char in chars:
            binary_str += CODES[char]
        return int(binary_str, 2)

    def add(self, key, value):
        idx = self.hashKey(key)
        if self.table[idx] is None:
            newLinkedList = LinkedList()
            newLinkedList.insertAtFirst(key, value)
            self.table[idx] = newLinkedList
        else:
            self.table[idx].insertAtLast(key, value)

    def get(self, key) -> any:
        idx = self.hashKey(key)
        if self.table[idx] is not None:
            return self.table[idx].find(key)

    def printAll(self):
        for i in range(len(self.table)):
            if self.table[i]:
                self.table[i].printAll(i)


x = "AGATTCGATTACAAGAGATTACAGATTACAAATTAGATTCGATTACAAGAGATTACAGATTACAAATTAGATTCGATTACAAGAGATTACAGATTACAAATTAGATTCGATTACAAGAGATTACAGATTACAAATT"
y = "AGATTCGCAGGATTACAAGATGATTACAATACAA"
ht = HashTable(math.ceil(len(x)/4))
for i in range(len(x) - 5):
    substr = x[i:i+6]
    ht.add(substr, i)
ht.printAll()