import variables


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
