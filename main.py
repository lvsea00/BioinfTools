from abc import ABC, abstractmethod

class BiologicalSequence(ABC):
    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)
    
    def __getitem__(self, index):
        return self.seq[index]

    def __str__(self):
        return f"{self.__class__.__name__}: {self.seq}"
    
    @abstractmethod
    def validate(self):
        pass


