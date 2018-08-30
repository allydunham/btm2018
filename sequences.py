#!/usr/bin/env python3
"""
Class creation exercise based on biological sequences
Implements a generic sequence class and specific subclasses for DNA/RNA/proteins
"""

class NucleotideError(Exception):
    """Error for a sequence containing a character outside the allowed alphabet"""
    pass

class RNA:
    """An RNA sequence"""
    alphabet = ('a', 'c', 'g', 'u')

    # Molecular weights from thermofisher
    base_weights = {'a':329.2, 'c':305.2, 'g':345.2, 'u':306.2}

    def __init__(self, seq):
        # Make sequence lower case
        seq = str(seq).lower()


        # Check all bases are part of the sequence alphabet
        for base in seq:
            if not base in self.alphabet:
                raise NucleotideError('{} is not in the alphabet {}'.format(base, self.alphabet))

        # Set the sequence
        self.seq = seq

    def __str__(self):
        return "5'-{}-3'".format(self.seq)

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, self.seq)

    def __len__(self):
        return len(self.seq)

    def molecular_weight(self):
        """Calculate the (approximate) molecular weight of the molecule, excluding ptential end phosphate"""
        weight = 0
        for base in self.seq:
            weight += self.base_weights[base]

        return weight

class ssDNA(RNA):
    """A ssDNA sequence"""
    alphabet = ('a', 'c', 'g', 't')

    # Molecular weights from thermofisher
    base_weights = {'a':313.2, 'c':289.2, 'g':329.2, 't':304.2}

    def __init__(self, seq):
        super().__init__(seq)

class dsDNA(ssDNA):
    """A dsDNA sequence"""
    base_complements = {'a':'t', 'c':'g', 'g':'c', 't':'a'}
    def __init__(self, seq):
        super().__init__(seq)
        self.complement = ''.join([self.base_complements[base] for base in self.seq])

    def __str__(self):
        foreward = super().__str__()
        reverse = "3'-{}-5'".format(self.complement)
        return '{}\n   {}   \n{}'.format(foreward, '|'*len(self), reverse)

    def molecular_weight(self):
        foreward_strand = super().molecular_weight()
        reverse_strand = 0
        for base in self.complement:
            reverse_strand += self.base_weights[base]

        return foreward_strand + reverse_strand

if __name__ == '__main__':
    try:
        bad_rna = RNA('acgacgdrcggcuacgu')
    except NucleotideError as error:
        print(error)

    try:
        bad_ssdna = ssDNA('sdfacgtgtcgatcgt')
    except NucleotideError as error:
        print(error)

    rna = RNA('acgacgucggcuacgu')
    ssdna = ssDNA('acgtgtcgatcgt')
    dsdna = dsDNA('acgtgtcgatcgt')