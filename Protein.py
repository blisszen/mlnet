# -*- coding: ms949 -*-

'''
Protein�� ���õ� ������ �۾����̴�.

'''

SingleLetters = [ 'A', 'C', 'D', 'E', 'T', 
                  'F', 'G', 'H', 'I', 'V', 
                  'K', 'L', 'M', 'N', 'W',
                  'P', 'Q', 'R', 'S', 'Y' ]

ThreeLetters = { 'A':'Ala', 'C':'Cys', 'D':'Asp', 'E':'Glu', 'T':'Thr',
                 'F':'Phe', 'G':'Gly', 'H':'His', 'I':'Ile', 'V':'Val',
                 'K':'Lys', 'L':'Leu', 'M':'Met', 'N':'Asn', 'W':'Trp',
                 'P':'Pro', 'Q':'Gln', 'R':'Arg', 'S':'Ser', 'Y':'Tyr'}


def cleaning(sequence):
	'''Single letter amino acid�� �ƴ� �͵��� �� ������'''
	seq = sequence.upper()

	r = ''
	for s in seq:
		if s in SingleLetters:
			r += s
	return r 

def convertSingleLetter2ThreeLetter(sequence):
	seq = cleaning(sequence)
	r = ''
	for s in seq:
		r += ThreeLetters[s]
	return r


def hasUnidentifiedResidues(seq):
	for s in seq:
		if not s in SingleLetters:
			#print '[Unidentified residue] ', s
			return True
	return False