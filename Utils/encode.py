# -*- coding: utf-8 -*-

"""
This module contains functions for manipulating and encoding
peptide sequences for use in machine learning.
"""

def make_aa_list():
    ''' Returns a list of amino acid codes'''

    return ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def make_aa_pairs_list():
    ''' Returns a list of all posible pairs of amino acids'''

    aal = make_aa_list()
    return [i+j for i in aal for j in aal]


def make_idx_dict(lst):
    ''' Returns a dict with lst entries as keys and indexes as values'''

    idx_dict = {ent : idx for idx, ent in enumerate(lst)}
    return idx_dict

def make_aa_binary_code_dict():
    ''' Makes a dict with amino acid letter codes as keys and 5 digit binary
    codes as values'''

    aal = make_aa_list()
    aa_dict = {aal[i]: list(format(i, '005b')) for i in range(len(aal))}
    return aa_dict

# ==================== Peptide encoders ====================

def encode_aa_counts(pep):
    ''' Encodes the amino acid sequence as counts. Returns a list of len==20'''

    pep = pep.upper()
    aa_idx_dict = make_idx_dict(make_aa_list())
    counts = [0] * len(aa_idx_dict)
    for i in pep:
        counts[aa_idx_dict[i]] += 1
    return counts

def encode_aa_pair_counts(pep):
    ''' Encodes the amino acid sequence in terms of counts of adjacent pairs of
    amino acids'''

    pep = pep.upper()
    pair_idx_dict = make_idx_dict(make_aa_pairs_list())
    counts = [0] * len(pair_idx_dict)
    for i in range(len(pep)-1):
        counts[pair_idx_dict[pep[i:i+2]]] += 1
    return counts

def encode_aa_binary(pep):
    ''' Encodes a peptide as a binary code list'''

    aabd = make_aa_binary_code_dict()
    pep = pep.upper()

    bin_pep = []
    for i in pep:
        bin_pep.append(aabd[i])
    return bin_pep
#===========================================================

def encode_peptides(pep, encode_pep):
    ''' Encodes the peptide or list of peptides with the given encoder'''

    if isinstance(pep, str):
        return encode_pep(pep)
    else:
        peps = []
        for i in pep:
            peps.append(encode_pep(i))
        return peps


