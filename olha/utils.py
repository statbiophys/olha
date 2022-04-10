#!/usr/bin/env python3

import itertools
import os
import re
import numpy as np

bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
aminoacids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
to_aminoacid = dict(zip(codons, aminoacids))
to_codons = {}
for a, c in zip(aminoacids, codons):
    if a in to_codons:
        to_codons[a].append(c)
    else:
        to_codons[a] = [c]

def n_to_aa(x):
    """ Transform the sequence x to its nucleotide form """
    return "".join([to_aminoacid[x[i:i+3]] for i in range(0, len(x), 3)])


def aa_to_n(x):
    """ Return all the possible versions of the aa sequence x
    in nucleotide form """
    for s in itertools.product(*(to_codons[a] for a in x)):
        yield "".join(s)


def gene_map(x, genomic_data):
    """ Map relatively standard gene name to
        the convention used in Alice.
        In case there's a doubt returns more than needed.
        It's far from perfect.
        @ Arguments:
        * x: V or J gene name, form: TYPE##-##*##, or TYPE##-##
        or TYPE[V,J]##, where ## can be interpreted as digits/letters
        and TYPE = "IGH","IGK","IGL" or "TRB"/"TRA"/"TRG"/"TRD"
        * genomic_data: the olga model
        @ Return:
        * list of str, V-gene names that Alice can understand
        @ Example:
        "IGHV07" -> ["IGHV7-34-1*01", "IGHV7-34-1*02", "IGHV7-4-1*01",
                     "IGHV7-4-1*02", "IGHV7-4-1*03","IGHV7-4-1*04",
                     "IGHV7-4-1*05", "IGHV7-40*01", "IGHV7-81*01"]
        "TRBV07-002*4 -> ["TRBV7-2*04"]
    """

    regex = (r"^(TRB|TRA|IGH|IGK|IGL|TRG|TRD)(?:\w+)?(V|D|J)"
             r"([\w-]+)?(?:/DV\d+)?(?:\*(\d+))?(?:/OR.*)?$")
    g = re.search(regex, x)

    chain = None
    gene_type = None
    gene_id = None
    allele = None

    if g is None:
        raise ValueError("Gene {} does not have a valid name".format(x))
    chain = g.group(1)
    gene_type = g.group(2)
    gene_id = g.group(3)
    allele = None if g.group(4) is None else int(g.group(4))

    if chain is None or gene_type is None:
        raise ValueError("Gene {} does not have a valid name".format(x))

    # check if gene_id contain something of the form
    # ##-## where ## is a digit or ##S##
    gene_id_1 = None
    gene_id_2 = None
    if gene_id is not None:
        g = re.search(r'(\d+)(?:[-S](\d+))?', gene_id)
        if g is not None:
            if g.span()[1] >= 3:
                gene_id_1 = int(g.group(1))
                gene_id_2 = int(g.group(2))
            else:
                gene_id_1 = int(g.group(1))

    key = chain + gene_type
    possible_genes = igor_genes(key, genomic_data)


    if allele is not None and gene_id_1 is not None and gene_id_2 is not None:
        guess = [a[0] for a in possible_genes if a[1] == gene_id_1
                 and a[2] == gene_id_2 and a[4] == allele]
        if guess != []:
            return guess
    if gene_id_1 is not None and gene_id_2 is not None:
        guess = [a[0] for a in possible_genes if a[1] == gene_id_1
                 and a[2] == gene_id_2]
        if guess != []:
            return guess
    if allele is not None and gene_id_1 is not None:
        guess = [a[0] for a in possible_genes if a[1] == gene_id_1
                 and a[4] == allele]
        if guess != []:
            return guess
    if gene_id_1 is not None:
        guess = [a[0] for a in possible_genes if a[1] == gene_id_1]
        if guess != []:
            return guess
    # if everything else failed return all the genes we have
    return [a[0] for a in possible_genes]


def igor_genes(key, genomic_data):
    """ Read the model directory and return all the gene matching key.
        If model and model_dir are None only consider human genes.

        It returns the full gene name, plus the gene family, its name, its allele
    """

    regex = (r"(\d+)(?:P)?(?:[\-S](\d+)(?:D)?(?:\-(\d+))?)?"
             r"(?:/DV\d+)?(?:-NL1)?(?:\*(\d+))?")  # match all the IGoR gene names


    lst = [tuple([gene[0]] + [None if a is None else int(a) for a in
                           re.search(key + regex, gene[0]).groups()])
           for gene in genomic_data.genV + genomic_data.genJ
           if re.search(key + regex, gene[0]) is not None]
    return lst


def igor_genes_to_idx(genomic_data, key):
    """ Return a dictionary that maps igor defined genes to their indices
        @ Arguments:
        * model: the place where model_params is located
        * key: the type of sequence "V" or "J"
    """
    dct = {}
    if key == "V":
        dct = {v[0]: ii for ii, v in enumerate(genomic_data.genV)}
    if key == "J":
        dct = {j[0]: ii for ii, j in enumerate(genomic_data.genJ)}
    return dct
