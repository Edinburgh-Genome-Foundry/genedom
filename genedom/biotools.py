import os
import re
from copy import deepcopy
import numpy as np

from snapgene_reader import snapgene_file_to_seqrecord

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO


complements_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}


def random_dna_sequence(length, probas=None, seed=None):
    """Return a random DNA sequence ("ATGGCGT...") with the specified length.

    Parameters
    ----------

    length
      Length of the DNA sequence.

    proba
      Frequencies for the different nucleotides, for instance
      ``probas={"A":0.2, "T":0.3, "G":0.3, "C":0.2}``.
      If not specified, all nucleotides are equiprobable (p=0.25).

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility

    """
    if seed is not None:
        np.random.seed(seed)
    if probas is None:
        sequence = np.random.choice(list("ATCG"), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p=probas)
    return "".join(sequence)


formats_dict = {
    '.fa': 'fasta',
    '.gb': 'genbank',
    '.gbk': 'genbank',
    '.dna': 'snapgene'
}


def load_record(filename, linear=True, name="unnamed", capitalize=True):
    no_extension, extension = os.path.splitext(filename)
    fmt = formats_dict[extension]
    if fmt == 'snapgene':
        record = snapgene_file_to_seqrecord(filename)
    else:
        record = SeqIO.read(filename, fmt)
    if capitalize:
        record.seq = record.seq.upper()
    record.linear = linear
    record.id = name
    record.name = name.replace(" ", "_")[:20]

    return record


def load_records(path, capitalize=True):
    if isinstance(path, (list, tuple)):
        return [record for p in path for record in load_records(p)]
    no_extension, extension = os.path.splitext(path)
    fmt = formats_dict[extension]
    if fmt == 'snapgene':
        records = [snapgene_file_to_seqrecord(path)]
    else:
        records = list(SeqIO.parse(path, fmt))
    for i, record in enumerate(records):
        if capitalize:
            record.seq = record.seq.upper()
        if str(record.id) in ['None', '', "<unknown id>", '.', ' ']:
            record.id = path.replace("/", "_").replace("\\", "_")
            if len(records) > 1:
                record.id += "_%04d" % i
    return records


def complement(sequence):
    return "".join(complements_dict[c] for c in sequence)

def reverse_complement(sequence):
    return complement(sequence)[::-1]


def sequence_to_record(sequence, features=()):
    return SeqRecord(Seq(sequence, alphabet=DNAAlphabet()),
                     features=list(features))


def annotate_record(seqrecord, location="full", feature_type="feature",
                    margin=0, **qualifiers):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
    """

    if location == "full":
        location = (margin, len(seqrecord) - margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )


def sanitize_string(string, max_length=15,
                    replacements=(("'", "p"), ("*", "s"), ("-", "_"))):
    for old, new in replacements:
        string = string.replace(old, new)
    string = re.sub(r'[^a-zA-Z\d\S]', '_', string)
    return string[:max_length]


def sanitize_and_uniquify(strings, max_length=15,
                          replacements=(("'", "p"), ("*", "s"), ("-", "_"))):
    dejavu = set()
    table = {}
    for string in strings:
        newstring = sanitize_string(string, max_length=max_length,
                                    replacements=replacements)
        i = 1
        while newstring in dejavu:
            i += 1
            newstring = newstring[:-1] + str(i)
        dejavu.add(newstring)
        table[string] = newstring
    return table

def write_record(record, target, fmt='genbank'):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes"""
    record = deepcopy(record)
    record.name = record.name[:20]
    if str(record.seq.alphabet.__class__.__name__) != 'DNAAlphabet':
        record.seq.alphabet = DNAAlphabet()
    if hasattr(target, 'open'):
        target = target.open('w')
    SeqIO.write(record, target, fmt)