import os
from dnachisel import (AllowPrimer, DnaOptimizationProblem, EnzymeSitePattern,
                       EnforceSequence, EnforceGCContent, AvoidPattern,
                       random_dna_sequence)
from collections import OrderedDict
from .biotools import sequence_to_record, annotate_record, write_record

class BarcodesCollection(OrderedDict):
    """Class representing a set of named barcode sequences.

    These barcodes are meant to be annealed with same-sequence primers for PCR
    or sequencing.
    
    The constructor taked a list [(name, barcode), ...] as an input.

    Use ``BarcodesCollection.from_specs(n_barcodes=25)`` to generate an
    instance with 25 compatible barcodes.
    """

    def __init__(self, barcodes):
        OrderedDict.__init__(self, barcodes)

    @staticmethod
    def from_specs(n_barcodes=96, barcode_length=20, spacer='AA',
                   forbidden_enzymes=('BsaI', 'BsmBI', 'BbsI'),
                   barcode_tmin=55, barcode_tmax=70,
                   other_primer_sequences=(), heterodim_tmax=5,
                   max_homology_length=10, include_spacers=True,
                   names_template="B_%03d"):
        """Return a BarcodesCollection object with compatible barcodes.

        Parameters
        ----------

        n_barcodes
          Number of barcodes to design
        
        barcode_length
          Length of each barcode
        
        spacer
          Spacer to place between each barcode during the optimization,
          ideally the same spacer that will be used when adding the barcode
          to a part.
        
        include_spacers
          Whether the spacers should be part of the final sequence of the
          barcodes (they still won't be considered part of the annealing
          primer and won't be used for melting temperature computations)

        forbidden_enzymes
          Name of enzymes whose sites should not be in the barcodes.

        barcode_tmin, barcode_tmax
          Interval of acceptable values for the melting temperature
        
        other_primer_sequences
          External sequences with which the primers should not anneal.
        
        heterodim_tmax
          Max acceptable melting temperature for the annealing of a barcode
          and one of the other_primer_sequences.
        
        max_homology_length
          Maximal homology between any two barcodes in the sequence.
        
        names_template
          The template used to name barcode number "i".
        """
        unit_length = barcode_length + len(spacer)
        seq_len = n_barcodes * unit_length
        units_coordinates = [(i, i + unit_length)
                            for i in range(0, seq_len, unit_length)]
        constraints = [
            AvoidPattern(EnzymeSitePattern(enzyme))
            for enzyme in forbidden_enzymes
        ]
        for start, end in units_coordinates:
            constraints += [
                AllowPrimer(
                    tmin=barcode_tmin,
                    tmax=barcode_tmax,
                    max_homology_length=max_homology_length,
                    avoid_heterodim_with=None,
                    max_heterodim_tm=5,
                    location=(start, end - len(spacer))
                ),
                EnforceSequence(spacer, location=(end - len(spacer), end)),
                EnforceGCContent(mini=0.4, maxi=0.6,
                                    location=(start, end - len(spacer)))   
            ]
        problem = DnaOptimizationProblem(
            sequence= random_dna_sequence(seq_len),
            constraints=constraints
        )
        problem.logger.ignored_bars.add('location')
        problem.resolve_constraints()

        barcodes = [problem.sequence[start: end]
                    for (start, end) in units_coordinates]
        if not include_spacers:
            barcodes = [b[:-len(spacer)] for b in barcodes]
        names = [(names_template % (i + 1)) for i in range(len(barcodes))]
        return BarcodesCollection(zip(names, barcodes))
        
        
    def to_sequences_list(self):
        """Return a list of sequences ["ATTG...", "TTCTGT..."]"""
        return list(self.values())

    def to_fasta(self, path=None):
        """Return (and optionally write) a fasta string of the barcodes."""
        fasta = "\n\n".join(["> %s\n%s" % (name, barcode)
                             for name, barcode in self.items()])
        if path is not None:
            with open(path, "w+") as f:
                f.write(fasta)
        else:
            return fasta
    
    def to_records(self, path=None):
        """Return (optionally write) individual Genbanks of the barcodes."""
        records = []
        for (name, barcode) in self.items():
            record = sequence_to_record(barcode)
            record.id = name
            annotate_record(record, label=name)
            records.append(record)
        if path is not None:
            for r in records:
                write_record(r, os.path.join(path, "%s.gb" % r.id))
        return records

