"""Defines central class BlockFinder."""
from collections import OrderedDict
import pandas
import numpy as np

from Bio import SeqIO, Restriction

from dnachisel import (AvoidPattern, reverse_complement, reverse_translate,
                       DnaOptimizationProblem, EnforceTranslation,
                       annotate_record, CodonOptimize, Location,
                       sequence_to_biopython_record, random_dna_sequence,
                       AvoidChanges, AvoidChanges)

from dnachisel.reports import (optimization_with_report,
                               SpecAnnotationsTranslator)


from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
import matplotlib.cm as cm


from .biotools import reverse_complement, sequence_to_record, annotate_record


def nan_to_empty_string(val):
    """Return the value unless it is NaN, then it returns an empty string."""
    return val if (isinstance(val, str) or not np.isnan(val)) else ''


class PartDomesticator:
    """Generic domesticator.

    Parameters
    ----------

    name
      Domesticator name as it will appear in reports etc.
    
    description
      Short domesticator description as it will appear in reports etc.
    
    left_flank
      String. Left addition to the sequence (homology arms, enzymes sites etc.)
      

    right_flank
      String. Right addition to the sequence (homology arms, enz. sites etc.)

    constraints
      Either Dnachisel constraints or functions (sequence => constraint) to be
      applied to the sequence for optimization


    objectives
      Either Dnachisel objectives or functions (sequence => objective) to be
      applied to the sequence for optimization.

    simultaneous_mutations
      Number of sequences mutations to be applied simulatenously during
      optimization. A larger number creates more noise but could allow to
      solve tougher problems.

    minimize_edits
      If true, the optimizer will attempt to minimize changes while making
      sure the constraints hold (each edit incurs a penalty of 1 in the
      total optimization score).

    logger
      A proglog logger or 'bar' or None for no logger at all.
    """

    def __init__(self, name='unnamed domesticator', left_flank='',
                 right_flank='', constraints=(), objectives=(),
                 description=None, simultaneous_mutations=1,
                 minimize_edits=True, logger=None):
        if isinstance(left_flank, str):
            left_flank = sequence_to_biopython_record(left_flank)
            annotate_record(left_flank, label='left flank')
        if isinstance(right_flank, str):
            right_flank = sequence_to_biopython_record(right_flank)
            annotate_record(right_flank, label='right flank')
        self.name = name
        self.constraints = constraints
        self.left_flank = left_flank
        self.right_flank = right_flank
        self.constraints = list(constraints)
        self.objectives = list(objectives)
        self.description = description
        self.logger = logger
        self.simultaneous_mutations = simultaneous_mutations
        self.minimize_edits = minimize_edits

    def domesticate(self, dna_sequence=None, protein_sequence=None,
                    is_cds=False, codon_optimization=None,
                    extra_constraints=(), extra_objectives=(),
                    final_record_target=None, edit=False,
                    report_target=None):
        """Domesticate a sequence.

        Parameters
        ----------

        dna_sequence
        
        protein_sequence
        
        is_cds


        codon_optimization
        
        extra_constraints
        
        extra_objectives
        
        final_record_target
        
        edit
        
        report_target
          Target for the sequence optimization report (a folder path, or a zip
          path)

        Returns
        -------

        final_record, edits_record, report_data, success, msg
        """
        if protein_sequence is not None:
            is_cds = True
            dna_sequence = reverse_translate(protein_sequence)
        constraints = [
            c(dna_sequence) if hasattr(c, '__call__') else c
            for c in list(extra_constraints) + self.constraints
        ]
        location = Location(len(self.left_flank),
                            len(self.left_flank) + len(dna_sequence))
        if is_cds:
            constraints.append(EnforceTranslation(location=location))
        objectives = [
            o(dna_sequence) if hasattr(o, '__call__') else o
            for o in list(extra_objectives) + self.objectives
        ]
        if codon_optimization:
            objectives.append(CodonOptimize(species=codon_optimization,
                                            location=location))
        if self.minimize_edits:
            objectives.append(AvoidChanges())

        extended_sequence = self.left_flank + dna_sequence + self.right_flank

        if (not is_cds) and (not edit):
            constraints.append(AvoidChanges())
        problem = DnaOptimizationProblem(
            extended_sequence,
            constraints=constraints,
            objectives=objectives,
            logger=self.logger
        )
        problem.n_mutations = self.simultaneous_mutations
        optimization_successful = True
        message = ""
        if report_target is not None:
            (success, message, report_data) = optimization_with_report(
                target=report_target,
                problem=problem,
                project_name=self.name
            )
            optimization_successful = success
        else:
            report_data = None
            try:
                problem.resolve_constraints()
                problem.optimize()
            except Exception as err:
                message = str(err)
                optimization_successful = False
                report_data = None
        final_record = problem.to_record(
            with_original_features=True,
            with_original_spec_features=False,
            with_constraints=False,
            with_objectives=False
        )
        edits_record = problem.to_record(
            with_original_features=True,
            with_original_spec_features=False,
            with_constraints=False,
            with_objectives=False,
            with_sequence_edits=True
        )
        if final_record_target is not None:
            SeqIO.write(final_record, final_record_target, 'genbank')

        return (final_record, edits_record, report_data,
                optimization_successful, message)


    def details_list(self):
        return [(label, value) for (label, value) in [
            ("Name", self.name),
            ("Description", self.description),
            ("Left addition", str(self.left_flank.seq)),
            ("Right addition", str(self.right_flank.seq)),
        ] if value not in (None, "")]

    def html_details(self):
        return "<br />".join([
            "<b>%s</b>: %s" % (name, value)
            for (name, value) in self.details_list()
        ])

    @staticmethod
    def plot_record(record, ax=None):
        translator = SpecAnnotationsTranslator()
        gr_record = translator.translate_record(record)
        return gr_record.plot(ax=ax)


class GoldenGateDomesticator(PartDomesticator):

    def __init__(self, left_overhang, right_overhang, left_addition='',
                 right_addition='', enzyme='BsmBI',
                 description='Golden Gate domesticator',
                 name='unnamed_domesticator', constraints=(), objectives=()):
        self.enzyme = enzyme
        self.left_overhang = left_overhang
        self.right_overhang = right_overhang
        self.enzyme_seq = Restriction.__dict__[enzyme].site
        left_flank = self.enzyme_seq + "A" + left_overhang + left_addition
        right_flank = (right_addition + right_overhang +
                       reverse_complement(self.enzyme_seq + "A"))
        constraints = list(constraints) + [
            lambda seq: AvoidPattern(
                enzyme=enzyme,
                location=Location(len(left_flank), len(left_flank) + len(seq)))
        ]
        PartDomesticator.__init__(self, left_flank=left_flank,
                                  right_flank=right_flank,
                                  constraints=constraints,
                                  objectives=objectives,
                                  description=description,
                                  name=name)

    def __repr__(self):
        return "GgDomesticator[%s](%s-%s)" % (self.enzyme, self.left_overhang,
                                              self.right_overhang)

    def __str__(self):
        return "GgDomesticator[%s](%s-%s)" % (self.enzyme, self.left_overhang,
                                              self.right_overhang)
    def details_list(self):
        return PartDomesticator.details_list(self) + [
            ("Enzyme", "%s (%s)" % (self.enzyme, self.enzyme_seq)),
            ("Left overhang", self.left_overhang),
            ("Right overhang", self.right_overhang)
        ]

    @staticmethod
    def from_spreadsheet(path=None, dataframe=None, name_prefix=''):
        if path is not None:
            if path.lower().endswith("csv"):
                dataframe = pandas.read_csv(path)
            else:
                dataframe = pandas.read_excel(path)

        return OrderedDict([
            (row.slot_name, GoldenGateDomesticator(
                left_overhang=row.left_overhang,
                right_overhang=row.right_overhang,
                left_addition=nan_to_empty_string(row.left_addition),
                right_addition=nan_to_empty_string(row.right_addition),
                enzyme=row.enzyme,
                description=row.description,
                name=name_prefix + row.slot_name
            ))
            for i, row in dataframe.iterrows()
        ])
