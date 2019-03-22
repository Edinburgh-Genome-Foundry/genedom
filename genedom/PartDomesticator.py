"""Defines central class BlockFinder."""
from collections import OrderedDict
import pandas
import numpy as np

from Bio import SeqIO, Restriction

from dnachisel import (AvoidPattern, reverse_complement, reverse_translate,
                       DnaOptimizationProblem, EnforceTranslation,
                       annotate_record, CodonOptimize, Location,
                       sequence_to_biopython_record, random_dna_sequence,
                       AvoidChanges, AvoidChanges, EnzymeSitePattern)

from dnachisel.reports import SpecAnnotationsTranslator

from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
import matplotlib.cm as cm


from .biotools import reverse_complement, sequence_to_record, annotate_record
from .StandardDomesticatorsSet import StandardDomesticatorsSet
from .DomesticationResult import DomesticationResult

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
                    final_record_target=None, edit=False, barcode='',
                    barcode_spacer='AA', report_target=None):
        """Domesticate a sequence.

        Parameters
        ----------

        dna_sequence
          The DNA sequence string to domesticate.
        
        protein_sequence
          Amino-acid sequence of the protein, which will be converted into
          a DNA sequence string.
        
        is_cds
          If True, sequence edits are restricted to synonymous mutations.

        codon_optimization
          Either None for no codon optimization or the name of an organism
          supported by DnaChisel.
        
        extra_constraints
          List of extra constraints to apply to the domesticated sequences.
          Each constraint is either a DnaChisel constraint or a function
          (dna_sequence => DnaChisel constraint).
        
        extra_objectives
          List of extra optimization objectives to apply to the domesticated
          sequences. Each objective is either a DnaChisel constraint or a
          function (dna_sequence => DnaChisel constraint).
        
        final_record_target
          Path to the file where to write the final genbank
        
        edit
          Turn to True to allow sequence edits (if it is false and no all
          constraints are originally satisfied, a failed domestication result
          (i.e. with attribute ``success`` set to False) will be returned.
        
        report_target
          Target for the sequence optimization report (a folder path, or a zip
          path)
        
        barcode
          A sequence of DNA that will be added to the left of the sequence once
          the domestication is done.

        barcode_spacer
          Nucleotides to be added between the barcode and the enzyme (optional,
          the idea here is that they will make sure to avoid the creation of
          unwanted cutting sites).

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
            (success, message, report_data) = problem.optimize_with_report(
                target=report_target,
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

        return DomesticationResult(problem.sequence_before, final_record,
                                   edits_record, report_data,
                                   optimization_successful, message)


    def details_list(self):
        """List of details for representing the domesticator in reports."""
        return [(label, value) for (label, value) in [
            ("Name", self.name),
            ("Description", self.description),
            ("Left addition", str(self.left_flank.seq)),
            ("Right addition", str(self.right_flank.seq)),
        ] if value not in (None, "")]

    def html_details(self):
        """HTML representation of the ``details_list``, for reports."""
        return "<br />".join([
            "<b>%s</b>: %s" % (name, value)
            for (name, value) in self.details_list()
        ])

    @staticmethod
    def plot_record(record, ax=None):
        """Plot the given record with a custom DnaFeaturesViewer plotter."""
        translator = SpecAnnotationsTranslator()
        gr_record = translator.translate_record(record)
        return gr_record.plot(ax=ax)


class GoldenGateDomesticator(PartDomesticator):
    """Special domesticator class for Golden-Gate standards
    
    Parameters
    ----------

    left_overhang
      4bp overhang to be added on the left
    
    right_overhang
      4bp overhang to be added on the right
    
    left_addition
      Extra sequence of DNA to be systematically added on the left of each part
      between the enzyme site and the rest of the sequence.
    
    right_addition
      Extra sequence to be systematically added on the right of each part
      between the enzyme site and the rest of the sequence.
    
    enzyme
      Enzyme used for the Golden Gate assembly. This enzyme will be added on
      the flanks of the sequence, and the internal sequence will be protected
      against sites from this enzyme during optimization.
    
    extra_avoided_sites
      Other enzymes from which the sequence should be protected during
      optimization in addition to the assembly ``enzyme``.
    
    description
      Description of the domesticator as it will appear in reports.
    
    name
      Name of the domesticator as it will appear in reports

    constraints
      Either Dnachisel constraints or functions (sequence => constraint) to be
      applied to the sequence for optimization

    objectives
      Either Dnachisel objectives or functions (sequence => objective) to be
      applied to the sequence for optimization.
    """

    def __init__(self, left_overhang, right_overhang, left_addition='',
                 right_addition='', enzyme='BsmBI', extra_avoided_sites=(),
                 description='Golden Gate domesticator',
                 name='unnamed_domesticator', constraints=(), objectives=()):
        self.enzyme = enzyme
        self.left_overhang = left_overhang
        left_overhang = sequence_to_biopython_record(left_overhang)
        self.right_overhang = right_overhang
        right_overhang = sequence_to_biopython_record(right_overhang)
        for seq in [left_overhang, right_overhang]:
            annotate_record(seq, label=str(seq.seq))
        enzyme_seq = Restriction.__dict__[enzyme].site
        enzyme_seq = sequence_to_biopython_record(enzyme_seq)
        annotate_record(enzyme_seq, label=enzyme)
        self.enzyme_seq = enzyme_seq
        left_flank = self.enzyme_seq + "A" + left_overhang + left_addition
        right_flank = (right_addition + right_overhang +
                       (self.enzyme_seq + "A").reverse_complement())
        constraints = list(constraints) + [
            (lambda seq: AvoidPattern(
                EnzymeSitePattern(enzyme),
                location=Location(len(left_flank),
                                  len(left_flank) + len(seq))))
            for enz in ([enzyme] + list(extra_avoided_sites))
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
    def standard_from_spreadsheet(path=None, dataframe=None, name_prefix=''):
        """Parse a spreadsheet into a standard with Golden Gate domesticators.

        The input should be a table with columns names as follows:
        slot_name, left_overhang, right_overhang, left_addition,
        right_addition, enzyme, extra_avoided_sites, description.

        Parameters
        ----------

        path
          Path to a CSV or XLS(X) file. A dataframe can be provided instead.
        
        dataframe
          A pandas Dataframe which can be provided instead of a path
          
        """
        if path is not None:
            if path.lower().endswith(".csv"):
                dataframe = pandas.read_csv(path)
            else:
                dataframe = pandas.read_excel(path)
        return StandardDomesticatorsSet(OrderedDict([
            (row.slot_name, GoldenGateDomesticator(
                left_overhang=row.left_overhang,
                right_overhang=row.right_overhang,
                left_addition=nan_to_empty_string(row.left_addition),
                right_addition=nan_to_empty_string(row.right_addition),
                enzyme=row.enzyme,
                extra_avoided_sites = [
                  e.strip() for e in row.extra_avoided_sites.split(',')
                ] if hasattr(row.extra_avoided_sites, 'split') else [],
                description=row.description,
                name=name_prefix + row.slot_name
            ))
            for i, row in dataframe.iterrows()
        ]))
