"""Defines central class BlockFinder."""
from collections import OrderedDict
import pandas
import numpy as np

from Bio import Restriction

from dnachisel import (
    AvoidPattern,
    annotate_record,
    Location,
    sequence_to_biopython_record,
    EnzymeSitePattern,
)

from ..StandardDomesticatorsSet import StandardDomesticatorsSet
from .PartDomesticator import PartDomesticator


def nan_to_empty_string(val):
    """Return the value unless it is NaN, then it returns an empty string."""
    return val if (isinstance(val, str) or not np.isnan(val)) else ""


class GoldenGateDomesticator(PartDomesticator):
    """Special domesticator class for Golden-Gate standards.

    Parameters
    ----------

    left_overhang
      4bp overhang to be added on the left.

    right_overhang
      4bp overhang to be added on the right.

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
      Name of the domesticator as it will appear in reports.

    constraints
      Either Dnachisel constraints or functions (sequence => constraint) to be
      applied to the sequence for optimization.

    objectives
      Either Dnachisel objectives or functions (sequence => objective) to be
      applied to the sequence for optimization.
    """

    def __init__(
        self,
        left_overhang,
        right_overhang,
        left_addition="",
        right_addition="",
        enzyme="BsmBI",
        extra_avoided_sites=(),
        description="Golden Gate domesticator",
        name="unnamed_domesticator",
        cds_by_default=False,
        constraints=(),
        objectives=(),
    ):
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
        right_flank = (
            right_addition
            + right_overhang
            + (self.enzyme_seq + "A").reverse_complement()
        )
        self.extra_avoided_sites = extra_avoided_sites
        constraints = list(constraints) + [
            (
                lambda seq: AvoidPattern(
                    EnzymeSitePattern(enzyme),
                    location=Location(len(left_flank), len(left_flank) + len(seq)),
                )
            )
            for enz in ([enzyme] + list(extra_avoided_sites))
        ]
        PartDomesticator.__init__(
            self,
            left_flank=left_flank,
            right_flank=right_flank,
            constraints=constraints,
            objectives=objectives,
            description=description,
            name=name,
            cds_by_default=cds_by_default,
        )

    def __repr__(self):
        return "GgDomesticator[%s](%s-%s)" % (
            self.enzyme,
            self.left_overhang,
            self.right_overhang,
        )

    def __str__(self):
        return "GgDomesticator[%s](%s-%s)" % (
            self.enzyme,
            self.left_overhang,
            self.right_overhang,
        )

    def details_list(self):
        result = PartDomesticator.details_list(self) + [
            ("Enzyme", "%s (%s)" % (self.enzyme, str(self.enzyme_seq.seq))),
            ("Left overhang", self.left_overhang),
            ("Right overhang", self.right_overhang),
        ]
        extra_sites = self.extra_avoided_sites
        if len(extra_sites):
            result += [("Other avoided sites", ", ".join(extra_sites))]
        return result

    @staticmethod
    def standard_from_spreadsheet(path=None, dataframe=None, name_prefix=""):
        """Parse a spreadsheet into a standard with Golden Gate domesticators.

        The input should be a table with columns names as follows:
        slot_name, left_overhang, right_overhang, left_addition,
        right_addition, enzyme, extra_avoided_sites, description.

        Parameters
        ----------

        path
          Path to a CSV or XLS(X) file. A dataframe can be provided instead.

        dataframe
          A pandas Dataframe which can be provided instead of a path.
        """
        if path is not None:
            if path.lower().endswith(".csv"):
                dataframe = pandas.read_csv(path)
            else:
                dataframe = pandas.read_excel(path)
        return StandardDomesticatorsSet(
            OrderedDict(
                [
                    (
                        row.slot_name,
                        GoldenGateDomesticator(
                            left_overhang=row.left_overhang,
                            right_overhang=row.right_overhang,
                            left_addition=nan_to_empty_string(row.left_addition),
                            right_addition=nan_to_empty_string(row.right_addition),
                            enzyme=row.enzyme,
                            extra_avoided_sites=[
                                e.strip() for e in row.extra_avoided_sites.split(",")
                            ]
                            if hasattr(row.extra_avoided_sites, "split")
                            else [],
                            description=row.description,
                            cds_by_default=(row.is_cds == "yes")
                            if hasattr(row, "is_cds")
                            else False,
                            name=name_prefix + row.slot_name,
                        ),
                    )
                    for i, row in dataframe.iterrows()
                ]
            )
        )
