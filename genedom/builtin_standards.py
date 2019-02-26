import os
from box import Box
from .PartDomesticator import GoldenGateDomesticator
this_dir = os.path.realpath(__file__)
standards_dir = os.path.join(os.path.dirname(this_dir), "assembly_standards")

BUILTIN_STANDARDS = Box({})
for fname in os.listdir(standards_dir):
    name, ext = os.path.splitext(fname)
    if ext != ".csv":
        continue
    path = os.path.join(standards_dir, fname)
    BUILTIN_STANDARDS[name] = GoldenGateDomesticator.standard_from_spreadsheet(
        path, name_prefix=name + "_")
