from copy import deepcopy

from Bio import SeqIO
import pandas
from proglog import MuteProgressBarLogger, TqdmProgressBarLogger

import flametree
from sequenticon import sequenticon

from .PartDomesticator import PartDomesticator
from .reports import domestication_report

def batch_domestication(records, domesticator, target, allow_edits=False,
                        domesticate_suffix="",
                        include_optimization_reports=True,
                        include_original_records=True,
                        logger="bar"):
    if logger == "bar":
        logger = TqdmProgressBarLogger(min_time_interval=0.2)
    elif logger is None:
        logger = MuteProgressBarLogger()

    root = flametree.file_tree(target, replace=True)
    domesticated_dir = root._dir("domesticated")
    if include_original_records:
        original_dir = root._dir("original")
    if include_optimization_reports:
        errors_dir = root._dir("error_reports")

    infos = []

    domesticators = set()
    nfails = 0
    for record in logger.iter_bar(record=records):
        record = deepcopy(record)
        original_id = record.id
        new_id = record.id + domesticate_suffix
        if isinstance(domesticator, PartDomesticator):
            record_domesticator = domesticator
        else:
            record_domesticator = domesticator(record)
        domesticators.add(record_domesticator)
        if include_optimization_reports:
            report_target = errors_dir._dir(new_id)
        else:
            report_target = None
        final, edits, report, success, msg = record_domesticator.domesticate(
            record, report_target=report_target, edit=allow_edits)
        if not success:
            nfails += 1
        final.id = new_id[:20]
        SeqIO.write(final, domesticated_dir._file(new_id + ".gb"), "genbank")
        if include_original_records:
            record.id = record.id[:20]
            SeqIO.write(final, original_dir._file(original_id + ".gb"),
                        "genbank")
        n_edits = sum([len(f) for f in edits.features
                       if f.qualifiers.get("is_edit", False)])
        added_bp = len(final) - len(record)
        before_seqicon = sequenticon(record, output_format="html_image")
        after_seqicon = sequenticon(final, output_format="html_image")
        columns = ["Record", "Domesticator", "Domesticated Record", "Added bp",
                   "Edited bp"]
        infos.append([
            before_seqicon + record.id,
            record_domesticator.name,
            ("Failed: " + msg) if not success else (after_seqicon + new_id),
            added_bp, n_edits])
    infos_dataframe = pandas.DataFrame(infos, columns=columns)
    domesticators = sorted(domesticators, key=lambda d: d.name)
    domestication_report(root._file("Report.pdf"), infos_dataframe,
                         domesticators)
    return nfails, root._close()
