from copy import deepcopy

from Bio import SeqIO
import pandas
from proglog import MuteProgressBarLogger, TqdmProgressBarLogger

import flametree
from sequenticon import sequenticon

from .PartDomesticator import PartDomesticator
from .reports import domestication_report
from .biotools import sanitize_and_uniquify

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
    domesticated_records = []
    columns = ["Record", "Ordering Name", "Domesticator",
               "Domesticated Record", "Added bp", "Edited bp"]
    for record in logger.iter_bar(record=records):
        record = deepcopy(record)
        original_id = record.id
        domesticated_id = record.id + domesticate_suffix
        domesticated_file_name = domesticated_id + ".gb"
        if isinstance(domesticator, PartDomesticator):
            record_domesticator = domesticator
        else:
            record_domesticator = domesticator(record)
        domesticators.add(record_domesticator)
        if include_optimization_reports:
            report_target = errors_dir._dir(record.id)
        else:
            report_target = None
        final, edits, report, success, msg = record_domesticator.domesticate(
            record, report_target=report_target, edit=allow_edits)
        if not success:
            nfails += 1
        final.original_id = original_id
        final.id = domesticated_id.replace(' ', '_')
        SeqIO.write(final, domesticated_dir._file(domesticated_file_name),
                    "genbank")
        domesticated_records.append(final)
        if include_original_records:
            record.id = record.id[:20]
            SeqIO.write(final, original_dir._file(original_id + ".gb"),
                        "genbank")
        n_edits = sum([len(f) for f in edits.features
                       if f.qualifiers.get("is_edit", False)])
        added_bp = len(final) - len(record)
        before_seqicon = sequenticon(record, output_format="html_image")
        after_seqicon = sequenticon(final, output_format="html_image")

        infos.append({
            "id": original_id,
            "Record": before_seqicon + original_id,
            "Domesticator": record_domesticator.name,
            "Domesticated Record": ("Failed: " + msg) if not success
                                   else (after_seqicon + domesticated_id),
            "Added bp": added_bp,
            "Edited bp": n_edits
        })
    sanitizing_table = sanitize_and_uniquify([info['id'] for info in infos])
    order_id_dataframe = pandas.DataFrame(list(sanitizing_table.items()),
                                          columns=["sequence", "order_id"])
    order_id_dataframe.to_csv(root._file("order_ids.csv").open('w'),
                              index=False)
    for info in infos:
        info['Order ID'] = sanitizing_table[info['id']]
    columns = ["Record", "Order ID", "Domesticator",
               "Domesticated Record", "Added bp", "Edited bp"]
    infos_dataframe = pandas.DataFrame(infos, columns=columns)
    infos_dataframe.sort_values('Order ID', inplace=True)
    domesticators = sorted(domesticators, key=lambda d: d.name)
    domestication_report(root._file("Report.pdf"), infos_dataframe,
                         domesticators)
    order_dir = root._dir('sequences_to_order', replace=True)
    for r in domesticated_records:
        r.id = sanitizing_table[r.original_id]
        r.name = ''
        r.description = ''
    SeqIO.write(domesticated_records, order_dir._file("sequences_to_order.fa"),
                "fasta")
    df = pandas.DataFrame.from_records(
        sorted([
            {'sequence': str(rec.seq).upper(),
             'length': len(rec),
             'sequence name': rec.id}
            for rec in domesticated_records
        ], key=lambda d: d['sequence name']),
        columns=['sequence name', 'length', 'sequence']
    )
    df.to_excel(order_dir._file("sequences_to_order.xls").open("wb"),
                index=False)
    df.to_csv(order_dir._file("all_domesticated_parts.csv").open("w"),
              index=False)
    return nfails, root._close()
