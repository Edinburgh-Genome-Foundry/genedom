from copy import deepcopy
import itertools

from Bio import SeqIO
import pandas
import proglog

import flametree
from sequenticon import sequenticon

from .PartDomesticator import PartDomesticator
from .reports import write_pdf_domestication_report
from .biotools import (sanitize_and_uniquify, sequence_to_record,
                       annotate_record, write_record)

def batch_domestication(records, target,
                        domesticator=None,
                        standard=None,
                        allow_edits=False,
                        domesticated_suffix="",
                        include_optimization_reports=True,
                        include_original_records=True,
                        barcodes=(), barcode_order='same_as_records',
                        barcode_spacer='AA', logger="bar"):
    """Domesticate a batch of parts according to some domesticator/standard.
    
    Examples
    --------

    >>> from genedom import BUILTIN_STANDARDS, batch_domestication
    >>> batch_domestication(some_records, standard=BUILTIN_STANDARDS.EMMA)

    Parameters
    ----------

    records
      List of Bioython records to be domesticated
    
    target
      Path to a folder, to a zip file, or "@memory" for in-memory report
      generatio (the raw binary data of a zip archive is then returned)

    domesticator
      Either a single domesticator, to be used for all parts in the batch, or
      a function f(record) => appropriate_domesticator. Note that a "standard"
      can be provided instead

    standard
      A StandardDomesticatorsSet object which will be used to attribute a
      specific domesticator to each part. See BUILTIN_STANDARDS for
      examples.

    allow_edits
      If False, sequences cannot be edited by the domesticator, only extended
      with flanks. If a sequence has for instance forbidden restriction sites,
      the domesticaton will fail for this sequence (and this will be noted in
      the report. 

    domesticated_suffix
      Suffix to give to the domesticated parts names to differentiate them from
      the original parts (this is optional).

    include_optimization_reports
      If yes, some genbanks and pdfs will be produced to show how each part
      was domesticated. This is in particular informative when a domestication
      fails and you want to understand why.

    include_original_records
      Will include the input records into the final report folder/archive,
      for traceability.

    barcodes
      Either a list [(barcode_name, barcode),...] or a dictionary {name: bc} or
      a BarcodesCollection instance. If any of this is provided, the final
      parts will have a barcode added on the left (this barcode will be
      "outside" the part and won't appear in final constructs, but can be used
      to check that the part is the one you think if your samples get mixed up).
      Note that if there are less barcodes than parts, the barcodes will cycle
      and several parts may get the same barcode (which is generally fine).
    
    barcode_order
      Either "same_as_records", or "by_size" if you want your barcodes to be
      attributed from the smallest to the longest part in the batch.

    barcode_spacer
      Sequence to appear between the barcode and the left flank of the
      domesticated part.
    
    logger
      Either "bar" or None for no logger or any Proglog ProgressBarLogger.

    """
    logger = proglog.default_bar_logger(logger, min_time_interval=0.2)
    root = flametree.file_tree(target, replace=True)
    domesticated_dir = root._dir("domesticated")
    if include_original_records:
        original_dir = root._dir("original")
    if include_optimization_reports:
        errors_dir = root._dir("error_reports")
    if standard is not None:
        domesticator = standard.record_to_domesticator
    
    if hasattr(barcodes, 'items'):
        barcodes = list(barcodes.items())
    if len(barcodes):
        barcodes = [b for b, r in zip(itertools.cycle(barcodes), records)]
    if barcode_order == 'by_size':
        lengths = [len(r) for r in records]
        barcodes = [b for _, b in sorted(zip(lengths, barcodes))]

    infos = []

    domesticators = set()
    nfails = 0
    domesticated_records = []
    columns = ["Record", "Ordering Name", "Domesticator",
               "Domesticated Record", "Added bp", "Edited bp"]

    # DOMESTICATE ALL PARTS, APPEND BARCODE, GATHER DATA
    for i, record in logger.iter_bar(record=list(enumerate(records))):
        record = deepcopy(record)
        original_id = record.id
        domesticated_id = record.id + domesticated_suffix
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
        if len(barcodes):
            barcode = barcodes[i]
            if not isinstance(barcode, str):
                barcode_id, barcode = barcode
                barcode_id = " " + barcode_id
            else:
                barcode_id = ""
            barcode = sequence_to_record(barcode)
            annotate_record(barcode, label="BARCODE" + barcode_id)
        else:
            barcode = None
        # final, edits, report, success, msg
        domestication_results = record_domesticator.domesticate(
            record, report_target=report_target, edit=allow_edits)
        if not domestication_results.success:
            nfails += 1
        if barcode is not None:
            domestication_results.record_after = (
                barcode + barcode_spacer + domestication_results.record_after)
        domestication_results.record_after.original_id = original_id
        domestication_results.record_after.id = domesticated_id.replace(' ', '_')
        SeqIO.write(domestication_results.record_after,
                    domesticated_dir._file(domesticated_file_name),
                    "genbank")
        domesticated_records.append(domestication_results.record_after)
        if include_original_records:
            write_record(domestication_results.record_after,
                         original_dir._file(original_id + ".gb"))
        n_edits = domestication_results.number_of_edits()
        added_bp = len(domestication_results.record_after) - len(record)
        before_seqicon = sequenticon(record, output_format="html_image")
        after_seqicon = sequenticon(domestication_results.record_after,
                                    output_format="html_image")

        infos.append({
            "id": original_id,
            "Record": before_seqicon + original_id,
            "Domesticator": record_domesticator.name,
            "Domesticated Record": ("Failed: " + domestication_results.message)
                                   if not domestication_results.success
                                   else (after_seqicon + domesticated_id),
            "Added bp": added_bp,
            "Edited bp": n_edits
        })
        if barcode is not None:
            infos[-1]['Barcode'] = barcode_id
    
    # WRITE PDF REPORT

    sanitizing_table = sanitize_and_uniquify([info['id'] for info in infos])
    order_id_dataframe = pandas.DataFrame(list(sanitizing_table.items()),
                                          columns=["sequence", "order_id"])
    order_id_dataframe.to_csv(root._file("order_ids.csv").open('w'),
                              index=False)
    for info in infos:
        info['Order ID'] = sanitizing_table[info['id']]
    columns = ["Record", "Order ID", "Domesticator", "Domesticated Record",
               "Added bp", "Edited bp"]
    if "Barcode" in infos[0]:
        columns.append("Barcode")
    infos_dataframe = pandas.DataFrame(infos, columns=columns)
    infos_dataframe.sort_values('Order ID', inplace=True)
    domesticators = sorted(domesticators, key=lambda d: d.name)
    write_pdf_domestication_report(root._file("Report.pdf"), infos_dataframe,
                             domesticators)
    
    # WRITE THE SEQUENCES TO ORDER AS FASTA

    order_dir = root._dir('sequences_to_order', replace=True)
    for r in domesticated_records:
        r.id = sanitizing_table[r.original_id]
        r.name = ''
        r.description = ''
    SeqIO.write(domesticated_records, order_dir._file("sequences_to_order.fa"),
                "fasta")

    # WRITE THE SEQUENCES TO ORDER AS EXCEL
    
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
