import os
import matplotlib
matplotlib.use("Agg")
from genedom import (load_records, batch_domestication, BUILTIN_STANDARDS,
                     BarcodesCollection)

DATA_DIR = os.path.join("tests", "data")

def test_domestication_batch(tmpdir):
    records = load_records(os.path.join(DATA_DIR, "example_sequences.fa"))
    output_target = os.path.join(str(tmpdir), "test_report")
    nfails, _ = batch_domestication(records, output_target,
                                    standard=BUILTIN_STANDARDS.EMMA,
                                    allow_edits=True)
    assert (nfails == 0)
    nfails, _ = batch_domestication(records, output_target,
                                    standard=BUILTIN_STANDARDS.EMMA,
                                    allow_edits=False)
    assert (nfails > 0)

def test_barcodes_collections(tmpdir):
    # TODO: test more internediary results
    barcodes = BarcodesCollection.from_specs(n_barcodes=10)
    fasta = barcodes.to_fasta(path=os.path.join(str(tmpdir), 'test.fa'))
    barcodes_records = barcodes.to_records()

    records = load_records(os.path.join(DATA_DIR, "example_sequences.fa"))
    output_target = os.path.join(str(tmpdir), "test_report")
    nfails, _ = batch_domestication(records, output_target,
                                    standard=BUILTIN_STANDARDS.EMMA,
                                    allow_edits=True,
                                    barcodes=barcodes)
    assert (nfails == 0)
    
