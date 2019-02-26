import os
import matplotlib
matplotlib.use("Agg")
from genedom import (load_records, batch_domestication, BUILTIN_STANDARDS,
                     BarcodesCollection, GoldenGateDomesticator,
                     random_dna_sequence, load_record)

DATA_DIR = os.path.join("tests", "data")
PARTS_DIR = os.path.join("tests", "data", "example_parts")

def test_basic_domestication():
    sequence = random_dna_sequence(2000, seed=123)
    domesticator = GoldenGateDomesticator("ATTC", "ATCG")
    domestication_results = domesticator.domesticate(sequence, edit=True)
    print (domestication_results.summary())
    assert domestication_results.summary().startswith('SUCCESS')

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
    records = [
        load_record(os.path.join(PARTS_DIR, filename), name=filename)
        for filename in os.listdir(PARTS_DIR)
    ]
    output_target = os.path.join(str(tmpdir), "test_report")
    nfails, _ = batch_domestication(records, output_target,
                                    standard=BUILTIN_STANDARDS.EMMA,
                                    allow_edits=True,
                                    barcodes=barcodes)
    assert (nfails == 0)
    
def test_BarcodesCollection(tmpdir):
    barcodes = BarcodesCollection.from_specs(n_barcodes=10)
    barcodes.to_fasta(os.path.join(str(tmpdir), 'test.fa'))
    barcodes.to_records(str(tmpdir))