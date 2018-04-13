import os
import matplotlib
matplotlib.use("Agg")
from genedom import load_records, batch_domestication, BUILTIN_DOMESTICATORS

DATA_DIR = os.path.join("tests", "data")

def test_domestication_batch(tmpdir):
    records = load_records(os.path.join(DATA_DIR, "example_sequences.fa"))

    def domesticator(record):
        """Find which domesticator to use for the given part.
        Here we use the fact that all records will have an
        ID of the form ``position_partname``. We extract the
        position and return the corresponding domesticator.
        """
        position = record.id.split("_")[0]
        return BUILTIN_DOMESTICATORS.EMMA[position]
    output_target = os.path.join(str(tmpdir), "test_report")
    nfails, _ = batch_domestication(records, domesticator, output_target,
                                    allow_edits=True)
    assert (nfails == 0)
    nfails, _ = batch_domestication(records, domesticator, output_target,
                                    allow_edits=False)
