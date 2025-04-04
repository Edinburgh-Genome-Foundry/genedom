.. raw:: html

    <p align="center">
    <img alt="logo" title="Genedom Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/genedom/master/docs/_static/images/logo.png" width="550">
    <br /><br />
    </p>

.. image:: https://github.com/Edinburgh-Genome-Foundry/genedom/actions/workflows/build.yml/badge.svg
    :target: https://github.com/Edinburgh-Genome-Foundry/genedom/actions/workflows/build.yml
    :alt: GitHub CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/genedom/badge.svg?branch=master
  :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/genedom?branch=master



Genedom is a Python package for managing the domestication of genetic parts,
which means modifying their sequences to make them compatible with a given
genetic assembly standard. Genedom binds together a
`sequence optimizer <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_,
information on the genetic standard, and a reporting routine to automate the
domestication of large batches in an easy and human-friendly way.

.. raw:: html

    <p align="center">
    <img alt="schema" title="schema" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/genedom/master/docs/_static/images/domestication_schema.png" width="800">
    <br /><br />
    </p>

Features include:

- User-defined part domesticators with extra nucleotides optionally added to both ends,
  hard constraints on the sequence (such as enforcing absence of a restriction
  site) and optimization objectives (such as codon optimization).
- Built-in pre-defined domesticators for popular genetic assembly standards
  (well, only EMMA at the moment).
- Generation of barcodes that can be added to the sequence
  (but won't be in final constructs). This allows easy verification of a sequence file
  or a DNA sample (e.g. in case of label mix-up).
- Routine for batch-domesticating sequences with report generation, including
  reports on each sequence optimization, spreadsheets of parts, ready-to-order FASTA
  and Genbank files of the parts, and a summary report to quickly verify everything,
  with a list of every domesticator used, for traceability.

An example summary report:

.. raw:: html

    <p align="center">
    <img alt="report" title="report" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/genedom/master/docs/_static/images/report_screenshot.png" width="600">
    <br /><br />
    </p>


You can also use Genedom online via EGF's `Domesticate Part Batches web app <https://cuba.genomefoundry.org/domesticate_part_batches>`_.

Usage examples
--------------

Simple domestication of one part
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from genedom import (GoldenGateDomesticator, random_dna_sequence,
                         write_record)
    sequence = random_dna_sequence(2000, seed=123)
    domesticator = GoldenGateDomesticator("ATTC", "ATCG", enzyme='BsmBI')
    domestication_results = domesticator.domesticate(sequence, edit=True)
    print (domestication_results.summary())
    write_record(domestication_results.record_after, 'domesticated.gb')


Generating a collection of 20bp barcodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(see docs for more options)

.. code:: python

    from genedom import BarcodesCollection

    barcodes_collection = BarcodesCollection.from_specs(
        n_barcodes=96, barcode_length=20,
        forbidden_enzymes=('BsaI', 'BsmBI', 'BbsI'))

    barcodes_collection.to_fasta('example_barcodes_collection.fa')


Domesticating a batch of parts with PDF report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from genedom import BUILTIN_STANDARDS, load_record, batch_domestication

    records = [
        load_record(filepath, name=filename)
        for filepath in records_filepaths
    ]
    barcodes_collection = BarcodesCollection.from_specs(n_barcodes=10)
    batch_domestication(records, 'domestication_report.zip',
                        barcodes=barcodes,  # optional
                        standard=BUILTIN_STANDARDS.EMMA)


Installation
------------

You can install Genedom through PIP:

.. code:: shell

    pip install genedom


License = MIT
-------------

Genedom is a free/libre and open-source software
`released on Github <https://github.com/Edinburgh-Genome-Foundry/genedom>`_ under
the MIT license (Copyright 2018 Edinburgh Genome Foundry, University of Edinburgh).
It was originally written by `Zulko <https://github.com/Zulko>`_ and is currently
being developed by `Peter Vegh <https://github.com/veghp>`_.
Everyone is welcome to `contribute <https://github.com/Edinburgh-Genome-Foundry/HowTo/blob/master/EGF/CONTRIBUTING.md>`_!
