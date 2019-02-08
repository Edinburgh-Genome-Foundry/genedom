.. raw:: html

    <p align="center">
    <img alt="logo" title="genedom Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/genedom/master/docs/logo.png" width="550">
    <br /><br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/genedom.svg?branch=master
  :target: https://travis-ci.org/Edinburgh-Genome-Foundry/genedom
  :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/genedom/badge.svg?branch=master
  :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/genedom?branch=master



GeneDom is a python library for managing the domestication of genetic parts
(i.e. the modification of their sequence so as to make them compatible with a
given genetic assembly standard). Genedom binds together a
`sequence optimizer <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_,
genetic standards informations, and a reporting routines, to automate the
domestication of large batches in an easy and human-friendly way.

.. raw:: html

    <p align="center">
    <img alt="schema" title="schema" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/genedom/master/docs/domestication_schema.png" width="800">
    <br /><br />
    </p>

Features include:

- Possibility to define parts domesticators with added right-hand and left-hand
  nucleotide, hard constraints on the sequence (such as absence of a restriction
  site) and optimization objectives (such as codon optimization).
- Built-in pre-defined domesticators for popular genetic assembly standards
  (well, only EMMA at the moment).
- Possibility to generate and attribute barcodes that will be added to the
  sequence (but won't be in final constructs) in order to easily check
  that this is the right part in the future in case of label mix-up. 
- Routine for mass-domesticating sequences with report generation, including
  reports on each sequence optimization, spreadsheets of parts, ready-to-order FASTA
  records of the parts, and a summary report to quickly verify everything, with a
  list of every domesticator used, for traceability.

Here is an example of summary report:

.. raw:: html

    <p align="center">
    <img alt="report" title="report" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/genedom/master/docs/report_screenshot.png" width="600">
    <br /><br />
    </p>

Genedom is still in development but works well enough for us to use it routinely. There will be much better docs later on.

Installation
-------------


You can install Genedom through PIP (coming soon)

.. code:: shell

    sudo pip install geneblocks

Alternatively, you can unzip the sources in a folder and type

.. code:: shell

    sudo python setup.py install


Licence
--------

Genedom is an open-source software originally written at the `Edinburgh Genome Foundry
<http://www.genomefoundry.org>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/genedom>`_ under the MIT licence (copyright Edinburgh Genome Foundry).
Everyone is welcome to contribute !
