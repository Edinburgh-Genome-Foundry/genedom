.. raw:: html

    <p align="center">
    <img alt="lala Logo" title="genedom Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/genedom/master/docs/logo.png" width="550">
    <br /><br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/genedom.svg?branch=master
  :target: https://travis-ci.org/Edinburgh-Genome-Foundry/genedom
  :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/genedom/badge.svg?branch=master
  :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/genedom?branch=master



GeneDom is a python library for managing the domestication of genetic parts so that they will comply to a genetic assembly standard. It puts together a `sequence optimizer <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_, genetic standards data, and a reporting system, to automate the domestication of large batches in an easy and human-friendly way.

Features include:
- Possibility to define parts domesticators with added right-hand and left-hand nucleotide, hard constraints on the sequence (such as absence of a restriction site) and optimization objectives (such as codon optimization).
- Built-in pre-defined domesticators for popular genetic assembly standards (well, only EMMA at the moment).
- Routine for mass-domesticating sequences with report generation

Genedom is very much in development. There will be much better docs later.

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
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/Geneblocks>`_ under the MIT licence (copyright Edinburgh Genome Foundry).
Everyone is welcome to contribute !
