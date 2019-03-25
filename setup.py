import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('genedom/version.py').read()) # loads __version__

setup(name='genedom',
      version=__version__,
      author='Zulko',
    description='Genetic parts standardization',
    long_description=open('pypi-readme.rst').read(),
    license='see LICENSE.txt',
    url='https://github.com/Edinburgh-Genome-Foundry/genedom',
    keywords="genetic DNA part standardization synthetic biology",
    packages= find_packages(exclude='docs'),
    include_package_data=True,
    install_requires=("dnachisel[reports]", "snapgene_reader", "pdf_reports",
                      "pandas", "dna_features_viewer", "python-box",
                      "openpyxl", "flametree", "sequenticon",
                      "dna_features_viewer"))
