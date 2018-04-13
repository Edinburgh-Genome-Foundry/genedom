import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('genedom/version.py').read()) # loads __version__

setup(name='genedom',
      version=__version__,
      author='Zulko',
    description='Genetic parts standardization',
    long_description=open('README.rst').read(),
    license='see LICENSE.txt',
    keywords="genetic DNA part standardization synthetic biology",
    packages= find_packages(exclude='docs'),
    install_requires("dnachisel", "snapgene_reader", "pdf_reports"))
