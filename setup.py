import ez_setup

ez_setup.use_setuptools()

from setuptools import setup, find_packages

version = {}
with open("genedom/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="genedom",
    version=version["__version__"],
    author="Zulko",
    description="Genetic part standardization",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    url="https://github.com/Edinburgh-Genome-Foundry/genedom",
    keywords="genetic DNA part standardization synthetic biology",
    packages=find_packages(exclude="docs"),
    include_package_data=True,
    install_requires=(
        "dnachisel[reports]",
        "snapgene_reader",
        "pdf_reports",
        "pandas",
        "dna_features_viewer",
        "python-box",
        "openpyxl",
        "flametree",
        "sequenticon",
    ),
)
