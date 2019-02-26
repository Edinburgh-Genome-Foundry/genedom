from genedom import BarcodesCollection

barcodes_collection = BarcodesCollection.from_specs(
    n_barcodes=96, barcode_length=20, spacer='AA',
    forbidden_enzymes=('BsaI', 'BsmBI', 'BbsI'),
    barcode_tmin=55, barcode_tmax=70,
    other_primer_sequences=(), heterodim_tmax=5,
    max_homology_length=10, include_spacers=True,
    names_template="B_%03d")

barcodes_collection.to_fasta('example_barcodes_collection.fa')