from genedom import GoldenGateDomesticator
from dnachisel import random_dna_sequence


sequence = random_dna_sequence(2000, seed=123)
domesticator = GoldenGateDomesticator("ATTC", "ATCG")
final_record, edits_record, data = domesticator.domesticate(sequence)

translator = SpecAnnotationsTranslator()
gr_record = translator.translate_record(edits_record)
gr_record.plot()
