from genedom import GoldenGateDomesticator, random_dna_sequence, write_record

sequence = random_dna_sequence(5000, seed=123)
domesticator = GoldenGateDomesticator("ATTC", "ATCG")
domestication_results = domesticator.domesticate(sequence, edit=True)
print (domestication_results.summary())
write_record(domestication_results.record_after, 'basic_example_after.gb')
write_record(domestication_results.edits_record, 'basic_example_edits.gb')