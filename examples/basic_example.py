from genedom import GoldenGateDomesticator, random_dna_sequence

sequence = random_dna_sequence(2000, seed=123)
domesticator = GoldenGateDomesticator("ATTC", "ATCG")
domestication_results = domesticator.domesticate(sequence, edit=True)
print (domestication_results.summary())
