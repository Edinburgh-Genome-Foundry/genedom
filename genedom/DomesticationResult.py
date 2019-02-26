class DomesticationResult:
    """Class to contain and represent one result of a part domestication.

    Parameters
    ----------
    record_before
      Biopython record of the sequence before it is domesticated. 
    
    record_after
      Biopython record of the sequence after it was domesticated
    
    edits_record
      Biopython record annotated with every mutation introduced by the
      domestication
    
    report_data
      Raw binary data of a PDF report
    
    success
      Boolean indicating whether the domestication was possible or some
      constraints could not be satisfied.
      
    message
      String containing some information on the reason for failure.
    """
    
    def __init__(self, record_before, record_after, edits_record,
                 report_data, success, message):
        self.record_before = record_before
        self.record_after = record_after
        self.edits_record = edits_record
        self.report_data = report_data
        self.success = success
        self.message = message
    
    def summary(self):
        """Return a string summarizing how the domestication went.

        Either "SUCCESS - 42bp edited" or "FAILURE - (some help message)"
        """
        if self.success:
            return "SUCCESS - %d nucleotides edited." % self.number_of_edits()
        else:
            return "FAILURE - %s" % self.message
    
    def number_of_edits(self):
        return sum([len(f) for f in self.edits_record.features
                    if f.qualifiers.get("is_edit", False)])