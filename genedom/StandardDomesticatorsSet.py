from box import Box

class StandardDomesticatorsSet:

    def __init__(self, domesticators):
        self.domesticators = domesticators
    
    def list_overhangs(self):
        domesticators = list(self.domesticators.values())
        overhangs = [domesticators[0].left_overhang]
        for d in domesticators:
            for o in (d.right_overhang, d.left_overhang):
                if o not in overhangs:
                    overhangs.append(o)
        return overhangs

    def record_to_domesticator(self, record):
        return self.domesticators[record.id.split("_")[0]]