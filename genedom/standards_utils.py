from .builtin_domesticators import BUILTIN_DOMESTICATORS

def list_standard_overhangs(standard_name):
    domesticators = list(BUILTIN_DOMESTICATORS[standard_name].values())
    overhangs = [domesticators[0].left_overhang]
    for d in domesticators:
        if d.right_overhang not in overhangs:
            overhangs.append(d.right_overhang)
    return overhangs
