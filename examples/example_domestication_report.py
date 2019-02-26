import os
from genedom import BUILTIN_STANDARDS, load_record, batch_domestication

records = [
    load_record(os.path.join('example_parts', filename), name=filename)
    for filename in os.listdir('example_parts')
]
batch_domestication(records, 'example_domestication_report.zip',
                    standard=BUILTIN_STANDARDS.EMMA)