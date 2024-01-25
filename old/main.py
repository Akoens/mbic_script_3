"""
    A Large File Analysis
"""

import timeit
from functools import reduce

import pandas as pd

MICROARRAY_PATH = "../data/normalized_microarray_donor9861/"
MICROARRAY_FILE = MICROARRAY_PATH + "MicroarrayExpression.csv"
ONTOLOGY_FILE = MICROARRAY_PATH + "Ontology.csv"
PA_CALL_FILE = MICROARRAY_PATH + "PACall.csv"
PROBES_FILE = MICROARRAY_PATH + "Probes.csv"
SAMPLE_ANNOT_FILE = MICROARRAY_PATH + "SampleAnnot.csv"

CHUNK_SIZE = 1000  # ~1.6 mb


def analyze_microarray(probes_filepath: str, microarray_filepath: str, sample_annot_filepath: str,
                       cutoff: float, define_probes_by: str, define_samples_by: str) -> None:
    # == Access CSV Data Files ==
    sample_annot_io = pd.read_csv(sample_annot_filepath)
    probes_io = pd.read_csv(probes_filepath, index_col=0, chunksize=CHUNK_SIZE)
    microarray_io = pd.read_csv(microarray_filepath, index_col=0, chunksize=CHUNK_SIZE)

    averages = pd.DataFrame(columns=["average"])
    samples = sample_annot_io.loc[:, [define_samples_by]].transpose()
    uniques = pd.DataFrame(samples)

    # == Chunk Process ==
    for probes_chunk, microarray_chunk in zip(probes_io, microarray_io):
        # Help interpreter with data type
        microarray_chunk: pd.DataFrame
        probes_chunk: pd.DataFrame

        chunk_average = pd.DataFrame(columns=["average"])
        chunk_average["average"] = microarray_chunk.mean(axis=1)
        chunk_average.index = microarray_chunk.index
        chunk_average.index.name = "probe_id"

        # Merge averages with probe Data
        chunk_probe_id = probes_chunk.loc[:, ['gene_id', define_probes_by]]

        # Add Averages to previously found averages.
        averages = reduce(lambda left, right: pd.merge(left, right),
                          [averages, chunk_probe_id, chunk_average])
        # Sort on highest average, drop duplicates.
        averages = averages.sort_values(by='average', ascending=False).drop_duplicates('gene_id') \
            .sort_index()

        # Filter on Cuttoff Values
        cuts = microarray_chunk[microarray_chunk > cutoff]
        # print(cuts)

        # Find Uniques

        # Find shared

    print(averages)
    # Use the cutoff value to ascertain the regions that contain the probed gene.

    # == Post Process ==


def main():
    analyze_microarray(PROBES_FILE, MICROARRAY_FILE, SAMPLE_ANNOT_FILE, 2.0,
                       "probe_name", "structure_acronym")


if __name__ == '__main__':
    timed = timeit.timeit(main, number=1)
    print(f"{round(timed, 2)} sec")
