"""
    This script is used to analyze a microarray dataset to find genes.

    Created on Wed 1-11-2023 09:21:00
    Author: A F Koens
    Version: 0.0.1
    Title: PDB Analyzer

    Usage:
    main.py <first_structure> <second_structure> <threshold> <micro-file> <probes-file>
    <sample-file>

"""

import csv
import itertools
from docopt import docopt
import logging
from collections import namedtuple
from enum import Enum
from pathlib import Path

# Can make this a DataClass as well.
Probe = namedtuple("Probe", ["probe_id", "probe_name", "gene_id",
                             "gene_symbol", "gene_name", "entrez_id", "chromosome"])


class ProbeFilter(Enum):
    """
    Standardizing the Probe column indexes.
    """
    PROBE_ID = 0
    PROBE_NAME = 1
    GENE_ID = 2
    GENE_SYMBOL = 3
    ENTREZ_ID = 4
    CHROMOSOME = 6


class SampleFilter(Enum):
    """
    Standardizing the Sample Annotation column indexes.
    """
    STRUCTURE_ID = 0
    SLAB_NUM = 1
    WELL_ID = 2
    SLAB_TYPE = 3
    STRUCTURE_ACRONYM = 4
    STRUCTURE_NAME = 5
    POLYGON_ID = 6
    MRI_VOXEL_X = 7
    MRI_VOXEL_Y = 8
    MRI_VOXEL_Z = 9
    MNI_X = 10
    MNI_Y = 11
    MNI_Z = 12


def calc_average(*values: float or int) -> float:
    """
    Calculate averages an iterable of numbers.
    :param values: Values to get the average of.
    :return: An average.
    """
    return sum(values) / len(values)


def get_sample_index(sample_file_path: Path, structures: list[str]) -> dict[str, list[int]]:
    """
    Find the indices of the structures in the sample file.
    :param sample_file_path: Path to Sample Annotation.
    :param structures: The structures to find the indices of.
    :return: A list of indicies found for each structure stored in a dict.
    """
    # Open File
    with sample_file_path.open(mode="r", encoding="UTF-8") as sample_file:
        # Pass through build in csv parser because csv has bad standards.
        sample_csv = csv.reader(sample_file)

        # Skip headers
        next(sample_csv)

        # Initialize dict.
        brain_region_indices = {structure: [] for structure in structures}

        # Loop over lines to find structures in structure acronyms.
        for line_index, sample_line in enumerate(sample_csv):
            for structure in structures:
                if structure == sample_line[SampleFilter.STRUCTURE_ACRONYM.value]:
                    brain_region_indices[structure].append(line_index - 1)

    return brain_region_indices


def analyze_microarray(microarray_file_path: Path, probes_file_path: Path,
                       sample_indices: dict[str, list[int]], threshold: float = 5.0):
    """
    Find the highest averaged probes with a value that is higher than the threshold for atl east
    one of the structures in the sample indices.
    :param microarray_file_path: Path to microarray file.
    :param probes_file_path: Path to probes file.
    :param sample_indices: Indicies of the structures gathered from the sample.
    :param threshold: A threshold value for constraining the output.
    :return: A dict of probes that are the highest average and a dict of array values for the
    samples.
    Both share a key for the gene id.
    """

    averages: dict[str, float] = {}
    probes: dict[str, Probe] = {}
    array: dict[str, list[list[float]]] = {}

    # Open Files.
    with (microarray_file_path.open(mode="r", encoding="UTF=8") as micro_file,
          probes_file_path.open(mode="r", encoding="UTF-8") as probes_file):

        # Convert to CSV object.
        micro_csv = csv.reader(micro_file)
        probes_csv = csv.reader(probes_file)

        # Skip headers
        next(probes_csv)

        # Compare data from lines.
        for micro_line, probes_line in zip(micro_csv, probes_csv):
            if micro_line[0] != probes_line[0]:
                logging.error("Error: Files are not ordered. Could not match Probe IDs")

            micro_values = list(map(float, micro_line[1:]))
            average = calc_average(*micro_values)

            probe = Probe(*probes_line)

            # Can be used to determine cutoff based on the array values. Not the average.
            brain_region_array_values = {region: [micro_values[index] for index in region_indices]
                                         for region, region_indices in sample_indices.items()}

            # Find the highest average and store the corresponding probe.
            if probe.gene_id in probes:
                highest_average = max(averages[probe.gene_id], average)
                if averages[probe.gene_id] < highest_average:
                    probes[probe.gene_id] = probe
                    averages[probe.gene_id] = highest_average
                    array[probe.gene_id] = [[micro_values[index] for index in indices] for indices
                                            in sample_indices.values()]
            # Add values if they pass the threshold.
            elif any(value >= threshold for value in
                     itertools.chain(*brain_region_array_values.values())):
                probes[probe.gene_id] = probe
                averages[probe.gene_id] = average
                array[probe.gene_id] = [[micro_values[index] for index in indices] for indices in
                                        sample_indices.values()]

    return probes, array


def find_uniques_shared(array: dict[str, list[list[float]]], structures: list[str],
                        threshold: float):
    """
    Finds the unique genes and the shared genes per structure,
    :param array:
    :param structures:
    :param threshold:
    :return:
    """
    shared = []
    unique = {structure: [] for structure in structures}
    for gene, regions in array.items():
        condition = [any(value >= threshold for value in values) for values in regions]
        if all(condition):
            shared.append(gene)
        else:
            index = condition.index(True)
            unique[structures[index]].append(gene)

    return unique, shared


def print_findings(probes: dict[str, Probe], uniques: dict[str, list[str]], shared: list[str],
                   probe_filter: ProbeFilter) -> None:
    """
    Prints the final findings in a human-readable way.
    :param probes: The probes that are the highest average above the threshold.
    :param uniques: The unique genes per structure.
    :param shared: The shared genes.
    :param probe_filter: The filter for the data that the user wants.
    :return: None
    """
    for structure, genes in uniques.items():
        print(
            f"For {structure},  {len(genes)} unique genes found: "
            f"{[probes[gene_id][probe_filter.value] for gene_id in genes]}")
    print(f"Shared genes are: {[probes[gene_id][probe_filter.value] for gene_id in shared]}\n")


def main() -> None:
    """
    Gather user commandline arguments and parse them to the analyzer function.
    :return: None
    """
    args = docopt(__doc__, version='PDB Analyzer 2.0')

    microarray_file_path = Path(args["<micro-file>"])
    probes_file_path = Path(args["<probes-file>"])
    sample_file_path = Path(args["<sample-file>"])
    structures = [args["<first_structure>"], args["<second_structure>"]]
    threshold = float(args["<threshold>"])

    print("Starting Analysis...\n")

    samples = get_sample_index(sample_file_path, structures)
    probes, array = analyze_microarray(microarray_file_path, probes_file_path, samples, threshold)
    uniques, shared = find_uniques_shared(array, structures, threshold)

    print_findings(probes, uniques, shared, ProbeFilter.GENE_SYMBOL)
    print("Analysis Complete.")


if __name__ == '__main__':
    main()
