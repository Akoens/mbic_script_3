import sys


def is_data_ordered(probes_filepath: str, microarray_filepath: str):
    match = True
    with (open(probes_filepath, 'r', encoding="UTF-8") as probes_file, \
          open(microarray_filepath, 'r', encoding="UTF-8") as microarray_file):
        # Skip headers
        probes_file.readline()

        # Loop over lines
        for probe_line, microarray_line in zip(probes_file, microarray_file):

            # Slit line into column segments.
            probe_line = probe_line.split(',')
            microarray_line = microarray_line.split(',')

            # If probe_id doesn't match
            if int(probe_line[0]) != int(microarray_line[0]) and match:
                print(probe_line[0], microarray_line[0])
                match = False
    print("Data is ordered." if match else "Data is NOT ordered.")


def print_csvfile_size_data(filepath: str):
    with open(filepath, 'r', encoding="UTF-8") as file:
        len_total = 0
        size_total = 0

        index = 0
        line = ""
        for index, line in enumerate(file):
            len_total += len(line)
            size_total += sys.getsizeof(line)

        cols = len(line.split(","))
        print(f"{filepath} \n"
              f"Contains: {index} lines \n"
              f"Contains: {cols} columns\n"
              f"Average line len: {round(len_total / index, 2)} Characters\n"
              f"Average line size: {round(size_total / index, 2)} Bytes\n"
              f"Total file size: {size_total} Bytes\n")
