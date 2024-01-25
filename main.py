from pathlib import Path


def main():
    ma_dir = "/data/normalized_microarray_donor9861/"
    microarray_file_path = Path(ma_dir + "MicroarrayExpression.csv")
    probes_file_path = Path(ma_dir + "Probes.csv")
    sample_file_path = Path(ma_dir + "SampleAnnot.csv")


if __name__ == '__main__':
    main()
