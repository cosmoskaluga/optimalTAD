import pandas as pd
import pytest
import filecmp
import os

def iterate_through_files(ref_path, test_path, relative_path):
    for file in os.listdir(ref_path):
        ref_file_path = os.path.join(ref_path, file)
        if os.path.isfile(ref_file_path):
            test_file_path = os.path.join(test_path, relative_path, file)

            if file == "amplitudes.csv":
                test_df = pd.read_csv(test_file_path)
                ref_df = pd.read_csv(ref_file_path)

                assert test_df.iloc[:, :2].equals(ref_df.iloc[:, :2])
            else:
                assert filecmp.cmp(test_file_path, ref_file_path, shallow=False)
        else:
            iterate_through_files(os.path.join(ref_path, file), test_path, os.path.join(relative_path, file))

def test_check_fly():
    os.system("optimalTAD check")

    ref_path = "./tests/reference_files/"
    test_path = "./optimalTAD/testouput/"
    relative_path = ""

    iterate_through_files(ref_path, test_path, relative_path)