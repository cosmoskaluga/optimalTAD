import pytest
import os
from optimalTAD.optimization.utils import check_chrnames, check_filenames


array_of_tests_check_chrnames = [(['chr2L', 'chrX', 'chrM'], ['chr3L', 'chr2L', 'chr3R'], ['chr2L']),
								 (['chr3L', 'chr2L', 'chr3R', 'chrX', 'chrM'], ['chr3L', 'chr2L', 'chr3R'], ['chr3L', 'chr2L', 'chr3R']),
								 (['chr3L', 'chr2L', 'chr3R'], ['chr3L', 'chr2L', 'chr3R', 'chrX', 'chrM'], ['chr3L', 'chr2L', 'chr3R']),
								 (['3L', '2L', '3R'], ['chr3L', 'chr2L', 'chrX', 'chrM'], ['chr3L', 'chr2L']),
								 (['chr3L', 'chr2L', 'chrX', 'chrM'], ['3L', '2L', '3R'], ['3L', '2L'])]

@pytest.mark.parametrize("labels_config, labels, expected_labels", array_of_tests_check_chrnames)
def test_check_chrnames(labels_config, labels, expected_labels):
	assert set(check_chrnames(labels_config, labels)) == set(expected_labels)

# def test_errors_check_chrnames(labels_config, labels, expected_labels):


array_of_tests_check_filenames = [(["hic/hic_sample_1.cool", "hic/hic_sample_2.cool"], ["chip/chipseq_sample_1.bedgraph", "chip/chipseq_sample_2.bedgraph"], ["hic/hic_sample_1.cool", "hic/hic_sample_2.cool"], ["chip/chipseq_sample_1.bedgraph", "chip/chipseq_sample_2.bedgraph"], ["hic_sample_1", "hic_sample_2"]),
								  (["hic/hic_sample_1.cool"], ["chip/chipseq_sample_1.bedgraph", "chip/chipseq_sample_2.bedgraph"], ["hic/hic_sample_1.cool", "hic/hic_sample_1.cool"], ["chip/chipseq_sample_1.bedgraph", "chip/chipseq_sample_2.bedgraph"], ["chipseq_sample_1", "chipseq_sample_2"]),
								  (["hic/hic_sample_1.cool", "hic/hic_sample_2.cool"], ["chip/chipseq_sample_1.bedgraph"], ["hic/hic_sample_1.cool", "hic/hic_sample_2.cool"], ["chip/chipseq_sample_1.bedgraph", "chip/chipseq_sample_1.bedgraph"], ["hic_sample_1", "hic_sample_2"])]

@pytest.mark.parametrize("hic_files, chipseq_files, expected_hic_files, expected_chipseq_files, expected_samplenames", array_of_tests_check_filenames)
def test_check_filenames(hic_files, chipseq_files, expected_hic_files, expected_chipseq_files, expected_samplenames):
	output_hic_files, output_chipseq_file, output_samplenames = check_filenames(hic_files, chipseq_files)

	assert all(output_hic_files == expected_hic_files)
	assert all(output_chipseq_file == expected_chipseq_files)
	assert all(output_samplenames == expected_samplenames)


def test_check_config_exists():
	assert os.path.exists("config.ini") == True