import pytest
import os
from optimalTAD.optimization.utils import check_chrnames 


def test_check_chrnames_intersection():
	labels_config = ['chr2L', 'chrX', 'chrM']
	labels = ['chr3L', 'chr2L', 'chr3R']
	expected = ['chr2L']
	assert check_chrnames(labels_config, labels) == expected


def test_check_chrnames_no_prefix():
	labels_config = ['chr2L', 'chrX', 'chrM']
	labels = ['3L', '2L', '4', 'X']
	expected = ['2L', 'X']
	assert not set(check_chrnames(labels_config, labels)) ^ set(expected)


def test_check_chrnames_chr_prefix():
	labels_config = ['3L', '2L', '4', 'X']
	labels = ['chr2L', 'chrX', 'chrM']
	expected = ['chr2L', 'chrX']
	assert not set(check_chrnames(labels_config, labels)) ^ set(expected)


def test_check_config_exists():
	assert os.path.exists("config.ini") == True