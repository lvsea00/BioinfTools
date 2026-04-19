import os
import pytest
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from main import filter_fastq


def test_filter_by_gc():
    """
    Check GC filtration validity
    """
    filter_fastq(
        input_fastq="example_fastq.fastq",
        output_fastq="gc_test.fastq",
        gc_bounds=(40, 60),
        length_bounds=(0, 2**32),
        quality_threshold=0
    )

    output_path = os.path.join('filtered', "gc_test.fastq")
    assert os.path.exists(output_path)

    records = list(SeqIO.parse(output_path, "fastq"))
    for record in records:
        gc = gc_fraction(record.seq) * 100
        assert 40 <= gc <= 60


def test_filter_by_length():
    """
    Check length filtration validity
    """
    filter_fastq(
        input_fastq="example_fastq.fastq",
        output_fastq="length_test.fastq",
        gc_bounds=(0, 100),
        length_bounds=(20, 40),
        quality_threshold=0
    )

    output_path = os.path.join('filtered', "length_test.fastq")
    assert os.path.exists(output_path)

    records = list(SeqIO.parse(output_path, "fastq"))
    for record in records:
        assert 20 <= len(record.seq) <= 40


def test_filter_by_quality():
    """
    Check quality filtration validity
    """
    filter_fastq(
        input_fastq="example_fastq.fastq",
        output_fastq="quality_test.fastq",
        gc_bounds=(0, 100),
        length_bounds=(0, 2**32),
        quality_threshold=30
    )

    output_path = os.path.join('filtered', "quality_test.fastq")
    assert os.path.exists(output_path)

    records = list(SeqIO.parse(output_path, "fastq"))
    for record in records:
        qualities = record.letter_annotations["phred_quality"]
        mean_quality = sum(qualities) / len(qualities)
        assert mean_quality > 30


def test_empty_result():
    """
    Check function work with empty result
    """
    filter_fastq(
        input_fastq="example_fastq.fastq",
        output_fastq="empty.fastq",
        gc_bounds=(0, 1),
        length_bounds=(10000, 2**32),
        quality_threshold=100
    )

    output_path = os.path.join('filtered', "empty.fastq")
    assert os.path.exists(output_path)

    records = list(SeqIO.parse(output_path, "fastq"))
    assert len(records) == 0


def test_read_file():
    """
    Check if file reading is correct
    """
    input_file = "example_fastq.fastq"
    assert os.path.exists(input_file)
    assert os.path.getsize(input_file) > 0, f"{input_file} is empty"


def test_write_file():
    """
    Check if file writing is correct
    """
    filter_fastq(
        input_fastq="example_fastq.fastq",
        output_fastq="write_test.fastq",
        gc_bounds=(30, 70),
        length_bounds=(20, 60),
        quality_threshold=20
    )

    output_path = os.path.join('filtered', "write_test.fastq")
    assert os.path.exists(output_path), "Output file is not created"
    assert os.path.getsize(output_path) > 0, "Output file is empty"


def test_error():
    """
    Check that exception is raised when file does not exist
    """
    no_fastq = "no_fastq_file"
    with pytest.raises(Exception) as e:
        filter_fastq(
            input_fastq=no_fastq,
            output_fastq="error_test.fastq",
            gc_bounds=(0, 100),
            length_bounds=(0, 2**32),
            quality_threshold=0
        )

    assert no_fastq in str(e.value)


def test_create_dir():
    """
    Check if output directory is created
    """
    filter_fastq(
        input_fastq="example_fastq.fastq",
        output_fastq="dir_test.fastq",
        gc_bounds=(0, 100),
        length_bounds=(0, 2**32),
        quality_threshold=0
    )

    assert os.path.exists('filtered')
    assert os.path.isdir('filtered')
