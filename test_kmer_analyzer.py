import pytest
from pythonbadcode import (
    validate_sequence,
    update_kmer_count,
    count_kmers_with_context,
    write_results_to_file
)

# testing the validate function
def test_normal_sequence_works():
    assert validate_sequence("ATGTCTGTCTGAA", 2) == True

def test_sequence_too_short():
    assert validate_sequence("AT", 5) == False

def test_sequence_has_numbers():
    assert validate_sequence("ATG1TC", 2) == False

def test_sequence_same_length_as_k():
    assert validate_sequence("AT", 2) == False

def test_empty_sequence():
    assert validate_sequence("", 2) == False

# testing kmer counting
def test_new_kmer_count_should_be_one():
    kmer_data = {}
    result = update_kmer_count(kmer_data, "AT", "G")
    assert result["AT"]["count"] == 1

def test_same_kmer_twice_gives_count_two():
    kmer_data = {}
    update_kmer_count(kmer_data, "AT", "G")
    update_kmer_count(kmer_data, "AT", "G")
    assert kmer_data["AT"]["count"] == 2

def test_different_next_chars_tracked():
    kmer_data = {}
    update_kmer_count(kmer_data, "AT", "G")
    update_kmer_count(kmer_data, "AT", "C")
    assert kmer_data["AT"]["next_chars"]["G"] == 1
    assert kmer_data["AT"]["next_chars"]["C"] == 1

def test_same_next_char_twice():
    kmer_data = {}
    update_kmer_count(kmer_data, "AT", "G")
    update_kmer_count(kmer_data, "AT", "G")
    assert kmer_data["AT"]["next_chars"]["G"] == 2

# testing the full sequence analysis
def test_tg_count_in_example_sequence():
    result = count_kmers_with_context("ATGTCTGTCTGAA", 2)
    assert result["TG"]["count"] == 3

def test_at_count_in_example_sequence():
    result = count_kmers_with_context("ATGTCTGTCTGAA", 2)
    assert result["AT"]["count"] == 1

def test_tg_followed_by_correct_chars():
    result = count_kmers_with_context("ATGTCTGTCTGAA", 2)
    assert result["TG"]["next_chars"]["T"] == 2
    assert result["TG"]["next_chars"]["A"] == 1

def test_last_kmer_not_included():
    result = count_kmers_with_context("ATGTCTGTCTGAA", 2)
    assert "AA" not in result

# testing the output file
def test_output_file_gets_created(tmp_path):
    kmer_data = {"AT": {"count": 1, "next_chars": {"G": 1}}}
    output = tmp_path / "output.txt"
    write_results_to_file(kmer_data, str(output))
    assert output.exists()

def test_output_has_correct_counts(tmp_path):
    kmer_data = {
        "AT": {"count": 1, "next_chars": {"G": 1}},
        "TG": {"count": 3, "next_chars": {"T": 2, "A": 1}},
    }
    output = tmp_path / "output.txt"
    write_results_to_file(kmer_data, str(output))
    content = output.read_text()
    assert "AT 1" in content
    assert "TG 3" in content

def test_output_is_alphabetical(tmp_path):
    kmer_data = {
        "TG": {"count": 1, "next_chars": {"A": 1}},
        "AT": {"count": 1, "next_chars": {"G": 1}},
    }
    output = tmp_path / "output.txt"
    write_results_to_file(kmer_data, str(output))
    lines = output.read_text().strip().split("\n")
    assert lines[0].startswith("AT")
    assert lines[1].startswith("TG")

# testing that multiple sequences get combined properly
def test_two_sequences_combined(tmp_path):
    seq_file = tmp_path / "seqs.txt"
    seq_file.write_text("ATGC\nATGT\n")
    output = tmp_path / "out.txt"

    import sys
    sys.argv = ["pythonbadcode.py", str(seq_file), "2", str(output)]
    
    from pythonbadcode import main
    main()

    content = output.read_text()
    # AT shows up once in each sequence so total should be 2
    assert "AT 2" in content
