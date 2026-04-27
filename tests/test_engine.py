import unittest

from src.engine import (
    process_sequence,
    clean_sequence,
    validate_sequence,
    detect_sequence_type,
    reverse_transcribe,
    get_reverse_complement,
    calculate_base_percentages,
)


class TestProcessSequence(unittest.TestCase):
    def test_template_dna_transcribes_correctly(self):
        result = process_sequence("TACGGC", strand_type="template")
        self.assertIsNone(result["error"])
        self.assertEqual(result["mrna"], "AUGCCG")

    def test_translation_stops_at_internal_stop_codon(self):
        result = process_sequence("ATGAAATAATTT", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertEqual(result["mrna"], "AUGAAAUAAUUU")
        self.assertEqual(result["proteinChain"], "Met-Lys")
        self.assertTrue(result["encounteredStop"])

    def test_batch_like_large_sequence_input_remains_valid(self):
        dna = "ATGC" * 200
        result = process_sequence(dna, strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertGreater(len(result["mrna"]), 0)

    def test_empty_input_returns_error(self):
        result = process_sequence("")
        self.assertIsNotNone(result["error"])
        self.assertEqual(result["error"], "DNA sequence is required.")

    def test_invalid_characters_returns_error(self):
        result = process_sequence("ATGXYZ")
        self.assertIsNotNone(result["error"])
        self.assertEqual(result["error"], "DNA must only contain A, T, C, and G.")

    def test_coding_strand_replaces_t_with_u(self):
        result = process_sequence("ATGCTT", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertEqual(result["mrna"], "AUGCUU")

    def test_gc_content_calculation(self):
        result = process_sequence("GGCC", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertAlmostEqual(result["gcContent"], 100.0)

    def test_gc_content_fifty_percent(self):
        result = process_sequence("AATTGGCC", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertAlmostEqual(result["gcContent"], 50.0)

    def test_molecular_weight_nonzero_for_single_codon(self):
        result = process_sequence("AAA", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertGreater(result["molecularWeight"], 0)

    def test_one_letter_protein_output(self):
        result = process_sequence("ATGAAATAATTT", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertEqual(result["oneLetterProtein"], "MK")

    def test_process_sequence_includes_reverse_transcription(self):
        result = process_sequence("ATGCTT", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertEqual(result["reverseTranscription"], "ATGCTT")

    def test_process_sequence_includes_reverse_complement(self):
        result = process_sequence("ATGCTT", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertEqual(result["reverseComplement"], "AAGCAT")

    def test_process_sequence_includes_base_percentages(self):
        result = process_sequence("AAGG", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertIn("basePercentages", result)
        self.assertAlmostEqual(result["basePercentages"]["A"], 50.0)
        self.assertAlmostEqual(result["basePercentages"]["G"], 50.0)


class TestProcessSequenceFasta(unittest.TestCase):
    def test_fasta_with_grch38_header(self):
        fasta = ">4 dna:chromosome chromosome:GRCh38:4:122868001:122946006:1\nATGCTT"
        result = process_sequence(fasta, strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertEqual(result["dna"], "ATGCTT")
        self.assertEqual(result["mrna"], "AUGCUU")

    def test_fasta_multiline_sequence(self):
        fasta = ">seq1\nATGC\nATGC"
        result = process_sequence(fasta, strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertEqual(result["dna"], "ATGCATGC")

    def test_plain_nucleotides_still_work(self):
        result = process_sequence("ATGCTT", strand_type="coding")
        self.assertIsNone(result["error"])
        self.assertEqual(result["dna"], "ATGCTT")

    def test_fasta_header_only_without_sequence_returns_error(self):
        fasta = ">4 dna:chromosome chromosome:GRCh38:4:122868001:122946006:1\n"
        result = process_sequence(fasta, strand_type="coding")
        self.assertIsNotNone(result["error"])


class TestCleanSequence(unittest.TestCase):
    def test_strips_whitespace_and_uppercases(self):
        self.assertEqual(clean_sequence("  atgc  "), "ATGC")

    def test_skips_fasta_header_lines(self):
        fasta = ">seq1\nATGC\nGGCC"
        self.assertEqual(clean_sequence(fasta), "ATGCGGCC")

    def test_removes_digits_and_special_chars(self):
        self.assertEqual(clean_sequence("ATG 123 C-T"), "ATGCT")

    def test_empty_fasta_only_headers(self):
        self.assertEqual(clean_sequence(">header\n"), "")

    def test_multiline_sequence(self):
        self.assertEqual(clean_sequence("ATG\nCCC\nTAG"), "ATGCCCTAG")


class TestValidateSequence(unittest.TestCase):
    def test_valid_dna(self):
        result = validate_sequence("ATGCATGC")
        self.assertTrue(result["valid"])
        self.assertIsNone(result["error"])

    def test_empty_sequence(self):
        result = validate_sequence("")
        self.assertFalse(result["valid"])
        self.assertIn("empty", result["error"])

    def test_too_short(self):
        result = validate_sequence("AT")
        self.assertFalse(result["valid"])
        self.assertIn("too short", result["error"])

    def test_both_t_and_u(self):
        result = validate_sequence("ATUGC")
        self.assertFalse(result["valid"])
        self.assertIn("both T and U", result["error"])

    def test_invalid_characters(self):
        result = validate_sequence("ATGXYZ")
        self.assertFalse(result["valid"])
        self.assertIn("Invalid characters", result["error"])

    def test_valid_rna(self):
        result = validate_sequence("AUGCUU")
        self.assertTrue(result["valid"])

    def test_n_allowed_in_dna(self):
        result = validate_sequence("ATGCN")
        self.assertTrue(result["valid"])

    def test_n_allowed_in_rna(self):
        result = validate_sequence("AUGCN")
        self.assertTrue(result["valid"])


class TestDetectSequenceType(unittest.TestCase):
    def test_dna_detected(self):
        self.assertEqual(detect_sequence_type("ATGCATGC"), "DNA")

    def test_rna_detected(self):
        self.assertEqual(detect_sequence_type("AUGCUU"), "RNA")

    def test_ambiguous_raises(self):
        with self.assertRaises(ValueError):
            detect_sequence_type("ATUGC")

    def test_only_gc_defaults_to_dna(self):
        self.assertEqual(detect_sequence_type("GGCCAACC"), "DNA")


class TestReverseTranscribe(unittest.TestCase):
    def test_u_replaced_by_t(self):
        self.assertEqual(reverse_transcribe("AUGCUU"), "ATGCTT")

    def test_no_u_unchanged(self):
        self.assertEqual(reverse_transcribe("AAGCAA"), "AAGCAA")

    def test_empty_string(self):
        self.assertEqual(reverse_transcribe(""), "")


class TestGetReverseComplement(unittest.TestCase):
    def test_known_sequence(self):
        self.assertEqual(get_reverse_complement("ATGCTT"), "AAGCAT")

    def test_palindrome(self):
        self.assertEqual(get_reverse_complement("ATAT"), "ATAT")

    def test_single_base(self):
        self.assertEqual(get_reverse_complement("A"), "T")
        self.assertEqual(get_reverse_complement("T"), "A")
        self.assertEqual(get_reverse_complement("G"), "C")
        self.assertEqual(get_reverse_complement("C"), "G")

    def test_n_preserved(self):
        self.assertEqual(get_reverse_complement("NATG"), "CATN")


class TestCalculateBasePercentages(unittest.TestCase):
    def test_equal_distribution(self):
        result = calculate_base_percentages("ATGC")
        self.assertAlmostEqual(result["A"], 25.0)
        self.assertAlmostEqual(result["T"], 25.0)
        self.assertAlmostEqual(result["G"], 25.0)
        self.assertAlmostEqual(result["C"], 25.0)

    def test_single_base(self):
        result = calculate_base_percentages("AAAA")
        self.assertEqual(result, {"A": 100.0})

    def test_raises_on_empty(self):
        with self.assertRaises(ValueError):
            calculate_base_percentages("")

    def test_percentages_sum_to_100(self):
        result = calculate_base_percentages("AATTTGGGCC")
        total = sum(result.values())
        self.assertAlmostEqual(total, 100.0, places=1)


if __name__ == "__main__":
    unittest.main()
