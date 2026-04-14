import unittest

from src.engine import process_sequence


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


if __name__ == "__main__":
    unittest.main()
