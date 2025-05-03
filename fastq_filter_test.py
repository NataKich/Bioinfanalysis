import unittest
import fastq_filter as genom
import tempfile
import os


class DnaTests(unittest.TestCase):
    def test_check_alphabet(self):
        dna = genom.DNASequence("ATGC")
        self.assertTrue(dna.check_alphabet())
        dna = genom.DNASequence("ATGCU")
        self.assertFalse(dna.check_alphabet())

    def test_complement(self):
        dna = genom.DNASequence("ATGC")
        compl = dna.complement()
        self.assertEqual(compl.seq, "TACG")

    def test_reverse(self):
        dna = genom.DNASequence("ATGC")
        rev = dna.reverse()
        self.assertEqual(rev.seq, "CGTA")


class RnaTests(unittest.TestCase):
    def test_check_alphabet(self):
        rna = genom.RNASequence("AUGC")
        self.assertTrue(rna.check_alphabet())
        rna = genom.RNASequence("ATGCU")
        self.assertFalse(rna.check_alphabet())

    def test_complement(self):
        rna = genom.RNASequence("AUGC")
        compl = rna.complement()
        self.assertEqual(compl.seq, "UACG")

    def test_reverse(self):
        rna = genom.RNASequence("AUGC")
        rev = rna.reverse()
        self.assertEqual(rev.seq, "CGUA")


class FilterTests(unittest.TestCase):
    seq = """@SRR1363257.37 GWZHISEQ01:153:C1W31ACXX:5:1101:14027:2198 length=101
GGTTGCAGATTCGCAGTGTCGCTGTTCCAGCGCATCACATCTTTGATGTTCACGCCGTGGCGTTTAGCAATGCTTGAAAGCGAATCGCCTTTGCCCACACG
+
@?:=:;DBFADH;CAECEE@@E:FFHGAE4?C?DE<BFGEC>?>FHE4BFFIIFHIBABEECA83;>>@>@CCCDC9@@CC08<@?@BB@9:CC#######
"""

    def test_filter_no_file(self):
        with self.assertRaises(Exception):
            # where is no such file
            genom.filter_fastq("---------------")

    def test_filter_id(self):
        name_in = tempfile.mktemp(suffix=".fastq")
        name_out = tempfile.mktemp(suffix=".fastq")
        try:
            with open(name_in, "w") as f:
                f.write(self.seq)
            genom.filter_fastq(name_in, name_out)
            self.assertTrue(os.path.isfile(name_out))

            with open(name_out, "r") as f:
                res = f.read()
            self.assertTrue(res.strip() == self.seq.strip())
        finally:
            if os.path.isfile(name_in):
                os.remove(name_in)
            if os.path.isfile(name_out):
                os.remove(name_out)
        
    def test_filter_all(self):
        name_in = tempfile.mktemp(suffix=".fastq")
        name_out = tempfile.mktemp(suffix=".fastq")
        try:
            with open(name_in, "w") as f:
                f.write(self.seq)
            genom.filter_fastq(name_in, name_out, length_bounds=(102, 1000))
            self.assertTrue(os.path.isfile(name_out))

            with open(name_out, "r") as f:
                res = f.read()
            self.assertTrue(res.strip() == "")
        finally:
            if os.path.isfile(name_in):
                os.remove(name_in)
            if os.path.isfile(name_out):
                os.remove(name_out)


if __name__ == "__main__":
    unittest.main()