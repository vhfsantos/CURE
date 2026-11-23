"""Regression checks for GeneRegion cases found in Poneroides analyses."""

import io
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / 'uce_kit'))

from parse_results import read_df
from gene_region_remainder import main as move_remainder
from gene_region_concat import run_job
from uce_kit import filter_uce_match_m8_file


class GeneRegionRegressionTests(unittest.TestCase):
    def test_parallel_job_preserves_parentheses_and_requires_output(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            source = folder / 'gene-Su(var)3-3'
            source.mkdir()
            executable = folder / 'fake_phyluce.py'
            executable.write_text(
                '#!/usr/bin/env python3\n'
                'import pathlib, sys\n'
                'args = sys.argv\n'
                'source = pathlib.Path(args[args.index("--alignments") + 1])\n'
                'output = pathlib.Path(args[args.index("--output") + 1])\n'
                'assert source.name == "gene-Su(var)3-3"\n'
                'output.mkdir(parents=True)\n'
                '(output / (output.name + ".nexus")).write_text("#NEXUS\\n")\n'
            )
            executable.chmod(0o755)
            output = folder / 'gene_Su(var)3_3'
            run_job((source, output, output.name), 'intron', executable, folder / 'logs')
            self.assertTrue((output / (output.name + '.nexus')).stat().st_size)

    def test_probe_metadata_is_removed_and_match_threshold_is_configurable(self):
        with tempfile.TemporaryDirectory() as directory:
            m8 = Path(directory) / 'matches.m8'
            m8.write_text(
                'uce-5_p1-design:test\tchr1\t98.0\t120\t0\t0\t1\t120\t200\t81\t1e-20\t210\n'
                'uce-5_p2-design:test\tchr1\t99.0\t120\t0\t0\t1\t120\t160\t41\t1e-20\t220\n'
            )
            output = io.StringIO()
            filter_uce_match_m8_file(str(m8), pct_match=99.0, fout=output)
            self.assertEqual(output.getvalue().splitlines()[0].split('\t')[:3],
                             ['uce-5', 'chr1', '41'])
            output = io.StringIO()
            filter_uce_match_m8_file(str(m8), pct_match=98.0, fout=output)
            self.assertEqual(output.getvalue().splitlines()[0].split('\t')[:3],
                             ['uce-5', 'chr1', '81'])

    def test_exon_parent_resolves_gene_without_gene_attribute(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder / 'genes.gff').write_text(
                'chr\tx\tgene\t1\t500\t.\t+\t.\tID=gene-test\n'
                'chr\tx\tmRNA\t1\t500\t.\t+\t.\tID=rna-test;Parent=gene-test\n'
                'chr\tx\texon\t10\t80\t.\t+\t.\tParent=rna-test\n'
            )
            (folder / 'exons.tsv').write_text('uce-1\tchr\t10\t80\tParent=rna-test\n')
            read_df(str(folder / 'exons.tsv'), '_', str(folder / 'out.txt'),
                    'exon', str(folder / 'genes.gff'))
            self.assertEqual((folder / 'out.txt').read_text(),
                             'uce-1.nexus\texon-rna-test-10-80\tgene-test\n')

    def test_mapped_unclassified_locus_is_not_intergenic(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            copied = folder / 'copied'
            copied.mkdir()
            (copied / 'uce-200.nexus').write_text('#NEXUS\n')
            (copied / 'uce-999.nexus').write_text('#NEXUS\n')
            locations = folder / 'locations.tsv'
            locations.write_text('uce-200\tchr1\t123\n')
            move_remainder(str(copied), str(folder), str(locations))
            self.assertTrue((folder / 'unclassified-mapped/uce-200.nexus').exists())
            self.assertTrue((folder / 'unassigned-regions/uce-999.nexus').exists())
            self.assertFalse((folder / 'intergenic-regions/uce-200.nexus').exists())


if __name__ == '__main__':
    unittest.main()
