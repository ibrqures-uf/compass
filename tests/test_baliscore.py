import tempfile
from pathlib import Path
from scoring.baliscore import BaliScorer
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO


def test_simple_sp_tc():
    # Create a tiny reference alignment of 3 sequences
    records_ref = [
        SeqRecord(Seq('ACGT'), id='s1'),
        SeqRecord(Seq('A-GT'), id='s2'),
        SeqRecord(Seq('AC-T'), id='s3')
    ]

    records_test = [
        SeqRecord(Seq('ACGT'), id='s1'),
        SeqRecord(Seq('A-GT'), id='s2'),
        SeqRecord(Seq('AC-T'), id='s3')
    ]

    aln_ref = MultipleSeqAlignment(records_ref)
    aln_test = MultipleSeqAlignment(records_test)

    with tempfile.TemporaryDirectory() as td:
        refp = Path(td) / 'ref.fasta'
        testp = Path(td) / 'test.fasta'
        AlignIO.write(aln_ref, str(refp), 'fasta')
        AlignIO.write(aln_test, str(testp), 'fasta')

        scorer = BaliScorer()
        res = scorer.score_alignment(str(testp), str(refp))
        assert 'sp_score' in res and 'tc_score' in res
        assert 0.0 <= res['sp_score'] <= 1.0
        assert 0.0 <= res['tc_score'] <= 1.0
