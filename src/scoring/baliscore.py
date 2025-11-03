from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np

class BaliScorer:
    def __init__(self):
        pass
    
    def calculate_sp_score(self, test_aln, ref_aln):
        """Calculate Sum-of-Pairs score"""
        test_pairs = self._get_alignment_pairs(test_aln)
        ref_pairs = self._get_alignment_pairs(ref_aln)
        
        correct = len(test_pairs & ref_pairs)
        total = len(ref_pairs)
        
        sp_score = correct / total if total > 0 else 0
        return sp_score
    
    def calculate_tc_score(self, test_aln, ref_aln):
        """Calculate Total Column score"""
        test_cols = self._get_columns(test_aln)
        ref_cols = self._get_columns(ref_aln)
        
        correct_cols = sum(1 for tc, rc in zip(test_cols, ref_cols) if tc == rc)
        total_cols = len(ref_cols)
        
        tc_score = correct_cols / total_cols if total_cols > 0 else 0
        return tc_score
    
    def _get_alignment_pairs(self, alignment):
        """Extract all aligned residue pairs"""
        pairs = set()
        seqs = [str(record.seq) for record in alignment]
        
        for i in range(len(alignment[0])):
            col = [seq[i] for seq in seqs]
            for j in range(len(col)):
                for k in range(j + 1, len(col)):
                    if col[j] != '-' and col[k] != '-':
                        pairs.add((j, k, i))
        
        return pairs
    
    def _get_columns(self, alignment):
        """Get column representations"""
        columns = []
        for i in range(len(alignment[0])):
            col = tuple(str(record.seq)[i] for record in alignment)
            columns.append(col)
        return columns
    
    def score_alignment(self, test_file, ref_file):
        """Score test alignment against reference"""
        # Try to read alignments using Biopython with format auto-detection
        def _read_try(path):
            for fmt in ('fasta', 'msf', 'clustal', 'stockholm', 'phylip', 'rsf'):
                try:
                    return AlignIO.read(str(path), fmt)
                except Exception:
                    continue
            raise ValueError(f"Unable to read alignment: {path}")

        # For test alignments, try multiple formats
        test_aln = _read_try(test_file)
        
        # For reference alignments, prefer MSF format first
        try:
            ref_aln = AlignIO.read(str(ref_file), "msf")
        except Exception:
            ref_aln = _read_try(ref_file)

        # Extract just the aligned sequences without gaps for comparison
        def strip_gaps(seq):
            return str(seq).replace('-', '')
            
        # Validate that sequences present in both alignments match their ungapped content
        ref_seqs = {r.id: strip_gaps(r.seq) for r in ref_aln}
        test_seqs = {r.id: strip_gaps(r.seq) for r in test_aln}
        
        common = set(ref_seqs.keys()) & set(test_seqs.keys())
        for cid in common:
            ref_ungapped = ref_seqs[cid]
            test_ungapped = test_seqs[cid]
            if ref_ungapped != test_ungapped:
                raise ValueError(f"Sequence content differs for ID '{cid}'")

        sp = self.calculate_sp_score(test_aln, ref_aln)
        tc = self.calculate_tc_score(test_aln, ref_aln)

        return {'sp_score': float(sp), 'tc_score': float(tc)}