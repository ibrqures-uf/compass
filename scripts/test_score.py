from scoring.baliscore import BaliScorer

if __name__ == '__main__':
    sc = BaliScorer()
    test_file = 'results/alignments/BB11001_mafft.fasta'
    # prefer MSF/RSF reference alignments from BAliBASE
    import os
    ref_msf = 'data/balibase/bb3_release/RV11/BB11001.msf'
    ref_rsf = 'data/balibase/bb3_release/RV11/BB11001.rsf'
    if os.path.exists(ref_msf):
        ref = ref_msf
    elif os.path.exists(ref_rsf):
        ref = ref_rsf
    else:
        ref = 'data/balibase/bb3_release/RV11/BB11001.tfa'

    res = sc.score_alignment(test_file, ref)
    print(res)
