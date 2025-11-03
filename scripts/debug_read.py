from Bio import AlignIO
from pathlib import Path
p = Path('data/balibase/bb3_release/RV11/BB11001.tfa')
print('Exists:', p.exists())
print('Size:', p.stat().st_size)
for fmt in ['fasta','clustal','stockholm','phylip','maf','msf']:
    try:
        print('Trying',fmt)
        aln = AlignIO.read(str(p), fmt)
        print('OK',fmt,'records',len(aln))
    except Exception as e:
        print('Error',fmt,':',e)
