from Bio import SeqIO, AlignIO

class FormatConverter:
    @staticmethod
    def to_fasta(input_file, output_file, input_format='clustal'):
        try:
            alignment = AlignIO.read(input_file, input_format)
            AlignIO.write(alignment, output_file, 'fasta')
            return True
        except:
            return False
    
    @staticmethod
    def validate_fasta(fasta_file):
        try:
            sequences = list(SeqIO.parse(fasta_file, 'fasta'))
            return len(sequences) > 0
        except:
            return False
