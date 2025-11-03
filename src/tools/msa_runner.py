import subprocess
import os
import time
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

class MSARunner:
    def __init__(self, output_dir="results/alignments"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def run_mafft(self, input_fasta, output_path):
        """Run MAFFT alignment"""
        # Use process redirection to avoid shell redirection which can be problematic in Docker
        cmd = ['mafft', '--maxiterate', '1000', '--globalpair', '--reorder', '--quiet', str(input_fasta)]
        stats = {}
        
        try:
            start_time = time.time()
            with open(output_path, 'w') as outfile:
                process = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
                stats['runtime'] = time.time() - start_time
                stats['return_code'] = process.returncode
                stats['peak_memory_mb'] = 0  # Hard to measure peak memory without psutil in Docker
        except subprocess.CalledProcessError as e:
            raise Exception(f"MAFFT failed with error: {e.stderr}")
            
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("MAFFT produced no output")
            
        return stats
    
    def run_muscle(self, input_fasta, output_path):
        """Run MUSCLE alignment"""
        cmd = ['muscle', '-in', str(input_fasta), '-out', str(output_path)]
        stats = {}
        
        try:
            start_time = time.time()
            process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            stats['runtime'] = time.time() - start_time
            stats['return_code'] = process.returncode
            stats['peak_memory_mb'] = 0
        except subprocess.CalledProcessError as e:
            raise Exception(f"MUSCLE failed with error: {e.stderr}")
        
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("MUSCLE produced no output")
        
        return stats
    
    def run_clustalo(self, input_fasta, output_path):
        """Run Clustal Omega"""
        cmd = ['clustalo', '-i', str(input_fasta), '-o', str(output_path), '--force']
        stats = {}
        
        try:
            start_time = time.time()
            process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            stats['runtime'] = time.time() - start_time
            stats['return_code'] = process.returncode
            stats['peak_memory_mb'] = 0
        except subprocess.CalledProcessError as e:
            raise Exception(f"Clustal Omega failed with error: {e.stderr}")
        
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("Clustal Omega produced no output")
        
        return stats
    
    def run_tcoffee(self, input_fasta, output_path):
        """Run T-Coffee"""
        cmd = ['t_coffee', str(input_fasta), '-output', 'fasta_aln', '-outfile', str(output_path), '-quiet']
        stats = {}
        
        try:
            start_time = time.time()
            process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            stats['runtime'] = time.time() - start_time
            stats['return_code'] = process.returncode
            stats['peak_memory_mb'] = 0
        except subprocess.CalledProcessError as e:
            raise Exception(f"T-Coffee failed with error: {e.stderr}")
        
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("T-Coffee produced no output")
        
        return stats
    
    def run_probcons(self, input_fasta, output_path):
        """Run ProbCons"""
        cmd = ['probcons', str(input_fasta)]
        stats = {}
        
        try:
            start_time = time.time()
            with open(output_path, 'w') as outfile:
                process = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
                stats['runtime'] = time.time() - start_time
                stats['return_code'] = process.returncode
                stats['peak_memory_mb'] = 0
        except subprocess.CalledProcessError as e:
            raise Exception(f"ProbCons failed with error: {e.stderr}")
        
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("ProbCons produced no output")
        
        return stats
    
    def run_all_tools(self, input_fasta, tool_name_prefix):
        """Run all MSA tools on input"""
        results = {}
        
        tools = {
            'mafft': self.run_mafft,
            'muscle': self.run_muscle,
            'clustalo': self.run_clustalo
            # Not installed in Docker image
            #'tcoffee': self.run_tcoffee,
            #'probcons': self.run_probcons
        }
        
        input_fasta = str(input_fasta)
        for tool, func in tools.items():
            output = self.output_dir / f"{tool_name_prefix}_{tool}.fasta"
            output = str(output)
            try:
                stats = func(input_fasta, output)
                results[tool] = {
                    'output': output,
                    'stats': stats,
                    'success': True
                }
            except Exception as e:
                results[tool] = {
                    'success': False, 
                    'error': str(e),
                    'output': output
                }
        
        return results