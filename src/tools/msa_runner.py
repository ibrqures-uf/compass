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
        
    def _run_with_memory(self, cmd, output_path=None, use_stdout_file=False):
        """Run a command with /usr/bin/time -v to capture peak memory usage."""
        import re
        time_cmd = ['/usr/bin/time', '-v'] + cmd
        stats = {}
        try:
            start_time = time.time()
            if use_stdout_file and output_path:
                with open(output_path, 'w') as outfile:
                    process = subprocess.run(time_cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
            else:
                process = subprocess.run(time_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            stats['runtime'] = time.time() - start_time
            stats['return_code'] = process.returncode
            # Parse peak memory from stderr
            mem_match = re.search(r'Maximum resident set size \(kbytes\): (\d+)', process.stderr)
            if mem_match:
                stats['peak_memory_mb'] = round(int(mem_match.group(1)) / 1024, 2)
            else:
                stats['peak_memory_mb'] = None
        except subprocess.CalledProcessError as e:
            raise Exception(f"Command failed with error: {e.stderr}")
        return stats, process

    def run_mafft(self, input_fasta, output_path):
        """Run MAFFT alignment"""
        cmd = ['mafft', '--maxiterate', '1000', '--globalpair', '--reorder', '--quiet', str(input_fasta)]
        stats, _ = self._run_with_memory(cmd, output_path, use_stdout_file=True)
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("MAFFT produced no output")
        return stats
    
    def run_muscle(self, input_fasta, output_path):
        """Run MUSCLE alignment"""
        cmd = ['muscle', '-in', str(input_fasta), '-out', str(output_path)]
        stats, _ = self._run_with_memory(cmd)
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("MUSCLE produced no output")
        return stats
    
    def run_clustalo(self, input_fasta, output_path):
        """Run Clustal Omega"""
        cmd = ['clustalo', '-i', str(input_fasta), '-o', str(output_path), '--force']
        stats, _ = self._run_with_memory(cmd)
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("Clustal Omega produced no output")
        return stats
    
    def run_tcoffee(self, input_fasta, output_path):
        """Run T-Coffee"""
        cmd = ['t_coffee', str(input_fasta), '-output', 'fasta_aln', '-outfile', str(output_path), '-quiet']
        stats, _ = self._run_with_memory(cmd)
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("T-Coffee produced no output")
        return stats
    
    def run_probcons(self, input_fasta, output_path):
        """Run ProbCons"""
        cmd = ['probcons', str(input_fasta)]
        stats, _ = self._run_with_memory(cmd, output_path, use_stdout_file=True)
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise Exception("ProbCons produced no output")
        return stats
    
    def run_all_tools(self, input_fasta, tool_name_prefix):
        """Run all MSA tools on input"""
        results = {}
        
        tools = {
            'mafft': self.run_mafft,
            'muscle': self.run_muscle,
            'clustalo': self.run_clustalo,
            'tcoffee': self.run_tcoffee,
            'probcons': self.run_probcons
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