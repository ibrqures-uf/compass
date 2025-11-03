import os
from pathlib import Path

def create_file(path, content):
    """Create file with content"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f"[OK] Created {path}")

def setup_complete_project():
    """Auto-generate entire project structure with code"""
    
    print("Setting up MSA Benchmark project...\n")
    
    # Create directories
    Path('data').mkdir(exist_ok=True)
    Path('results/alignments').mkdir(parents=True, exist_ok=True)
    Path('results/figures').mkdir(parents=True, exist_ok=True)
    print("[OK] Created directories\n")
    
    # utils/resource_monitor.py
    resource_monitor_code = '''import time
import psutil
import subprocess

class ResourceMonitor:
    def run_with_monitoring(self, command):
        start_time = time.time()
        initial_memory = psutil.virtual_memory().used / (1024 ** 2)
        
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        peak_memory = initial_memory
        try:
            ps_process = psutil.Process(process.pid)
            while process.poll() is None:
                try:
                    mem_info = ps_process.memory_info()
                    current_mem = mem_info.rss / (1024 ** 2)
                    peak_memory = max(peak_memory, current_mem)
                except:
                    break
                time.sleep(0.1)
        except:
            peak_memory = 0
        
        stdout, stderr = process.communicate()
        end_time = time.time()
        
        return {
            'runtime': end_time - start_time,
            'peak_memory_mb': max(0, peak_memory - initial_memory),
            'return_code': process.returncode
        }
'''
    create_file(Path('utils/resource_monitor.py'), resource_monitor_code)
    
    # analysis/aggregator.py
    aggregator_code = '''import pandas as pd

class ResultAggregator:
    def __init__(self):
        self.results = []
    
    def add_result(self, result):
        self.results.append(result)
    
    def to_dataframe(self):
        return pd.DataFrame(self.results)
    
    def save_csv(self, filepath):
        df = self.to_dataframe()
        df.to_csv(filepath, index=False)
        print(f"Results saved to {filepath}")
'''
    create_file(Path('analysis/aggregator.py'), aggregator_code)
    
    # utils/format_converter.py
    converter_code = '''from Bio import SeqIO, AlignIO

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
'''
    create_file(Path('utils/format_converter.py'), converter_code)
    
    print("\n" + "="*60)
    print("Project setup complete!")
    print("="*60)
    print("\nNext steps:")
    print("1. pip install -r requirements.txt")
    print("2. python main.py")

if __name__ == "__main__":
    setup_complete_project()