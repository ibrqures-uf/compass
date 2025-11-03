import os
import time
import logging
from pathlib import Path
from typing import Dict, List, Any
from utils.resource_monitor import ResourceMonitor
from utils.format_converter import FormatConverter
from scoring.baliscore import BaliScorer

class BenchmarkRunner:
    def __init__(self):
        self.logger = logging.getLogger("BenchmarkRunner")
        self.resource_monitor = ResourceMonitor()
        self.format_converter = FormatConverter()
        self.scorer = BaliScorer()
        
        # Initialize paths
        self.results_dir = Path("results")
        self.alignments_dir = self.results_dir / "alignments"
        self.balibase_dir = Path("data/balibase/bb3_release")
        
        # Create necessary directories
        self.alignments_dir.mkdir(parents=True, exist_ok=True)
        
        # Available tools
        self.tools = {
            "mafft": self._run_mafft,
            "muscle": self._run_muscle,
            "clustalo": self._run_clustalo,
            "tcoffee": self._run_tcoffee,
            "probcons": self._run_probcons
        }
    
    def _run_mafft(self, input_file: Path) -> Dict[str, Any]:
        """Run MAFFT alignment"""
        output_file = self.alignments_dir / f"{input_file.stem}_mafft.fasta"
        command = f"mafft --auto {input_file} > {output_file}"
        
        start_time = time.time()
        self.resource_monitor.start()
        
        success = os.system(command) == 0
        
        runtime = time.time() - start_time
        peak_memory = self.resource_monitor.stop()
        
        return {
            "success": success,
            "output": output_file if success else None,
            "stats": {
                "runtime": runtime,
                "peak_memory_mb": peak_memory
            }
        }
    
    def _run_muscle(self, input_file: Path) -> Dict[str, Any]:
        """Run MUSCLE alignment"""
        output_file = self.alignments_dir / f"{input_file.stem}_muscle.fasta"
        command = f"muscle -in {input_file} -out {output_file}"
        
        start_time = time.time()
        self.resource_monitor.start()
        
        success = os.system(command) == 0
        
        runtime = time.time() - start_time
        peak_memory = self.resource_monitor.stop()
        
        return {
            "success": success,
            "output": output_file if success else None,
            "stats": {
                "runtime": runtime,
                "peak_memory_mb": peak_memory
            }
        }
    
    def _run_clustalo(self, input_file: Path) -> Dict[str, Any]:
        """Run Clustal Omega alignment"""
        output_file = self.alignments_dir / f"{input_file.stem}_clustalo.fasta"
        command = f"clustalo -i {input_file} -o {output_file}"
        
        start_time = time.time()
        self.resource_monitor.start()
        
        success = os.system(command) == 0
        
        runtime = time.time() - start_time
        peak_memory = self.resource_monitor.stop()
        
        return {
            "success": success,
            "output": output_file if success else None,
            "stats": {
                "runtime": runtime,
                "peak_memory_mb": peak_memory
            }
        }
    
    def _run_tcoffee(self, input_file: Path) -> Dict[str, Any]:
        """Run T-Coffee alignment"""
        output_file = self.alignments_dir / f"{input_file.stem}_tcoffee.fasta"
        command = f"t_coffee {input_file} -output fasta -outfile {output_file}"
        
        start_time = time.time()
        self.resource_monitor.start()
        
        success = os.system(command) == 0
        
        runtime = time.time() - start_time
        peak_memory = self.resource_monitor.stop()
        
        return {
            "success": success,
            "output": output_file if success else None,
            "stats": {
                "runtime": runtime,
                "peak_memory_mb": peak_memory
            }
        }
    
    def _run_probcons(self, input_file: Path) -> Dict[str, Any]:
        """Run ProbCons alignment"""
        output_file = self.alignments_dir / f"{input_file.stem}_probcons.fasta"
        command = f"probcons {input_file} > {output_file}"
        
        start_time = time.time()
        self.resource_monitor.start()
        
        success = os.system(command) == 0
        
        runtime = time.time() - start_time
        peak_memory = self.resource_monitor.stop()
        
        return {
            "success": success,
            "output": output_file if success else None,
            "stats": {
                "runtime": runtime,
                "peak_memory_mb": peak_memory
            }
        }
    
    def run_all_tools(self, input_file: Path, sequence_id: str) -> Dict[str, Dict[str, Any]]:
        """Run all available MSA tools on a single input file"""
        results = {}
        
        for tool_name, tool_func in self.tools.items():
            try:
                self.logger.info(f"Running {tool_name} on {sequence_id}")
                results[tool_name] = tool_func(input_file)
            except Exception as e:
                self.logger.error(f"Error running {tool_name}: {str(e)}")
                results[tool_name] = {
                    "success": False,
                    "output": None,
                    "stats": {
                        "runtime": 0,
                        "peak_memory_mb": 0
                    }
                }
        
        return results
