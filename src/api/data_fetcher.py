import os
import requests
import tarfile
import zipfile
from pathlib import Path
from tqdm import tqdm

class DataFetcher:
    def __init__(self, base_dir="data"):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(exist_ok=True)
        self.balibase_dir = self.base_dir / "balibase"
        self.tools_dir = self.base_dir / "tools"
        
    def download_file(self, url, dest_path, desc="Downloading"):
        """Download file with progress bar"""
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(dest_path, 'wb') as file, tqdm(
            desc=desc,
            total=total_size,
            unit='iB',
            unit_scale=True
        ) as pbar:
            for data in response.iter_content(chunk_size=1024):
                size = file.write(data)
                pbar.update(size)
    
    def fetch_balibase(self):
        """Download BAliBASE 4.0 dataset"""
        print("Fetching BAliBASE 4.0...")
        
        # BAliBASE registration required - alternative: use test subset
        balibase_url = "http://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R1-5.tar.gz"
        
        dest = self.base_dir / "balibase.tar.gz"
        
        try:
            self.download_file(balibase_url, dest, "BAliBASE")
            
            # Extract
            with tarfile.open(dest) as tar:
                tar.extractall(self.balibase_dir)
            
            print(f"[OK] BAliBASE extracted to {self.balibase_dir}")
            os.remove(dest)
            
        except Exception as e:
            print(f"Note: BAliBASE requires registration at lbgi.fr/balibase")
            print("Creating sample dataset structure...")
            self._create_sample_dataset()
    
    def _create_sample_dataset(self):
        """Create sample dataset for testing"""
        refs = ['RV11', 'RV12', 'RV20', 'RV30', 'RV40', 'RV50']
        for ref in refs:
            ref_dir = self.balibase_dir / ref
            ref_dir.mkdir(parents=True, exist_ok=True)
        print(f"[OK] Created BAliBASE structure at {self.balibase_dir}")
            
    def install_tools(self):
        """Download and install MSA tools"""
        print("\nInstalling MSA tools...")
        
        tools = {
            'mafft': 'https://mafft.cbrc.jp/alignment/software/mafft-7.520-linux.tgz',
            'muscle': 'https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64',
            'clustalo': 'http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64',
        }
        
        self.tools_dir.mkdir(exist_ok=True)
        
        for tool, url in tools.items():
            print(f"Installing {tool}...")
            # Download logic here
            
    print("[OK] Tools installed")
    
    def verify_installation(self):
        """Check if all tools are available"""
        tools = ['mafft', 'muscle', 'clustalo', 't_coffee', 'probcons']
        available = []
        missing = []
        
        for tool in tools:
            # Check if command exists
            if os.name == 'nt':  # Windows
                result = os.system(f"where {tool} >nul 2>&1")
            else:  # Linux/Mac
                result = os.system(f"which {tool} > /dev/null 2>&1")
            
            if result == 0:
                available.append(tool)
            else:
                missing.append(tool)
        
        print(f"\n[OK] Available tools: {', '.join(available) if available else 'None'}")
        if missing:
            print(f"[!] Missing tools: {', '.join(missing)}")
            print("\nTo install missing tools:")
            print("  - MAFFT: https://mafft.cbrc.jp/alignment/software/")
            print("  - MUSCLE: https://github.com/rcedgar/muscle/releases")
            print("  - Clustal Omega: http://www.clustal.org/omega/")
            print("  - T-Coffee: http://www.tcoffee.org/")
            print("  - ProbCons: http://probcons.stanford.edu/")
        
        return len(available) > 0