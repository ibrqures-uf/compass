import os
import sys
import yaml
import logging
import requests
import platform
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass

@dataclass
class ToolConfig:
    version: str
    url: str
    command: str

@dataclass
class DatasetConfig:
    version: str
    url: str
    refs: List[str]

class MSAToolManager:
    def __init__(self, config_path: str = "config.yaml"):
        self.config = self._load_config(config_path)
        self.tools_dir = Path("tools")
        self.data_dir = Path("data")
        self.setup_logging()
        
    def _load_config(self, config_path: str) -> dict:
        """Load configuration from YAML file"""
        with open(config_path) as f:
            return yaml.safe_load(f)
    
    def setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger("MSAToolManager")
    
    def get_system_info(self) -> Dict[str, str]:
        """Get system information"""
        return {
            "os": platform.system(),
            "architecture": platform.machine(),
            "python_version": sys.version,
            "platform": platform.platform()
        }
    
    def check_tool_availability(self, tool_name: str) -> bool:
        """Check if a tool is available in system PATH"""
        tool_cmd = self.config["tools"][tool_name]["command"]
        try:
            subprocess.run([tool_cmd, "--version"], 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE)
            return True
        except FileNotFoundError:
            return False
    
    def install_tool(self, tool_name: str) -> bool:
        """Install MSA tool"""
        tool_config = self.config["tools"][tool_name]
        self.logger.info(f"Installing {tool_name} version {tool_config['version']}")
        
        # Create tools directory
        self.tools_dir.mkdir(exist_ok=True)
        
        # Download tool
        response = requests.get(tool_config["url"], stream=True)
        if response.status_code != 200:
            self.logger.error(f"Failed to download {tool_name}")
            return False
            
        # Save and process based on tool
        tool_path = self.tools_dir / f"{tool_name}_v{tool_config['version']}"
        with open(tool_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                
        # Make executable
        os.chmod(tool_path, 0o755)
        
        return True
    
    def setup_tools(self) -> Dict[str, bool]:
        """Setup all required MSA tools"""
        results = {}
        for tool_name in self.config["tools"]:
            if not self.check_tool_availability(tool_name):
                results[tool_name] = self.install_tool(tool_name)
            else:
                results[tool_name] = True
        return results
    
    def fetch_dataset(self, dataset_name: str = "balibase") -> bool:
        """Fetch and setup dataset"""
        dataset_config = self.config["datasets"][dataset_name]
        self.logger.info(f"Fetching {dataset_name} version {dataset_config['version']}")
        
        dataset_dir = self.data_dir / dataset_name
        dataset_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Download dataset
            response = requests.get(dataset_config["url"], stream=True)
            if response.status_code != 200:
                self.logger.error(f"Failed to download {dataset_name}")
                return False
                
            # Save dataset
            dataset_file = dataset_dir / f"{dataset_name}.tar.gz"
            with open(dataset_file, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
                    
            # Extract dataset
            subprocess.run(["tar", "-xzf", str(dataset_file), "-C", str(dataset_dir)])
            dataset_file.unlink()  # Remove the tar.gz file
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error fetching dataset: {str(e)}")
            return False
    
    def verify_setup(self) -> Dict[str, bool]:
        """Verify all tools and datasets are properly set up"""
        verification = {
            "tools": {},
            "datasets": {},
            "directories": {}
        }
        
        # Check tools
        for tool_name in self.config["tools"]:
            verification["tools"][tool_name] = self.check_tool_availability(tool_name)
        
        # Check datasets
        for dataset_name in self.config["datasets"]:
            dataset_dir = self.data_dir / dataset_name
            verification["datasets"][dataset_name] = dataset_dir.exists()
        
        # Check directories
        for dir_name in ["results_dir", "alignments_dir", "figures_dir"]:
            dir_path = Path(self.config["output"][dir_name])
            verification["directories"][dir_name] = dir_path.exists()
        
        return verification
