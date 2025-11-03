import shutil
import os
from pathlib import Path

class ToolFinder:
    @staticmethod
    def find_msa_tools():
        """Find installed MSA tools and return their paths"""
        tools = {
            'mafft': shutil.which('mafft'),
            'muscle': shutil.which('muscle'),
            'clustalo': shutil.which('clustalo'),
            't_coffee': shutil.which('t_coffee'),
            'probcons': shutil.which('probcons')
        }
        return {name: path for name, path in tools.items() if path is not None}

    @staticmethod
    def verify_tool(tool_name):
        """Check if a specific tool is available"""
        return shutil.which(tool_name) is not None

    @staticmethod
    def get_tool_path(tool_name):
        """Get the full path of a tool if available"""
        return shutil.which(tool_name)
