import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path

class ResultVisualizer:
    def __init__(self, output_dir="results/figures"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        sns.set_style("whitegrid")
        
    def plot_accuracy_comparison(self, df):
        """Plot SP and TC scores for all tools"""
        if df.empty:
            print("  [!] No data to plot for accuracy comparison")
            return
            
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # SP scores
        sns.boxplot(data=df, x='tool', y='sp_score', ax=axes[0])
        axes[0].set_title('Sum-of-Pairs Score by Tool')
        axes[0].set_ylabel('SP Score')
        axes[0].tick_params(axis='x', rotation=45)
        
        # TC scores
        sns.boxplot(data=df, x='tool', y='tc_score', ax=axes[1])
        axes[1].set_title('Total Column Score by Tool')
        axes[1].set_ylabel('TC Score')
        axes[1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'accuracy_comparison.png', dpi=300)
        plt.close()
        print(f"  [OK] Saved accuracy_comparison.png")
    
    def plot_efficiency_comparison(self, df):
        """Plot runtime and memory usage"""
        if df.empty:
            print("  [!] No data to plot for efficiency comparison")
            return
            
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Runtime
        sns.boxplot(data=df, x='tool', y='runtime_sec', ax=axes[0])
        axes[0].set_title('Runtime by Tool')
        axes[0].set_ylabel('Runtime (seconds)')
        axes[0].tick_params(axis='x', rotation=45)
        
        # Memory
        sns.boxplot(data=df, x='tool', y='memory_mb', ax=axes[1])
        axes[1].set_title('Peak Memory Usage by Tool')
        axes[1].set_ylabel('Memory (MB)')
        axes[1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'efficiency_comparison.png', dpi=300)
        plt.close()
        print(f"  [OK] Saved efficiency_comparison.png")
    
    def plot_by_reference_set(self, df):
        """Plot performance across reference sets"""
        if df.empty:
            print("  [!] No data to plot for reference set comparison")
            return
            
        fig, ax = plt.subplots(figsize=(12, 6))
        
        sns.lineplot(data=df, x='reference_set', y='sp_score', 
                     hue='tool', marker='o', ax=ax)
        ax.set_title('SP Score Across BAliBASE Reference Sets')
        ax.set_ylabel('SP Score')
        ax.set_xlabel('Reference Set')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'performance_by_refset.png', dpi=300)
        plt.close()
        print(f"  [OK] Saved performance_by_refset.png")