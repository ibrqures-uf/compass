import sys
from pathlib import Path

# Add current directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

import argparse
import pandas as pd
from tqdm import tqdm
from api.data_fetcher import DataFetcher
from tools.msa_runner import MSARunner
from scoring.baliscore import BaliScorer
from analysis.aggregator import ResultAggregator
from analysis.visualizer import ResultVisualizer

def main():
    print("=" * 60)
    print("MSA Benchmarking Pipeline")
    print("=" * 60)
    
    # Step 1: Fetch data and tools
    print("\n[1] Fetching data and tools...")
    fetcher = DataFetcher()
    fetcher.fetch_balibase()
    fetcher.install_tools()
    
    if not fetcher.verify_installation():
        print("\nWarning: Some tools are missing. Install manually and re-run.")
        print("Continuing with available tools...")
    
    # Step 2: Run alignments
    print("\n[2] Running MSA tools on BAliBASE...")
    runner = MSARunner()
    scorer = BaliScorer()
    
    results = []
    balibase_dir = Path("data/balibase/bb3_release")
    
    # Count total files first
    all_files = []
    for ref_set in ['RV11', 'RV12', 'RV20', 'RV30', 'RV40', 'RV50']:
        ref_dir = balibase_dir / ref_set
        if ref_dir.exists():
            all_files.extend(list(ref_dir.glob("*.tfa")))
    
    total_files = len(all_files)
    if total_files == 0:
        print("\n[!] No BAliBASE data found!")
        print("Please download BAliBASE dataset from: http://www.lbgi.fr/balibase")
        print(f"Extract it to: {balibase_dir}")
        return
    
    # Optionally limit files for a quick smoke test via BENCH_LIMIT env var
    import os
    bench_limit = int(os.getenv('BENCH_LIMIT', '0'))
    allowed_names = None
    if bench_limit > 0:
        all_files = all_files[:bench_limit]
        allowed_names = set(p.name for p in all_files)
        total_files = len(all_files)

    print(f"\nFound {total_files} sequence files to process")
    print("="*60)
    
    # Process with progress bar
    with tqdm(total=total_files, desc="Overall Progress", unit="file") as pbar:
        for ref_set in ['RV11', 'RV12', 'RV20', 'RV30', 'RV40', 'RV50']:
            ref_dir = balibase_dir / ref_set
            if not ref_dir.exists():
                continue
            
            files = list(ref_dir.glob("*.tfa"))
            if not files:
                continue
            
            tqdm.write(f"\n{'='*60}")
            tqdm.write(f"Processing {ref_set} ({len(files)} files)")
            tqdm.write(f"{'='*60}")
            
            for input_file in files:
                # If BENCH_LIMIT was set, skip files not in the allowed set
                if allowed_names is not None and input_file.name not in allowed_names:
                    continue
                tqdm.write(f"\n  File: {input_file.name}")
                
                # Run all tools
                tool_results = runner.run_all_tools(input_file, input_file.stem)

                # Locate the reference alignment file for scoring.
                # BAliBASE includes reference alignments in .msf (MSF) or .rsf files;
                # prefer .msf, then .rsf, then fall back to the .tfa if nothing else.
                ref_alignment = input_file.with_suffix('.msf')
                if not ref_alignment.exists():
                    ref_alignment = input_file.with_suffix('.rsf')
                if not ref_alignment.exists():
                    # fallback (may be unaligned) but keep previous behavior
                    ref_alignment = input_file

                # Score each tool
                for tool, result in tool_results.items():
                    if not result['success']:
                        tqdm.write(f"    [X] {tool} failed")
                        continue
                    
                    try:
                        # Calculate scores
                        scores = scorer.score_alignment(
                            result['output'],
                            ref_alignment
                        )
                        
                        results.append({
                            'reference_set': ref_set,
                            'sequence': input_file.stem,
                            'tool': tool,
                            'sp_score': scores['sp_score'],
                            'tc_score': scores['tc_score'],
                            'runtime_sec': result['stats']['runtime'],
                            'memory_mb': result['stats']['peak_memory_mb']
                        })
                        
                        tqdm.write(f"    [OK] {tool}: SP={scores['sp_score']:.3f}, TC={scores['tc_score']:.3f}, Time={result['stats']['runtime']:.2f}s")
                        
                    except Exception as e:
                        tqdm.write(f"    [X] {tool} scoring failed: {e}")
                
                pbar.update(1)
    
    # Step 3: Aggregate results
    print("\n" + "="*60)
    print("[3] Aggregating results...")
    df = pd.DataFrame(results)
    
    # Check if we have any results
    if df.empty:
        print("[!] No results generated")
        print("\nPlease install at least one MSA tool:")
        print("  - MAFFT: https://mafft.cbrc.jp/alignment/software/")
        print("  - MUSCLE: https://github.com/rcedgar/muscle/releases")
        print("  - Clustal Omega: http://www.clustal.org/omega/")
        print("\nExiting...")
        return
    
    df.to_csv("results/benchmark_results.csv", index=False)
    print(f"  [OK] Saved {len(results)} results to results/benchmark_results.csv")
    print(f"  [OK] Processed {len(df['sequence'].unique())} sequences")
    print(f"  [OK] Used tools: {', '.join(df['tool'].unique())}")
    
    # Step 4: Generate visualizations
    print("\n[4] Generating visualizations...")
    viz = ResultVisualizer()
    viz.plot_accuracy_comparison(df)
    viz.plot_efficiency_comparison(df)
    viz.plot_by_reference_set(df)
    print("  [OK] Plots saved to results/figures/")
    
    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    summary = df.groupby('tool').agg({
        'sp_score': ['mean', 'std'],
        'tc_score': ['mean', 'std'],
        'runtime_sec': ['mean', 'sum'],
        'memory_mb': ['mean', 'max']
    }).round(3)
    print(summary)
    
    print("\n" + "=" * 60)
    print("Benchmarking complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()