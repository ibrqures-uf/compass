from pathlib import Path
from tools.msa_runner import MSARunner
from utils.tool_finder import find_executable, suggest_install_command
from tools.benchmark_runner import discover_families

def main():
    data_dir = Path('data/balibase/bb3_release')
    families = discover_families(data_dir)
    if not families:
        print('No families found under', data_dir)
        return

    fam = families[0]
    print('Selected family:', fam)

    binary_map = {
        'mafft': 'mafft',
        'muscle': 'muscle',
        'clustalo': 'clustalo',
        'tcoffee': 't_coffee',
        'probcons': 'probcons'
    }

    available = []
    for t, exe in binary_map.items():
        if find_executable(exe):
            available.append(t)

    if not available:
        print('No MSA tools found on PATH. Detected none of:', list(binary_map.values()))
        print('Install tools (recommended via WSL + conda/bioconda) or run locally with installed binaries.')
        return

    print('Available tools:', available)
    # Respect WSL if the environment wants to use it (detect via env var or default False)
    use_wsl = False
    runner = MSARunner(use_wsl=use_wsl)
    for tool in available:
        out = Path('results/alignments') / fam.parent.name
        out.mkdir(parents=True, exist_ok=True)
        out_path = out / f"{fam.stem}_{tool}.fasta"
        print(f"Running {tool} on {fam} -> {out_path}")
        res = runner.run_tool(tool, str(fam), str(out_path))
        print('Result:', res)

if __name__ == '__main__':
    main()
