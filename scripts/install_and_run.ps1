<#
Simple PowerShell wrapper to run the WSL setup script for this project.

IMPORTANT: Run this script from PowerShell. If WSL is not installed it will attempt to
install it (requires Administrator privileges). After a fresh WSL install you may need
to reboot and run this script again.

This script is intentionally simple and uses an explicit project path. If your project
is located elsewhere, edit the $projectRoot variable below.
#>

Set-StrictMode -Version Latest

function info($m){ Write-Host "[info] $m" -ForegroundColor Cyan }
function warn($m){ Write-Host "[warn] $m" -ForegroundColor Yellow }
function err($m){ Write-Host "[error] $m" -ForegroundColor Red }

# Edit this if your project is in a different path
$projectRoot = 'C:/Users/iqure/Desktop/bio_informatics'

info "Using project root: $projectRoot"

try {
    wsl --status > $null 2>&1
    $wslInstalled = $true
} catch {
    $wslInstalled = $false
}

if (-not $wslInstalled) {
    info "WSL not detected. Attempting to install WSL (requires Admin)..."
    try {
        Start-Process -FilePath wsl -ArgumentList "--install -d Ubuntu" -Verb RunAs -Wait
        info "Requested WSL install. After installation you may need to reboot and set up an Ubuntu user. Re-run this script afterward."
        exit 0
    } catch {
        err "Automatic WSL installation failed. Please run 'wsl --install -d Ubuntu' in an elevated PowerShell and re-run this script."
        exit 1
    }
}

info "WSL is installed. Invoking WSL setup script (this may take 10-30 minutes)."

# Run the WSL setup script inside WSL. The script path is relative to the project root.
$wslCommand = "cd /mnt/c/Users/iqure/Desktop/bio_informatics && bash scripts/wsl_setup.sh"

info "Running: wsl bash -lc \"$wslCommand\""
try {
    wsl bash -lc "$wslCommand"
    info "WSL setup script finished. Check results/benchmark_results.csv in the project root after the script finishes."
} catch {
    err "WSL setup script failed. Check WSL output above."
    exit 1
}
