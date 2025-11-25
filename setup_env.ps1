param(
    [string]$venvName = "venv",
    [string]$pythonExe = "python"
)

Write-Host "Creating virtual environment '$venvName' using $pythonExe"
$pythonExe -m venv $venvName
if ($LASTEXITCODE -ne 0) { Write-Error "Failed to create venv"; exit 1 }

Write-Host "Activating virtual environment..."
Write-Host "Run:`n  .\\$venvName\\Scripts\\Activate.ps1` to activate in this session"

Write-Host "Upgrading pip..."
& $pythonExe -m pip install --upgrade pip
Write-Host "Installing requirements..."
& $pythonExe -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { Write-Error "pip install failed"; exit 1 }
Write-Host "Environment ready. Activate with: .\$venvName\Scripts\Activate.ps1"
