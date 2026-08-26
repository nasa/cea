#Requires -Version 5.1

<#
.SYNOPSIS
Build and install CEA and its test dependencies on 64-bit Windows with Intel oneAPI.

.DESCRIPTION
The script performs a complete native and Python build:

* initializes the Intel oneAPI environment when necessary;
* installs Python build/test prerequisites into the selected Python environment;
* downloads, builds, and installs pFUnit;
* builds CEA with the Fortran, C, Python, and MATLAB-compatible interfaces;
* runs the CTest suite, installs the Python package, and runs pytest.

Excel and C++ interfaces are disabled. Existing build directories are reused so
that rerunning the script is incremental.

.PARAMETER Python
Python command or absolute executable path. Defaults to "python".

.PARAMETER BuildDir
Native build and dependency directory. Defaults to build-win-oneapi in the CEA
source tree.

.PARAMETER InstallPrefix
CEA installation directory. Defaults to <BuildDir>\install.

.PARAMETER OneApiSetVars
Optional absolute path to Intel oneAPI setvars.bat or oneapi-vars.bat. Both the
component and unified standard installation layouts are searched when omitted.

.PARAMETER M4Executable
Optional absolute path to GNU M4. The active PATH and common MSYS2, Git, and
Conda installation locations are searched when omitted.

.PARAMETER Parallel
Maximum number of parallel build jobs. Defaults to 4.

.PARAMETER SkipPythonPrerequisites
Do not install or update Python build/test prerequisites.

.PARAMETER SkipTests
Build and install CEA without running CTest or pytest. CEA is still configured
with CEA_BUILD_TESTING=ON so the test executables are built.

.EXAMPLE
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\install_windows_oneapi.ps1

.EXAMPLE
.\scripts\install_windows_oneapi.ps1 -Python C:\Python314\python.exe -Parallel 8
#>

[CmdletBinding()]
param(
    [string]$Python = "python",
    [string]$BuildDir = "",
    [string]$InstallPrefix = "",
    [string]$OneApiSetVars = "",
    [string]$M4Executable = "",
    [ValidateRange(1, 128)]
    [int]$Parallel = 4,
    [switch]$SkipPythonPrerequisites,
    [switch]$SkipTests
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = "Stop"

$RepoRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot ".."))
if (-not $BuildDir) {
    $BuildDir = Join-Path $RepoRoot "build-win-oneapi"
}
$BuildDir = [System.IO.Path]::GetFullPath($BuildDir)

if (-not $InstallPrefix) {
    $InstallPrefix = Join-Path $BuildDir "install"
}
$InstallPrefix = [System.IO.Path]::GetFullPath($InstallPrefix)

$PFUnitVersion = "4.19.0"
$PFUnitSourceDir = Join-Path $BuildDir "_deps\pfunit-src"
$PFUnitBuildDir = Join-Path $BuildDir "_deps\pfunit-build"
$PFUnitInstallDir = Join-Path $BuildDir "_deps\pfunit-install"
$CEABuildDir = Join-Path $BuildDir "cea"

function Invoke-External {
    param(
        [Parameter(Mandatory = $true)]
        [string]$FilePath,
        [string[]]$ArgumentList = @(),
        [string]$WorkingDirectory = ""
    )

    if ($WorkingDirectory) {
        Push-Location $WorkingDirectory
    }
    try {
        Write-Host "> $FilePath $($ArgumentList -join ' ')" -ForegroundColor DarkGray
        & $FilePath @ArgumentList
        if ($LASTEXITCODE -ne 0) {
            throw "Command failed with exit code ${LASTEXITCODE}: $FilePath"
        }
    }
    finally {
        if ($WorkingDirectory) {
            Pop-Location
        }
    }
}

function Find-OneApiSetVars {
    param([string]$RequestedPath)

    if ($RequestedPath) {
        if (-not (Test-Path -LiteralPath $RequestedPath -PathType Leaf)) {
            throw "Intel oneAPI environment script was not found: $RequestedPath"
        }
        return [System.IO.Path]::GetFullPath($RequestedPath)
    }

    $InstallRoots = @()
    if (${env:ProgramFiles(x86)}) {
        $InstallRoots += Join-Path ${env:ProgramFiles(x86)} "Intel\oneAPI"
    }
    if ($env:ProgramFiles) {
        $InstallRoots += Join-Path $env:ProgramFiles "Intel\oneAPI"
    }

    foreach ($InstallRoot in $InstallRoots) {
        $Candidate = Join-Path $InstallRoot "setvars.bat"
        if (Test-Path -LiteralPath $Candidate -PathType Leaf) {
            return $Candidate
        }
    }

    # Intel's unified directory layout stores oneapi-vars.bat below a toolkit
    # version directory, for example oneAPI\2025.3\oneapi-vars.bat.
    foreach ($InstallRoot in $InstallRoots) {
        if (-not (Test-Path -LiteralPath $InstallRoot -PathType Container)) {
            continue
        }
        $VersionDirectories = Get-ChildItem -LiteralPath $InstallRoot -Directory |
            Sort-Object -Property Name -Descending
        foreach ($VersionDirectory in $VersionDirectories) {
            $Candidate = Join-Path $VersionDirectory.FullName "oneapi-vars.bat"
            if (Test-Path -LiteralPath $Candidate -PathType Leaf) {
                return $Candidate
            }
        }
    }
    return $null
}

function Import-OneApiEnvironment {
    param([string]$SetVarsPath)

    $Command = "call `"$SetVarsPath`" intel64 >nul && set"
    $EnvironmentLines = & $env:ComSpec /d /s /c $Command
    if ($LASTEXITCODE -ne 0) {
        throw "Intel oneAPI environment setup failed: $SetVarsPath"
    }

    foreach ($Line in $EnvironmentLines) {
        $Separator = $Line.IndexOf("=")
        if ($Separator -le 0) {
            continue
        }
        $Name = $Line.Substring(0, $Separator)
        $Value = $Line.Substring($Separator + 1)
        Set-Item -Path "Env:$Name" -Value $Value
    }
}

function Require-Command {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Name,
        [string]$HelpText = ""
    )

    $Command = Get-Command $Name -ErrorAction SilentlyContinue
    if (-not $Command) {
        if ($HelpText) {
            throw "Required command '$Name' was not found. $HelpText"
        }
        throw "Required command '$Name' was not found."
    }
    Write-Host "Found $Name at $($Command.Source)"
}

function Find-M4Executable {
    param([string]$RequestedPath)

    if ($RequestedPath) {
        if (-not (Test-Path -LiteralPath $RequestedPath -PathType Leaf)) {
            throw "GNU M4 executable was not found: $RequestedPath"
        }
        return [System.IO.Path]::GetFullPath($RequestedPath)
    }

    foreach ($Name in @("gm4", "m4", "m4.exe")) {
        $Command = Get-Command $Name -ErrorAction SilentlyContinue
        if ($Command -and $Command.Source) {
            return $Command.Source
        }
    }

    $Candidates = @(
        "C:\msys64\usr\bin\m4.exe",
        "C:\msys32\usr\bin\m4.exe"
    )
    if ($env:ProgramFiles) {
        $Candidates += Join-Path $env:ProgramFiles "Git\usr\bin\m4.exe"
    }
    if (${env:ProgramFiles(x86)}) {
        $Candidates += Join-Path ${env:ProgramFiles(x86)} "Git\usr\bin\m4.exe"
    }
    if ($env:CONDA_PREFIX) {
        $Candidates += Join-Path $env:CONDA_PREFIX "Library\usr\bin\m4.exe"
    }

    foreach ($Candidate in $Candidates) {
        if (Test-Path -LiteralPath $Candidate -PathType Leaf) {
            return [System.IO.Path]::GetFullPath($Candidate)
        }
    }
    return $null
}

function Disable-GFTLDependencySelfTests {
    param([string]$PFUnitSource)

    $CMakeFile = Join-Path $PFUnitSource (
        "extern\fArgParse\extern\gFTL-shared\extern\gFTL\include\v2\parameters\CMakeLists.txt"
    )
    if (-not (Test-Path -LiteralPath $CMakeFile -PathType Leaf)) {
        throw "Unable to locate the bundled gFTL CMake file: $CMakeFile"
    }

    $Content = [System.IO.File]::ReadAllText($CMakeFile)
    $TestDirective = "add_subdirectory(tests EXCLUDE_FROM_ALL)"
    $GuardedDirective = @(
        "if (ENABLE_TESTS)",
        "  $TestDirective",
        "endif ()"
    ) -join [Environment]::NewLine

    if ($Content.Contains($GuardedDirective)) {
        return
    }
    if (-not $Content.Contains($TestDirective)) {
        throw "The bundled gFTL test configuration has an unexpected format: $CMakeFile"
    }

    # gFTL configures its parameter-generation self-tests unconditionally. On
    # Windows those self-tests require AWK and cpp even when pFUnit was
    # configured with ENABLE_TESTS=OFF. Guard only that dependency test folder;
    # the gFTL libraries and all downstream CEA tests remain enabled.
    $Content = $Content.Replace($TestDirective, $GuardedDirective)
    [System.IO.File]::WriteAllText($CMakeFile, $Content)
    Write-Host "Disabled unused gFTL dependency self-tests in $CMakeFile" -ForegroundColor Cyan
}

if (-not (Test-Path -LiteralPath (Join-Path $RepoRoot "CMakeLists.txt") -PathType Leaf)) {
    throw "CEA source root could not be determined from $PSScriptRoot."
}

Write-Host "CEA source:      $RepoRoot" -ForegroundColor Cyan
Write-Host "Build directory: $BuildDir" -ForegroundColor Cyan
Write-Host "Install prefix:  $InstallPrefix" -ForegroundColor Cyan

$Ifx = Get-Command "ifx" -ErrorAction SilentlyContinue
$Cl = Get-Command "cl" -ErrorAction SilentlyContinue
if (-not $Ifx -or -not $Cl) {
    $SetVars = Find-OneApiSetVars -RequestedPath $OneApiSetVars
    if (-not $SetVars) {
        throw "ifx/cl are unavailable and neither setvars.bat nor oneapi-vars.bat was found in the standard Intel oneAPI locations. Install the oneAPI HPC Toolkit and Visual Studio C++ tools, or pass the environment script's full path with -OneApiSetVars."
    }
    Write-Host "Loading Intel oneAPI environment from $SetVars" -ForegroundColor Cyan
    Import-OneApiEnvironment -SetVarsPath $SetVars
}

Require-Command -Name "ifx" -HelpText "Install the Intel oneAPI HPC Toolkit."
Require-Command -Name "cl" -HelpText "Install the Visual Studio Desktop development with C++ workload."
Require-Command -Name "git" -HelpText "Install Git for Windows."

$PythonExecutable = (& $Python -c "import sys; print(sys.executable)") | Select-Object -Last 1
if ($LASTEXITCODE -ne 0 -or -not $PythonExecutable) {
    throw "Unable to run the requested Python interpreter: $Python"
}
$PythonExecutable = $PythonExecutable.Trim()

$PythonInfoJson = & $PythonExecutable -c "import json, struct, sys, sysconfig; print(json.dumps({'version': list(sys.version_info[:3]), 'bits': struct.calcsize('P') * 8, 'implementation': sys.implementation.name, 'free_threaded': bool(sysconfig.get_config_var('Py_GIL_DISABLED'))}))"
if ($LASTEXITCODE -ne 0) {
    throw "Unable to inspect Python: $PythonExecutable"
}
$PythonInfo = $PythonInfoJson | ConvertFrom-Json
if ($PythonInfo.version[0] -ne 3 -or $PythonInfo.version[1] -lt 11) {
    throw "CEA requires Python 3.11 or newer; found $($PythonInfo.version -join '.')."
}
if ($PythonInfo.bits -ne 64) {
    throw "A 64-bit Python installation is required for this Windows build."
}
if ($PythonInfo.implementation -ne "cpython") {
    Write-Warning "CEA's compiled Python binding is primarily tested with CPython."
}
if ($PythonInfo.free_threaded) {
    Write-Warning "This is a free-threaded Python build; use the standard GIL-enabled CPython build if compilation or runtime tests fail."
}
Write-Host "Using Python $($PythonInfo.version -join '.') at $PythonExecutable" -ForegroundColor Cyan

if (-not $SkipPythonPrerequisites) {
    Invoke-External -FilePath $PythonExecutable -ArgumentList @(
        "-m", "pip", "install", "--upgrade",
        "pip", "cmake>=3.24", "ninja", "setuptools",
        "scikit-build-core>=0.11", "cython", "numpy>=2.4", "pytest"
    )
}

$PythonScriptsDir = (& $PythonExecutable -c "import sysconfig; print(sysconfig.get_path('scripts'))") |
    Select-Object -Last 1
if ($LASTEXITCODE -ne 0 -or -not $PythonScriptsDir) {
    throw "Unable to determine the Python scripts directory."
}
$PythonScriptsDir = $PythonScriptsDir.Trim()
if (Test-Path -LiteralPath $PythonScriptsDir -PathType Container) {
    $env:PATH = "$PythonScriptsDir;$env:PATH"
}

Require-Command -Name "cmake"
Require-Command -Name "ninja"

$ResolvedM4 = Find-M4Executable -RequestedPath $M4Executable
if (-not $ResolvedM4) {
    throw @"
GNU M4 is required to build pFUnit's gFTL dependency but was not found.

Install MSYS2, then run:
  & C:\msys64\usr\bin\pacman.exe -S --needed m4

Rerun this installer with:
  -M4Executable C:\msys64\usr\bin\m4.exe
"@
}
$M4Version = & $ResolvedM4 --version 2>&1 | Out-String
if ($LASTEXITCODE -ne 0) {
    throw "Unable to run M4: $ResolvedM4"
}
if ($M4Version -notmatch "GNU M4") {
    throw "pFUnit requires GNU M4, but '$ResolvedM4' did not identify itself as GNU M4."
}
Write-Host "Using GNU M4 at $ResolvedM4" -ForegroundColor Cyan

New-Item -ItemType Directory -Force $BuildDir | Out-Null

if (-not (Test-Path -LiteralPath $PFUnitSourceDir -PathType Container)) {
    New-Item -ItemType Directory -Force (Split-Path $PFUnitSourceDir -Parent) | Out-Null
    Invoke-External -FilePath "git" -ArgumentList @(
        "clone", "--branch", "v$PFUnitVersion", "--depth", "1",
        "--recurse-submodules",
        "https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git",
        $PFUnitSourceDir
    )
}
elseif (-not (Test-Path -LiteralPath (Join-Path $PFUnitSourceDir ".git") -PathType Container)) {
    throw "Existing pFUnit source directory is not a Git checkout: $PFUnitSourceDir"
}

Invoke-External -FilePath "git" -WorkingDirectory $PFUnitSourceDir -ArgumentList @(
    "submodule", "update", "--init", "--recursive"
)
Disable-GFTLDependencySelfTests -PFUnitSource $PFUnitSourceDir

$PFUnitConfigureArgs = @(
    "-S", $PFUnitSourceDir,
    "-B", $PFUnitBuildDir,
    "-G", "Ninja",
    "-DCMAKE_BUILD_TYPE=Release",
    "-DCMAKE_INSTALL_PREFIX=$PFUnitInstallDir",
    "-DCMAKE_Fortran_COMPILER=ifx",
    "-DCMAKE_C_COMPILER=cl",
    "-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded",
    "-DPython_EXECUTABLE=$PythonExecutable",
    "-DM4=$ResolvedM4",
    "-DENABLE_TESTS=OFF",
    "-DSKIP_MPI=ON",
    "-DSKIP_OPENMP=ON",
    "-DSKIP_FHAMCREST=ON",
    "-DSKIP_ROBUST=ON"
)
Invoke-External -FilePath "cmake" -ArgumentList $PFUnitConfigureArgs
Invoke-External -FilePath "cmake" -ArgumentList @("--build", $PFUnitBuildDir, "--parallel", "$Parallel")
Invoke-External -FilePath "cmake" -ArgumentList @("--install", $PFUnitBuildDir)

$PFUnitConfig = Get-ChildItem -Path $PFUnitInstallDir -Recurse -Filter "PFUNITConfig.cmake" |
    Select-Object -First 1
if (-not $PFUnitConfig) {
    throw "pFUnit installed, but PFUNITConfig.cmake was not found below $PFUnitInstallDir."
}
$PFUnitDir = $PFUnitConfig.Directory.FullName
Write-Host "Using pFUnit package at $PFUnitDir" -ForegroundColor Cyan

$CEAConfigureArgs = @(
    "-S", $RepoRoot,
    "-B", $CEABuildDir,
    "-G", "Ninja",
    "-DCMAKE_BUILD_TYPE=Release",
    "-DCMAKE_INSTALL_PREFIX=$InstallPrefix",
    "-DCMAKE_Fortran_COMPILER=ifx",
    "-DCMAKE_C_COMPILER=cl",
    "-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded",
    "-DPython3_EXECUTABLE=$PythonExecutable",
    "-DPYTHON_EXECUTABLE=$PythonExecutable",
    "-DPFUNIT_DIR=$PFUnitDir",
    "-DCEA_BUILD_TESTING=ON",
    "-DCEA_ENABLE_BIND_C=ON",
    "-DCEA_ENABLE_BIND_CXX=OFF",
    "-DCEA_ENABLE_BIND_PYTHON=ON",
    "-DCEA_ENABLE_BIND_MATLAB=ON",
    "-DCEA_ENABLE_BIND_EXCEL=OFF"
)
Invoke-External -FilePath "cmake" -ArgumentList $CEAConfigureArgs
Invoke-External -FilePath "cmake" -ArgumentList @("--build", $CEABuildDir, "--parallel", "$Parallel")

$RegisteredTests = & ctest --test-dir $CEABuildDir -N 2>&1 | Out-String
if ($LASTEXITCODE -ne 0) {
    throw "Unable to list CEA tests."
}
if ($RegisteredTests -notmatch "cea_core_test") {
    throw "cea_core_test was not registered; CMake did not load the pFUnit installation."
}

if (-not $SkipTests) {
    Invoke-External -FilePath "ctest" -ArgumentList @(
        "--test-dir", $CEABuildDir, "--output-on-failure"
    )
}

Invoke-External -FilePath "cmake" -ArgumentList @("--install", $CEABuildDir)

$PreviousCMakeArgs = [Environment]::GetEnvironmentVariable("CMAKE_ARGS", "Process")
try {
    $env:CMAKE_ARGS = @(
        "-DCMAKE_Fortran_COMPILER=ifx",
        "-DCMAKE_C_COMPILER=cl",
        "-DCEA_BUILD_TESTING=OFF",
        "-DCEA_ENABLE_BIND_C=ON",
        "-DCEA_ENABLE_BIND_CXX=OFF",
        "-DCEA_ENABLE_BIND_PYTHON=ON",
        "-DCEA_ENABLE_BIND_MATLAB=ON",
        "-DCEA_ENABLE_BIND_EXCEL=OFF"
    ) -join " "

    Invoke-External -FilePath $PythonExecutable -WorkingDirectory $RepoRoot -ArgumentList @(
        "-m", "pip", "install", "--no-build-isolation", "--no-deps",
        "--force-reinstall", "."
    )
}
finally {
    if ($null -eq $PreviousCMakeArgs) {
        Remove-Item Env:CMAKE_ARGS -ErrorAction SilentlyContinue
    }
    else {
        $env:CMAKE_ARGS = $PreviousCMakeArgs
    }
}

Invoke-External -FilePath $PythonExecutable -ArgumentList @(
    "-c", "import cea, cea.matlab; print('CEA', cea.__version__, '- Python and MATLAB wrappers loaded')"
)

if (-not $SkipTests) {
    Invoke-External -FilePath $PythonExecutable -WorkingDirectory $RepoRoot -ArgumentList @(
        "-m", "pytest", "source\bind\python\tests"
    )
}

Write-Host ""
Write-Host "CEA installation completed successfully." -ForegroundColor Green
Write-Host "Native install: $InstallPrefix"
Write-Host "Python:         $PythonExecutable"
Write-Host ""
Write-Host "For a PowerShell session, expose the native CEA executable with:"
Write-Host "  `$env:PATH = `"$InstallPrefix\bin;$InstallPrefix\lib;`$env:PATH`""
Write-Host "  `$env:CEA_EXE = `"$CEABuildDir\source\cea.exe`""
Write-Host ""
Write-Host "For an Intel oneAPI cmd.exe session, use:"
Write-Host "  set `"PATH=$InstallPrefix\bin;$InstallPrefix\lib;%PATH%`""
Write-Host "  set `"CEA_EXE=$CEABuildDir\source\cea.exe`""
