#!/usr/bin/env python3
"""
Verification script for ArtFlow C# implementation
Checks that all required files are present and have proper structure
"""

import os
import json

def check_file_exists(filepath, description):
    if os.path.exists(filepath):
        print(f"✓ {description}: {filepath}")
        return True
    else:
        print(f"✗ {description}: {filepath} NOT FOUND")
        return False

def check_project_structure():
    print("Checking ArtFlow Project Structure")
    print("==================================\n")
    
    base_dir = "/Users/sspicer/work/ArtFlow"
    
    # Check solution file
    check_file_exists(f"{base_dir}/ArtFlow.sln", "Solution file")
    
    # Check project files
    projects = ["ArtFlow.Core", "ArtFlow.Solver", "ArtFlow.IO", "ArtFlow.Console"]
    for project in projects:
        print(f"\n{project}:")
        check_file_exists(f"{base_dir}/{project}/{project}.csproj", "Project file")
    
    # Check core source files
    print("\nCore Components:")
    core_files = [
        "Node.cs", "Element.cs", "ArterialSegment.cs", "Material.cs",
        "BoundaryConditions.cs", "BranchPoint.cs", "Stenosis.cs"
    ]
    for file in core_files:
        check_file_exists(f"{base_dir}/ArtFlow.Core/{file}", file)
    
    # Check solver files
    print("\nSolver Components:")
    solver_files = [
        "FiniteElementSolver.cs", "ElementAssembler.cs", "SolverSettings.cs"
    ]
    for file in solver_files:
        check_file_exists(f"{base_dir}/ArtFlow.Solver/{file}", file)
    
    # Check IO files
    print("\nIO Components:")
    io_files = ["SimulationConfig.cs", "ConfigLoader.cs"]
    for file in io_files:
        check_file_exists(f"{base_dir}/ArtFlow.IO/{file}", file)
    
    # Check console files
    print("\nConsole Components:")
    check_file_exists(f"{base_dir}/ArtFlow.Console/Program.cs", "Program.cs")
    
    # Check examples
    print("\nExample Configurations:")
    example_files = ["carotid_test.json", "bifurcation_example.json"]
    for file in example_files:
        if check_file_exists(f"{base_dir}/examples/{file}", file):
            # Validate JSON
            try:
                with open(f"{base_dir}/examples/{file}", 'r') as f:
                    json.load(f)
                print(f"  ✓ {file} is valid JSON")
            except json.JSONDecodeError as e:
                print(f"  ✗ {file} has invalid JSON: {e}")
    
    print("\n" + "="*40)
    print("Verification complete!")
    print("\nTo build the project, run:")
    print("  cd /Users/sspicer/work/ArtFlow")
    print("  dotnet build")
    print("\nTo run an example:")
    print("  dotnet run --project ArtFlow.Console run examples/carotid_test.json")

if __name__ == "__main__":
    check_project_structure()
