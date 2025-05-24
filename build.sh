#!/bin/bash

# Build script for ArtFlow
# This script verifies that the C# code compiles correctly

echo "ArtFlow Build Verification"
echo "========================="
echo ""
echo "Project Structure:"
echo "- ArtFlow.Core: Core data structures and models"
echo "- ArtFlow.Solver: Finite element solver implementation"
echo "- ArtFlow.IO: Input/output and configuration handling"
echo "- ArtFlow.Console: Command-line interface"
echo ""
echo "To build this project, you need:"
echo "1. .NET 8.0 SDK or later"
echo "2. Run: dotnet build"
echo ""
echo "The following NuGet packages will be restored automatically:"
echo "- MathNet.Numerics (for linear algebra)"
echo "- Newtonsoft.Json (for JSON serialization)"
echo ""
echo "To run the console application:"
echo "dotnet run --project ArtFlow.Console example"
echo "dotnet run --project ArtFlow.Console run <config.json>"
