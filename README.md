# ArtFlow - 1D Finite Element Blood Flow Simulation

ArtFlow is a C# implementation of the one-dimensional finite element method for simulating blood flow in arterial networks, based on the paper "A One-dimensional Finite Element Method for Simulation-based Medical Planning for Cardiovascular Disease" by Wan et al. (2002).

## Overview

This implementation provides a space-time finite element solver for the 1D blood flow equations, featuring:

- **Discontinuous Galerkin method in time** with piecewise constant approximation
- **Continuous piecewise linear elements in space**
- **Modified Newton-Raphson solver** for the nonlinear system
- **Olufsen's constitutive model** for arterial wall elasticity
- **Support for various boundary conditions**: prescribed flow rate, pressure, and resistance
- **Branch point handling** with Lagrange multipliers for mass conservation and pressure continuity
- **Stenosis modeling** based on empirical pressure loss formulas

## Architecture

The solution is organized into four main projects:

### ArtFlow.Core
Contains the fundamental data structures:
- `Node`: Mesh nodes with area (S) and flow rate (Q) degrees of freedom
- `Element`: Linear finite elements connecting nodes
- `ArterialSegment`: Represents arterial vessels with material properties
- `Material`: Arterial wall properties using Olufsen's parameters
- `BoundaryCondition`: Base class for different BC types
- `BranchPoint`: Junction points enforcing conservation laws
- `Stenosis`: Models for arterial narrowing

### ArtFlow.Solver
Implements the finite element method:
- `FiniteElementSolver`: Main solver class managing the global system
- `ElementAssembler`: Computes element matrices and residuals
- `SolverSettings`: Configuration for numerical parameters

The solver implements:
1. Quasi-linear formulation of the 1D blood flow equations
2. Stabilized finite element method with optional SUPG
3. Time integration using backward Euler
4. Newton-Raphson iteration for nonlinear convergence

### ArtFlow.IO
Handles input/output operations:
- `SimulationConfig`: JSON-serializable configuration structure
- `ConfigLoader`: Loads configurations and builds solver instances

### ArtFlow.Console
Command-line interface for running simulations

## Mathematical Formulation

The solver implements the 1D blood flow equations in quasi-linear form:

```
∂U/∂t + A(U) ∂U/∂z - K ∂²U/∂z² = G(U)
```

Where:
- U = [S, Q]ᵀ (cross-sectional area and flow rate)
- A = advection matrix
- K = diffusion matrix  
- G = source terms

The constitutive relationship (Olufsen model):
```
p = p₀ + (4/3)(Eh/r₀)(1 - √(S₀/S))
```

## Usage

### Creating a Configuration File

Configuration files are JSON documents specifying:
- Solver settings (time step, element size, etc.)
- Arterial segments (geometry and properties)
- Boundary conditions
- Branch points (if any)

Example structure:
```json
{
  "Name": "Carotid Artery Simulation",
  "SolverSettings": {
    "ElementSize": 1.0,
    "TimeStepSize": 0.001,
    "TotalTime": 3.0,
    "MaxNewtonIterations": 10,
    "NewtonTolerance": 1e-6
  },
  "Segments": [
    {
      "Id": 0,
      "Name": "Common Carotid",
      "Length": 16.0,
      "InitialRadius": 0.4,
      "InitialPressure": 100.0,
      "Material": {
        "K1": 2e7,
        "K2": -22.53,
        "K3": 8.65e5
      }
    }
  ],
  "BoundaryConditions": [
    {
      "NodeId": 0,
      "Type": "FlowRate",
      "Parameters": {
        "type": "cardiac",
        "value": 5.0
      }
    },
    {
      "NodeId": 1,
      "Type": "Resistance", 
      "Parameters": {
        "resistance": 1.0
      }
    }
  ]
}
```

### Running a Simulation

```bash
# Generate example configuration
dotnet run --project ArtFlow.Console example

# Run simulation
dotnet run --project ArtFlow.Console run carotid_example.json
```

### Boundary Conditions

#### Flow Rate BC
Prescribes volumetric flow rate Q(t):
```json
{
  "NodeId": 0,
  "Type": "FlowRate",
  "Parameters": {
    "type": "constant|sine|cardiac",
    "value": 5.0
  }
}
```

#### Pressure BC
Prescribes pressure p(t):
```json
{
  "NodeId": 1,
  "Type": "Pressure",
  "Parameters": {
    "value": 100.0
  }
}
```

#### Resistance BC
Enforces p - QR = 0:
```json
{
  "NodeId": 1,
  "Type": "Resistance",
  "Parameters": {
    "resistance": 1.0
  }
}
```

### Material Properties

The default material parameters (K1, K2, K3) are based on experimental data:
- K1 = 2×10⁷ g·s⁻²·cm⁻¹
- K2 = -22.53 cm⁻¹
- K3 = 8.65×10⁵ g·s⁻²·cm⁻¹

These determine the elastic modulus: Eh/r₀ = K1·exp(K2·r₀) + K3

## Implementation Notes

1. **Mesh Generation**: Elements are automatically generated based on the specified element size
2. **Time Stepping**: Uses backward Euler with configurable time step
3. **Stabilization**: Optional SUPG stabilization for advection-dominated flows
4. **Output**: Results are currently written to console; extend the solver for file output as needed

## Limitations and Future Work

- Branch point constraints are partially implemented
- Stenosis model integration needs completion
- Full stabilization terms following the paper need implementation
- Result visualization and output to standard formats (VTK, CSV) should be added
- Impedance boundary conditions could be implemented for more realistic outflow

## Building the Project

Requirements:
- .NET 8.0 SDK or later
- MathNet.Numerics package (automatically restored)
- Newtonsoft.Json package (automatically restored)

```bash
# Build solution
dotnet build

# Run tests (when implemented)
dotnet test

# Run console application
dotnet run --project ArtFlow.Console
```

## References

Wan, J., Steele, B., Spicer, S. A., Strohband, S., Feijóo, G. R., Hughes, T. J., & Taylor, C. A. (2002). A one-dimensional finite element method for simulation-based medical planning for cardiovascular disease. Computer methods in biomechanics and biomedical engineering, 5(3), 195-206.
