using ArtFlow.Core;
using ArtFlow.Solver;
using Newtonsoft.Json;

namespace ArtFlow.IO;

/// <summary>
/// Configuration for a simulation run
/// </summary>
public class SimulationConfig
{
    public string Name { get; set; } = "Simulation";
    public SolverSettings SolverSettings { get; set; } = new();
    public List<SegmentConfig> Segments { get; set; } = new();
    public List<BranchConfig> Branches { get; set; } = new();
    public List<BoundaryConditionConfig> BoundaryConditions { get; set; } = new();
}

public class SegmentConfig
{
    public int Id { get; set; }
    public string Name { get; set; } = "";
    public double Length { get; set; }
    public double InitialRadius { get; set; }
    public double InitialPressure { get; set; } = 100.0; // mmHg
    public MaterialConfig? Material { get; set; }
}

public class MaterialConfig
{
    public double K1 { get; set; } = 2e7;
    public double K2 { get; set; } = -22.53;
    public double K3 { get; set; } = 8.65e5;
}

public class BranchConfig
{
    public int Id { get; set; }
    public List<int> InletNodeIds { get; set; } = new();
    public List<int> OutletNodeIds { get; set; } = new();
}

public class BoundaryConditionConfig
{
    public int NodeId { get; set; }
    public string Type { get; set; } = "FlowRate"; // FlowRate, Pressure, Resistance
    public object? Parameters { get; set; }
}
