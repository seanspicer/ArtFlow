namespace ArtFlow.Solver;

/// <summary>
/// Settings for the finite element solver
/// </summary>
public class SolverSettings
{
    public double ElementSize { get; set; } = 1.0; // cm
    public double TimeStepSize { get; set; } = 0.001; // seconds
    public double TotalTime { get; set; } = 1.0; // seconds
    public int MaxNewtonIterations { get; set; } = 10;
    public double NewtonTolerance { get; set; } = 1e-6;
    public int OutputInterval { get; set; } = 100; // Output every N time steps
    
    // Stabilization parameters
    public bool UseStabilization { get; set; } = true;
    public double StabilizationParameter { get; set; } = 0.5;
}
