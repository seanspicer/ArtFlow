namespace ArtFlow.Core;

/// <summary>
/// Represents a branch point where multiple segments connect
/// Enforces pressure continuity and mass conservation
/// </summary>
public class BranchPoint
{
    public int Id { get; set; }
    public List<int> InletNodeIds { get; set; }
    public List<int> OutletNodeIds { get; set; }
    
    // Lagrange multipliers for constraints
    public double FlowRateLagrangeMultiplier { get; set; }
    public List<double> PressureLagrangeMultipliers { get; set; }
    
    public BranchPoint(int id)
    {
        Id = id;
        InletNodeIds = new List<int>();
        OutletNodeIds = new List<int>();
        PressureLagrangeMultipliers = new List<double>();
    }
    
    public int TotalConstraints => InletNodeIds.Count + OutletNodeIds.Count - 1;
}
