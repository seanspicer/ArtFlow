namespace ArtFlow.Core;

/// <summary>
/// Represents a node in the finite element mesh
/// </summary>
public class Node
{
    public int Id { get; set; }
    public double Position { get; set; } // z-coordinate along vessel
    public double Area { get; set; } // Cross-sectional area S
    public double FlowRate { get; set; } // Volumetric flow rate Q
    
    // Previous time step values
    public double AreaPrev { get; set; }
    public double FlowRatePrev { get; set; }
    
    public Node(int id, double position)
    {
        Id = id;
        Position = position;
    }
}
