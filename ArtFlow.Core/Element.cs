namespace ArtFlow.Core;

/// <summary>
/// Represents a finite element connecting two nodes
/// </summary>
public class Element
{
    public int Id { get; set; }
    public Node Node1 { get; set; }
    public Node Node2 { get; set; }
    
    public double Length => Math.Abs(Node2.Position - Node1.Position);
    
    public Element(int id, Node node1, Node node2)
    {
        Id = id;
        Node1 = node1;
        Node2 = node2;
    }
    
    /// <summary>
    /// Get shape function value for node at local coordinate xi (-1 to 1)
    /// </summary>
    public double ShapeFunction(int nodeIndex, double xi)
    {
        if (nodeIndex == 0) // Node1
            return 0.5 * (1.0 - xi);
        else // Node2
            return 0.5 * (1.0 + xi);
    }
    
    /// <summary>
    /// Get shape function derivative with respect to xi
    /// </summary>
    public double ShapeFunctionDerivative(int nodeIndex, double xi)
    {
        if (nodeIndex == 0) // Node1
            return -0.5;
        else // Node2
            return 0.5;
    }
}
