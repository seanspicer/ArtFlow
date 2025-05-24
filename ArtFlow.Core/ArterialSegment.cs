namespace ArtFlow.Core;

/// <summary>
/// Represents an arterial segment in the vascular network
/// </summary>
public class ArterialSegment
{
    public int Id { get; set; }
    public string Name { get; set; }
    public double Length { get; set; }
    public double InitialRadius { get; set; }
    public double InitialArea => Math.PI * InitialRadius * InitialRadius;
    public Material Material { get; set; }
    public List<Node> Nodes { get; set; }
    public List<Element> Elements { get; set; }
    
    // Physical properties
    public double BloodDensity { get; set; } = 1.055; // g/cm^3
    public double KinematicViscosity { get; set; } = 0.04; // cm^2/s
    
    // Profile parameters for velocity distribution
    public double ProfileParameter { get; set; } = 2.0; // n=2 for Poiseuille
    
    // Initial conditions
    public double InitialPressure { get; set; }
    
    public ArterialSegment(int id, string name, double length, double initialRadius)
    {
        Id = id;
        Name = name ?? $"Segment_{id}";
        Length = length;
        InitialRadius = initialRadius;
        Material = new Material();
        Nodes = new List<Node>();
        Elements = new List<Element>();
    }
    
    /// <summary>
    /// Calculate delta parameter for velocity profile
    /// </summary>
    public double GetDelta()
    {
        return 1.0 / (1.0 + ProfileParameter);
    }
    
    /// <summary>
    /// Calculate N parameter for viscous losses
    /// </summary>
    public double GetN()
    {
        return -2.0 * (ProfileParameter + 2.0) * Math.PI * KinematicViscosity;
    }
    
    /// <summary>
    /// Calculate pressure from area using Olufsen's constitutive equation
    /// </summary>
    public double CalculatePressure(double area, double z)
    {
        double r0 = InitialRadius;
        double S0 = InitialArea;
        double Eh_r0 = Material.GetElasticModulus(r0);
        
        return InitialPressure + (4.0 / 3.0) * Eh_r0 * (1.0 - Math.Sqrt(S0 / area));
    }
    
    /// <summary>
    /// Calculate derivative of pressure with respect to area
    /// </summary>
    public double CalculateDPressureDArea(double area, double z)
    {
        double r0 = InitialRadius;
        double S0 = InitialArea;
        double Eh_r0 = Material.GetElasticModulus(r0);
        
        return (4.0 / 3.0) * Eh_r0 * (0.5 * Math.Sqrt(S0) / Math.Pow(area, 1.5));
    }
}
