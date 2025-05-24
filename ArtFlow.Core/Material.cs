namespace ArtFlow.Core;

/// <summary>
/// Material properties for arterial wall
/// Based on Olufsen's constitutive relationship
/// </summary>
public class Material
{
    public double K1 { get; set; } = 2e7; // g s^-2 cm^-1
    public double K2 { get; set; } = -22.53; // cm^-1
    public double K3 { get; set; } = 8.65e5; // g s^-2 cm^-1
    
    /// <summary>
    /// Calculate Eh/r0 based on initial radius
    /// </summary>
    public double GetElasticModulus(double r0)
    {
        return K1 * Math.Exp(K2 * r0) + K3;
    }
}
