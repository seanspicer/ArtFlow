namespace ArtFlow.Core;

/// <summary>
/// Represents a stenosis (narrowing) in an arterial segment
/// </summary>
public class Stenosis
{
    public double StartPosition { get; set; }
    public double EndPosition { get; set; }
    public double Length => EndPosition - StartPosition;
    public double MinimumRadius { get; set; }
    public double ProximalRadius { get; set; }
    
    public Stenosis(double startPos, double endPos, double minRadius, double proxRadius)
    {
        StartPosition = startPos;
        EndPosition = endPos;
        MinimumRadius = minRadius;
        ProximalRadius = proxRadius;
    }
    
    /// <summary>
    /// Calculate modified viscous loss parameter N for stenosis
    /// Based on Seeley and Young empirical formula
    /// </summary>
    public double CalculateStenosisN(double Q0, double S0, double S1, double L, double D0)
    {
        double V0 = Q0 / S0;
        double Re0 = V0 * D0 / 0.04; // kinematic viscosity = 0.04
        
        double Kt = 1.52;
        double Kv = 32.0 * L / D0 * Math.Pow(S0 / S1, 2);
        
        double term1 = -S1 * S1 * Q0 * Q0 / (S0 * S0 * Q0 * L);
        double term2 = Kv / Re0 + 0.5 * Kt * Math.Pow(S0 / S1 - 1.0, 2);
        
        return term1 * term2;
    }
}
