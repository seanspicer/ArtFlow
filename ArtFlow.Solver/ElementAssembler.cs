using ArtFlow.Core;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace ArtFlow.Solver;

/// <summary>
/// Assembles element matrices and residuals for the finite element method
/// </summary>
public class ElementAssembler
{
    private readonly ArterialSegment _segment;
    private readonly SolverSettings _settings;
    
    public ElementAssembler(ArterialSegment segment, SolverSettings settings)
    {
        _segment = segment;
        _settings = settings;
    }
    
    /// <summary>
    /// Compute element matrix and residual vector
    /// </summary>
    public (Matrix<double>, Vector<double>) ComputeElementMatrixAndResidual(Element element, double dt)
    {
        var Ke = Matrix<double>.Build.Dense(4, 4);
        var Re = Vector<double>.Build.Dense(4);
        
        // Get nodal values
        var U = GetNodalValues(element);
        var Un = GetNodalValuesPrev(element);
        
        // Gauss quadrature points (2-point rule)
        double[] xi = { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
        double[] w = { 1.0, 1.0 };
        
        for (int g = 0; g < 2; g++)
        {
            // Evaluate shape functions at Gauss point
            double N1 = element.ShapeFunction(0, xi[g]);
            double N2 = element.ShapeFunction(1, xi[g]);
            double dN1dxi = element.ShapeFunctionDerivative(0, xi[g]);
            double dN2dxi = element.ShapeFunctionDerivative(1, xi[g]);
            
            // Transform derivatives to physical coordinates
            double J = element.Length / 2.0; // Jacobian
            double dN1dz = dN1dxi / J;
            double dN2dz = dN2dxi / J;

            // Interpolate solution at Gauss point
            double S = N1 * U[0] + N2 * U[2];
            double Q = N1 * U[1] + N2 * U[3];
            double Sn = N1 * Un[0] + N2 * Un[2];
            double Qn = N1 * Un[1] + N2 * Un[3];
            
            // Calculate matrices A, K, G
            var (A, K, G) = CalculateMatrices(S, Q);
            
            // Shape function arrays
            var N = Vector<double>.Build.DenseOfArray(new[] { N1, 0, N2, 0, 0, N1, 0, N2 });
            var B = Vector<double>.Build.DenseOfArray(new[] { dN1dz, 0, dN2dz, 0, 0, dN1dz, 0, dN2dz });
            
            // Time term: (U - Un) / dt
            var timeTerm = (U - Un) / dt;
            
            // Advection term: A * dU/dz
            var dUdz = Vector<double>.Build.DenseOfArray(new[] {
                dN1dz * U[0] + dN2dz * U[2],
                dN1dz * U[1] + dN2dz * U[3]
            });
            var advectionTerm = A * dUdz;
            
            // Diffusion term: K * d²U/dz²
            // For linear elements, second derivative is zero
            
            // Source term G
            
            // Assemble element contributions
            double wJ = w[g] * J;
            
            // Add contributions to element matrix and residual
            AddElementContributions(Ke, Re, N, B, A, K, timeTerm, advectionTerm, G, wJ, dt);
        }
        
        // Add stabilization if enabled
        if (_settings.UseStabilization)
        {
            AddStabilization(element, Ke, Re, U, dt);
        }
        
        return (Ke, Re);
    }

    private (Matrix<double>, Matrix<double>, Vector<double>) CalculateMatrices(double S, double Q)
    {
        // Matrix A (quasi-linear form)
        var A = Matrix<double>.Build.Dense(2, 2);
        double delta = _segment.GetDelta();
        double rho = _segment.BloodDensity;
        double dpdS = _segment.CalculateDPressureDArea(S, 0); // TODO: pass actual z
        
        A[0, 0] = 0;
        A[0, 1] = 1;
        A[1, 0] = -(1 + delta) * Q * Q / (S * S) + S / rho * dpdS;
        A[1, 1] = (1 + delta) * 2 * Q / S;
        
        // Matrix K (diffusion)
        var K = Matrix<double>.Build.Dense(2, 2);
        K[0, 0] = 0;
        K[0, 1] = 0;
        K[1, 0] = 0;
        K[1, 1] = _segment.KinematicViscosity;
        
        // Vector G (source terms)
        var G = Vector<double>.Build.Dense(2);
        double N = _segment.GetN();
        G[0] = 0; // No outflow (c = 0)
        G[1] = 0 + N * Q / S; // External force f = 0, viscous term
        
        return (A, K, G);
    }
    
    private Vector<double> GetNodalValues(Element element)
    {
        return Vector<double>.Build.DenseOfArray(new[] {
            element.Node1.Area,
            element.Node1.FlowRate,
            element.Node2.Area,
            element.Node2.FlowRate
        });
    }
    
    private Vector<double> GetNodalValuesPrev(Element element)
    {
        return Vector<double>.Build.DenseOfArray(new[] {
            element.Node1.AreaPrev,
            element.Node1.FlowRatePrev,
            element.Node2.AreaPrev,
            element.Node2.FlowRatePrev
        });
    }

    private void AddElementContributions(Matrix<double> Ke, Vector<double> Re, 
        Vector<double> N, Vector<double> B, Matrix<double> A, Matrix<double> K, 
        Vector<double> timeTerm, Vector<double> advectionTerm, Vector<double> G, 
        double wJ, double dt)
    {
        // Implementation of element matrix assembly
        // This is simplified - full implementation would include all terms
        
        // Time term contribution
        for (int i = 0; i < 4; i += 2)
        {
            for (int j = 0; j < 4; j += 2)
            {
                Ke[i, j] += N[i] * N[j] / dt * wJ;
                Ke[i+1, j+1] += N[i+1] * N[j+1] / dt * wJ;
            }
        }
        
        // Advection term contribution
        // TODO: Complete implementation
        
        // Add to residual
        for (int i = 0; i < 4; i++)
        {
            Re[i] += N[i] * timeTerm[i/2] * wJ;
        }
    }
    
    private void AddStabilization(Element element, Matrix<double> Ke, Vector<double> Re, Vector<double> U, double dt)
    {
        // SUPG stabilization
        // tau = stabilization parameter
        double h = element.Length;
        double tau = _settings.StabilizationParameter * h * h / dt;
        
        // TODO: Implement stabilization terms following the paper
    }
}
