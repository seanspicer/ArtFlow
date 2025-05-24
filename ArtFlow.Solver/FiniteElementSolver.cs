using ArtFlow.Core;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace ArtFlow.Solver;

/// <summary>
/// Space-time finite element solver for 1D blood flow
/// Based on discontinuous Galerkin in time, continuous piecewise linear in space
/// </summary>
public class FiniteElementSolver
{
    private readonly SolverSettings _settings;
    private readonly List<ArterialSegment> _segments;
    private readonly List<BranchPoint> _branchPoints;
    private readonly List<BoundaryCondition> _boundaryConditions;
    
    // Global system
    private int _totalDofs;
    private Matrix<double> _globalMatrix;
    private Vector<double> _globalResidual;
    private Vector<double> _solution;
    
    public FiniteElementSolver(SolverSettings settings)
    {
        _settings = settings;
        _segments = new List<ArterialSegment>();
        _branchPoints = new List<BranchPoint>();
        _boundaryConditions = new List<BoundaryCondition>();
    }
    
    public void AddSegment(ArterialSegment segment)
    {
        _segments.Add(segment);
    }
    
    public void AddBranchPoint(BranchPoint branch)
    {
        _branchPoints.Add(branch);
    }
    
    public void AddBoundaryCondition(BoundaryCondition bc)
    {
        _boundaryConditions.Add(bc);
    }

    /// <summary>
    /// Initialize the mesh for all segments
    /// </summary>
    public void InitializeMesh()
    {
        foreach (var segment in _segments)
        {
            int numElements = (int)(segment.Length / _settings.ElementSize);
            if (numElements < 1) numElements = 1;
            
            double dz = segment.Length / numElements;
            
            // Create nodes
            for (int i = 0; i <= numElements; i++)
            {
                var node = new Node(segment.Nodes.Count, i * dz)
                {
                    Area = segment.InitialArea,
                    FlowRate = 0.0,
                    AreaPrev = segment.InitialArea,
                    FlowRatePrev = 0.0
                };
                segment.Nodes.Add(node);
            }
            
            // Create elements
            for (int i = 0; i < numElements; i++)
            {
                var element = new Element(i, segment.Nodes[i], segment.Nodes[i + 1]);
                segment.Elements.Add(element);
            }
        }
        
        CalculateTotalDofs();
    }
    
    private void CalculateTotalDofs()
    {
        _totalDofs = 0;
        foreach (var segment in _segments)
        {
            _totalDofs += segment.Nodes.Count * 2; // S and Q for each node
        }
        
        // Add Lagrange multipliers for branch points
        foreach (var branch in _branchPoints)
        {
            _totalDofs += branch.TotalConstraints;
        }
    }

    /// <summary>
    /// Solve the system for one time step
    /// </summary>
    public void SolveTimeStep(double currentTime, double dt)
    {
        _globalMatrix = Matrix<double>.Build.Sparse(_totalDofs, _totalDofs);
        _globalResidual = Vector<double>.Build.Dense(_totalDofs);
        _solution = Vector<double>.Build.Dense(_totalDofs);
        
        // Initialize solution vector with current values
        InitializeSolutionVector();
        
        // Newton-Raphson iteration
        for (int iter = 0; iter < _settings.MaxNewtonIterations; iter++)
        {
            AssembleGlobalSystem(currentTime, dt);
            ApplyBoundaryConditions(currentTime + dt);
            
            // Solve linear system
            var deltaU = _globalMatrix.Solve(-_globalResidual);
            _solution += deltaU;
            
            // Update nodal values
            UpdateNodalValues();
            
            // Check convergence
            if (deltaU.L2Norm() < _settings.NewtonTolerance)
            {
                break;
            }
        }
        
        // Update previous time step values
        UpdatePreviousValues();
    }

    private void InitializeSolutionVector()
    {
        int dofIndex = 0;
        foreach (var segment in _segments)
        {
            foreach (var node in segment.Nodes)
            {
                _solution[dofIndex++] = node.Area;
                _solution[dofIndex++] = node.FlowRate;
            }
        }
    }
    
    private void UpdateNodalValues()
    {
        int dofIndex = 0;
        foreach (var segment in _segments)
        {
            foreach (var node in segment.Nodes)
            {
                node.Area = _solution[dofIndex++];
                node.FlowRate = _solution[dofIndex++];
            }
        }
    }
    
    private void UpdatePreviousValues()
    {
        foreach (var segment in _segments)
        {
            foreach (var node in segment.Nodes)
            {
                node.AreaPrev = node.Area;
                node.FlowRatePrev = node.FlowRate;
            }
        }
    }
            }
        }
        
        // Add branch point constraints
        foreach (var branch in _branchPoints)
        {
            ApplyBranchConstraints(branch);
        }
    }
    
    private void AddElementContribution(Element element, Matrix<double> Ke, Vector<double> Re)
    {
        // Map local DOFs to global DOFs
        var globalDofs = GetGlobalDofs(element);
        
        for (int i = 0; i < 4; i++)
        {
            _globalResidual[globalDofs[i]] += Re[i];
            
            for (int j = 0; j < 4; j++)
            {
                _globalMatrix[globalDofs[i], globalDofs[j]] += Ke[i, j];
            }
        }
    }
    
    private int[] GetGlobalDofs(Element element)
    {
        // Each node has 2 DOFs (S and Q)
        // Find global node indices
        int node1GlobalIndex = GetGlobalNodeIndex(element.Node1);
        int node2GlobalIndex = GetGlobalNodeIndex(element.Node2);
        
        return new int[]
        {
            node1GlobalIndex * 2,     // S1
            node1GlobalIndex * 2 + 1, // Q1
            node2GlobalIndex * 2,     // S2
            node2GlobalIndex * 2 + 1  // Q2
        };
    }
    
    private int GetGlobalNodeIndex(Node node)
    {
        int index = 0;
        foreach (var segment in _segments)
        {
            var nodeIndex = segment.Nodes.IndexOf(node);
            if (nodeIndex >= 0)
                return index + nodeIndex;
            index += segment.Nodes.Count;
        }
        return -1;
    }

    private void ApplyBranchConstraints(BranchPoint branch)
    {
        // Implementation of Lagrange multiplier constraints for branch points
        // This enforces pressure continuity and mass conservation
        
        // TODO: Implement branch constraint assembly
        // This involves modifying the global matrix and residual
        // to include Lagrange multiplier terms
    }
    
    private void ApplyBoundaryConditions(double time)
    {
        foreach (var bc in _boundaryConditions)
        {
            int nodeGlobalIndex = GetGlobalNodeIndexById(bc.NodeId);
            if (nodeGlobalIndex < 0) continue;
            
            switch (bc.Type)
            {
                case BoundaryConditionType.FlowRate:
                    ApplyFlowRateBC(nodeGlobalIndex, bc as FlowRateBoundaryCondition, time);
                    break;
                case BoundaryConditionType.Pressure:
                    ApplyPressureBC(nodeGlobalIndex, bc as PressureBoundaryCondition, time);
                    break;
                case BoundaryConditionType.Resistance:
                    ApplyResistanceBC(nodeGlobalIndex, bc as ResistanceBoundaryCondition, time);
                    break;
            }
        }
    }
    
    private int GetGlobalNodeIndexById(int nodeId)
    {
        int index = 0;
        foreach (var segment in _segments)
        {
            foreach (var node in segment.Nodes)
            {
                if (node.Id == nodeId)
                    return index;
                index++;
            }
        }
        return -1;
    }

    private void ApplyFlowRateBC(int nodeIndex, FlowRateBoundaryCondition bc, double time)
    {
        int qDof = nodeIndex * 2 + 1; // Q DOF
        
        // Clear row and column
        for (int i = 0; i < _totalDofs; i++)
        {
            _globalMatrix[qDof, i] = 0;
            _globalMatrix[i, qDof] = 0;
        }
        
        // Set diagonal to 1
        _globalMatrix[qDof, qDof] = 1;
        
        // Set residual
        double prescribedQ = bc.FlowRateFunction(time);
        _globalResidual[qDof] = _solution[qDof] - prescribedQ;
    }
    
    private void ApplyPressureBC(int nodeIndex, PressureBoundaryCondition bc, double time)
    {
        // Convert pressure to area BC
        int sDof = nodeIndex * 2; // S DOF
        
        // Find the segment containing this node
        ArterialSegment? segment = null;
        Node? node = null;
        foreach (var seg in _segments)
        {
            node = seg.Nodes.FirstOrDefault(n => GetGlobalNodeIndex(n) == nodeIndex);
            if (node != null)
            {
                segment = seg;
                break;
            }
        }
        
        if (segment == null || node == null) return;
        
        // Calculate area from pressure using inverse constitutive relation
        double prescribedP = bc.PressureFunction(time);
        double targetArea = CalculateAreaFromPressure(segment, prescribedP, node.Position);
        
        // Apply as area BC
        for (int i = 0; i < _totalDofs; i++)
        {
            _globalMatrix[sDof, i] = 0;
            _globalMatrix[i, sDof] = 0;
        }
        _globalMatrix[sDof, sDof] = 1;
        _globalResidual[sDof] = _solution[sDof] - targetArea;
    }

    private void ApplyResistanceBC(int nodeIndex, ResistanceBoundaryCondition bc, double time)
    {
        // Resistance BC: p - QR = 0
        // This modifies both S and Q equations
        
        // TODO: Implement resistance BC following the paper's approach
        // This involves modifying the Jacobian matrix appropriately
    }
    
    private double CalculateAreaFromPressure(ArterialSegment segment, double pressure, double z)
    {
        // Inverse of Olufsen's constitutive equation
        // p = p0 + (4/3) * (Eh/r0) * (1 - sqrt(S0/S))
        // Solve for S
        
        double r0 = segment.InitialRadius;
        double S0 = segment.InitialArea;
        double Eh_r0 = segment.Material.GetElasticModulus(r0);
        double p0 = segment.InitialPressure;
        
        double term = 1.0 - (pressure - p0) / ((4.0 / 3.0) * Eh_r0);
        return S0 / (term * term);
    }
    
    /// <summary>
    /// Run the complete simulation
    /// </summary>
    public void RunSimulation()
    {
        double currentTime = 0;
        int timeStep = 0;
        
        InitializeMesh();
        
        while (currentTime < _settings.TotalTime)
        {
            double dt = _settings.TimeStepSize;
            
            SolveTimeStep(currentTime, dt);
            
            currentTime += dt;
            timeStep++;
            
            // Output results periodically
            if (timeStep % _settings.OutputInterval == 0)
            {
                OutputResults(currentTime, timeStep);
            }
        }
    }
    
    private void OutputResults(double time, int timeStep)
    {
        // TODO: Implement result output
        Console.WriteLine($"Time step {timeStep}, t = {time:F4}");
    }
}
