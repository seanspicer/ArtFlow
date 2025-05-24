namespace ArtFlow.Core;

/// <summary>
/// Base class for boundary conditions
/// </summary>
public abstract class BoundaryCondition
{
    public int NodeId { get; set; }
    public abstract BoundaryConditionType Type { get; }
    
    protected BoundaryCondition(int nodeId)
    {
        NodeId = nodeId;
    }
}

public enum BoundaryConditionType
{
    FlowRate,
    Pressure,
    Resistance
}

/// <summary>
/// Prescribed flow rate boundary condition
/// </summary>
public class FlowRateBoundaryCondition : BoundaryCondition
{
    public override BoundaryConditionType Type => BoundaryConditionType.FlowRate;
    public Func<double, double> FlowRateFunction { get; set; }
    
    public FlowRateBoundaryCondition(int nodeId, Func<double, double> flowRateFunc) : base(nodeId)
    {
        FlowRateFunction = flowRateFunc;
    }
}

/// <summary>
/// Prescribed pressure boundary condition
/// </summary>
public class PressureBoundaryCondition : BoundaryCondition
{
    public override BoundaryConditionType Type => BoundaryConditionType.Pressure;
    public Func<double, double> PressureFunction { get; set; }
    
    public PressureBoundaryCondition(int nodeId, Func<double, double> pressureFunc) : base(nodeId)
    {
        PressureFunction = pressureFunc;
    }
}

/// <summary>
/// Resistance boundary condition (p - QR = 0)
/// </summary>
public class ResistanceBoundaryCondition : BoundaryCondition
{
    public override BoundaryConditionType Type => BoundaryConditionType.Resistance;
    public double Resistance { get; set; }
    
    public ResistanceBoundaryCondition(int nodeId, double resistance) : base(nodeId)
    {
        Resistance = resistance;
    }
}
