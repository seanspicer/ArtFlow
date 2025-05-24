using ArtFlow.Core;
using ArtFlow.Solver;
using Newtonsoft.Json;

namespace ArtFlow.IO;

/// <summary>
/// Loads simulation configuration from JSON files
/// </summary>
public class ConfigLoader
{
    public static SimulationConfig LoadFromFile(string filePath)
    {
        var json = File.ReadAllText(filePath);
        return JsonConvert.DeserializeObject<SimulationConfig>(json) 
            ?? throw new InvalidOperationException("Failed to load configuration");
    }
    
    public static void SaveToFile(SimulationConfig config, string filePath)
    {
        var json = JsonConvert.SerializeObject(config, Formatting.Indented);
        File.WriteAllText(filePath, json);
    }
    
    /// <summary>
    /// Build solver from configuration
    /// </summary>
    public static FiniteElementSolver BuildSolver(SimulationConfig config)
    {
        var solver = new FiniteElementSolver(config.SolverSettings);
        
        // Add segments
        foreach (var segConfig in config.Segments)
        {
            var segment = new ArterialSegment(
                segConfig.Id,
                segConfig.Name,
                segConfig.Length,
                segConfig.InitialRadius
            )
            {
                InitialPressure = segConfig.InitialPressure
            };
            
            if (segConfig.Material != null)
            {
                segment.Material = new Material
                {
                    K1 = segConfig.Material.K1,
                    K2 = segConfig.Material.K2,
                    K3 = segConfig.Material.K3
                };
            }
                    K1 = segConfig.Material.K1,
                    K2 = segConfig.Material.K2,
                    K3 = segConfig.Material.K3
                };
            }
            
            solver.AddSegment(segment);
        }
        
        // Add branch points
        foreach (var branchConfig in config.Branches)
        {
            var branch = new BranchPoint(branchConfig.Id)
            {
                InletNodeIds = branchConfig.InletNodeIds,
                OutletNodeIds = branchConfig.OutletNodeIds
            };
            solver.AddBranchPoint(branch);
        }
        
        // Add boundary conditions
        foreach (var bcConfig in config.BoundaryConditions)
        {
            var bc = CreateBoundaryCondition(bcConfig);
            if (bc != null)
                solver.AddBoundaryCondition(bc);
        }
        
        return solver;
    }
    
    private static BoundaryCondition? CreateBoundaryCondition(BoundaryConditionConfig config)
    {
        switch (config.Type.ToLower())
        {
            case "flowrate":
                return CreateFlowRateBC(config);
            case "pressure":
                return CreatePressureBC(config);
            case "resistance":
                return CreateResistanceBC(config);
            default:
                return null;
        }
    }

    private static FlowRateBoundaryCondition? CreateFlowRateBC(BoundaryConditionConfig config)
    {
        if (config.Parameters is Newtonsoft.Json.Linq.JObject jObj)
        {
            var type = jObj["type"]?.ToString() ?? "constant";
            var value = jObj["value"]?.ToObject<double>() ?? 0.0;
            
            Func<double, double> flowFunc = type switch
            {
                "constant" => t => value,
                "sine" => t => value * Math.Sin(2 * Math.PI * t),
                "cardiac" => t => CreateCardiacWaveform(t, value),
                _ => t => value
            };
            
            return new FlowRateBoundaryCondition(config.NodeId, flowFunc);
        }
        return null;
    }
    
    private static PressureBoundaryCondition? CreatePressureBC(BoundaryConditionConfig config)
    {
        if (config.Parameters is Newtonsoft.Json.Linq.JObject jObj)
        {
            var value = jObj["value"]?.ToObject<double>() ?? 100.0;
            return new PressureBoundaryCondition(config.NodeId, t => value);
        }
        return null;
    }
    
    private static ResistanceBoundaryCondition? CreateResistanceBC(BoundaryConditionConfig config)
    {
        if (config.Parameters is Newtonsoft.Json.Linq.JObject jObj)
        {
            var resistance = jObj["resistance"]?.ToObject<double>() ?? 1.0;
            return new ResistanceBoundaryCondition(config.NodeId, resistance);
        }
        return null;
    }
    
    private static double CreateCardiacWaveform(double t, double meanFlow)
    {
        // Simple cardiac waveform approximation
        double period = 1.0; // 1 second cardiac cycle
        double phase = (t % period) / period;
        
        if (phase < 0.3) // Systole
            return meanFlow * (1.0 + 2.0 * Math.Sin(Math.PI * phase / 0.3));
        else // Diastole
            return meanFlow * 0.5;
    }
}
