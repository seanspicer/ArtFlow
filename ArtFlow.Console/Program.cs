using ArtFlow.Core;
using ArtFlow.Solver;
using ArtFlow.IO;

namespace ArtFlow.Console;

class Program
{
    static void Main(string[] args)
    {
        System.Console.WriteLine("ArtFlow - 1D Blood Flow Simulation");
        System.Console.WriteLine("==================================");
        
        if (args.Length < 1)
        {
            ShowUsage();
            return;
        }
        
        try
        {
            switch (args[0].ToLower())
            {
                case "run":
                    if (args.Length < 2)
                    {
                        System.Console.WriteLine("Error: Please specify a configuration file");
                        return;
                    }
                    RunSimulation(args[1]);
                    break;
                    
                case "example":
                    CreateExampleConfig();
                    break;
                    
                default:
                    ShowUsage();
                    break;
            }
        }
        catch (Exception ex)
        {
            System.Console.WriteLine($"Error: {ex.Message}");
        }
    }

    static void ShowUsage()
    {
        System.Console.WriteLine("Usage:");
        System.Console.WriteLine("  ArtFlow.Console run <config.json>  - Run simulation from config file");
        System.Console.WriteLine("  ArtFlow.Console example            - Create example configuration");
    }
    
    static void RunSimulation(string configFile)
    {
        System.Console.WriteLine($"Loading configuration from {configFile}...");
        var config = ConfigLoader.LoadFromFile(configFile);
        
        System.Console.WriteLine($"Building solver for simulation '{config.Name}'...");
        var solver = ConfigLoader.BuildSolver(config);
        
        System.Console.WriteLine("Running simulation...");
        solver.RunSimulation();
        
        System.Console.WriteLine("Simulation complete!");
    }
    
    static void CreateExampleConfig()
    {
        // Create example configuration for carotid artery
        var config = new SimulationConfig
        {
            Name = "Carotid Artery Example",
            SolverSettings = new SolverSettings
            {
                ElementSize = 1.0,
                TimeStepSize = 0.001,
                TotalTime = 3.0, // 3 cardiac cycles
                OutputInterval = 100
            },
            Segments = new List<SegmentConfig>
            {
                new SegmentConfig
                {
                    Id = 0,
                    Name = "Common Carotid",
                    Length = 16.0, // cm
                    InitialRadius = 0.4, // cm
                    InitialPressure = 100.0 // mmHg
                }
            },
            BoundaryConditions = new List<BoundaryConditionConfig>
            {
                new BoundaryConditionConfig
                {
                    NodeId = 0,
                    Type = "FlowRate",
                    Parameters = new { type = "cardiac", value = 5.0 } // ml/s
                },
                new BoundaryConditionConfig
                {
                    NodeId = 1, // Will be set after mesh generation
                    Type = "Resistance",
                    Parameters = new { resistance = 1.0 }
                }
            }
        };
        
        string filename = "carotid_example.json";
        ConfigLoader.SaveToFile(config, filename);
        System.Console.WriteLine($"Example configuration saved to {filename}");
    }
}
