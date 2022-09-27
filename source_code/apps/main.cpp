#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <numeric>

#include <boost/program_options.hpp>


#include "optimizer/Optimizer.hpp"


void declareOptions(boost::program_options::options_description& desc)
{
    desc.add_options()
        // General Settings:
        ("help", "produce help message")
        ("logs", boost::program_options::bool_switch()->default_value(false),
            "Write program logs.")

        // File Paths and settings:
        ("nodeCityFile", boost::program_options::value<std::string>(),
            "Path to input file with node city pair data.")
        ("edgesFile", boost::program_options::value<std::string>(),
            "Path to input file with edge data.")
        ("citiesFile", boost::program_options::value<std::string>(),
            "Path to input file with city data.")
        ("addedEdgesFile", boost::program_options::value<std::string>(),
            "Path to output file with added edges data.")
        ("preparedNodesFile", boost::program_options::value<std::string>(),
            "Path to output file with prepared combined nodes city data.")
        ("preparedEdgesFile", boost::program_options::value<std::string>(),
            "Path to output file with prepared edges data.")
        ("modelFile", boost::program_options::value<std::string>(),
            "Path to output file with best model stations representation.")
        ("savePrepared", boost::program_options::bool_switch()->default_value(false),
            "Save preprocessed graph in the given file.")
        ("saveBestModel", boost::program_options::bool_switch()->default_value(false),
            "Save best model stations to the given file.")

        // Simulation Parameters:
        ("simulationTime", boost::program_options::value<int>()->default_value(5000),
            "Time of the simulation in minutes.")
        ("segmentLength", boost::program_options::value<double>()->default_value(1),
            "Length of the edge segments in kilometers.")
        ("numClosestStations", boost::program_options::value<int>()->default_value(1),
            "How many nearest charging stations consider during finding the best one to visit.")
        ("carConsumption", boost::program_options::value<double>()->default_value(0.002),
            "Consumption of the car battery level per minute.")
        ("numStations", boost::program_options::value<int>()->default_value(300),
            "Total number of charging stations in the map.")
        ("stationCapacity", boost::program_options::value<int>()->default_value(4),
            "Capacity of the charging station (how many vehicles it can serve at once).")
        ("exponentialLambdaCities", boost::program_options::value<double>()->default_value(0.01),
            "Parameter of the exponential disitribution for generating final city based \
            on the traveled distance (mean distance in kilometers is (1/exponentialLambdaCities).")
        ("exponentialLambdaDepartures", boost::program_options::value<double>()->default_value(10),
            "Parameter of the exponential distribution of next car departure time \
            (mean in minutes is (1/exponentialLambdaDepartures)).")
        ("endCityRatio", boost::program_options::value<double>()->default_value(1),
            "How much important is the distance of the final city in comparison with city population \
            (it is the parameter in linear combination of distance and population likelihood). \
            If `< 1` - distance is less important, `== 1` - both equaly important, \
            `> 1' - distance is more important).")
        ("batteryTresholdLambda", boost::program_options::value<double>()->default_value(20),
            "Parameter of the exponential distribution which is used to stochasticaly decide \
            whether go to the charging station or not (used when battery level between top and \
            bottom tresholds (when not decision to go charging is not certain). Value of the parameter \
            means that mean battery level to go to the charging station is \
            (1/batteryTresholdLambda) more than the bottom level.")
        ("carBatteryMean", boost::program_options::value<double>()->default_value(0.9),
            "Mean of the battery level at the start.")
        ("carBatteryDeviation", boost::program_options::value<double>()->default_value(0.2),
            "Standard deviation of the battery level at the start.")
        ("carStartBatteryBottomLimit", boost::program_options::value<double>()->default_value(0.5),
            "Minimal battery level of the car at the begining.")
        ("chargingTreshold", boost::program_options::value<double>()->default_value(0.5),
            "If battery level smaller than given theshold, then always consider going charging.")
        ("notChargingTreshold", boost::program_options::value<double>()->default_value(0.9),
            "If battery level above the given treshold, then never consider going charging.")
        ("batteryTolerance", boost::program_options::value<double>()->default_value(0.05),
            "The minimal expected battery level in the end of the path \
            when it is possible not to go charging (based on the current traffic).")
        ("carVelocity", boost::program_options::value<double>()->default_value(1),
            "Average speed of the car when there is no traffic (in kilometers per minute (1/60 * km/hr)")
        ("chargingWaitingTime", boost::program_options::value<double>()->default_value(40),
            "Waiting time in minutes to fully charge car battery (from level 0 to level 1).")
        ("meanChargingLevel", boost::program_options::value<double>()->default_value(0.5),
            "Mean battery level to be charged (to estimate waiting time in the queue).")
        
        // Optimizer Parameters:
        ("lossDifferenceTreshold", boost::program_options::value<double>()->default_value(10),
            "Minimal difference of loss to still continue with iteration (for some optimization algorithms).")
        ("stationNumberParameter", boost::program_options::value<double>()->default_value(100),
            "Multiplication parameter for number of stations for loss computing.")
        ("runDownParameter", boost::program_options::value<double>()->default_value(100),
            "Multiplication parameter for number of run down battery cars for loss computing.")
        ("durationParameter", boost::program_options::value<double>()->default_value(0.001),
            "Multiplication parameter for average traveling duration for loss computing.")
        ("batteryDifferenceParameter", boost::program_options::value<double>()->default_value(1),
            "Multiplication parameter for average battery difference between start and end of the route \
            for loss computing.")
        ("waitingTimesParameter", boost::program_options::value<double>()->default_value(1),
            "Multiplication parameter for average waiting times in the charging stations for loss computing.")

        // Greedy Parameters:
        ("greedyMaxIterations", boost::program_options::value<int>()->default_value(5),
            "Maximal number of iterations for the greedy algorithm.")
        ("greedyNumThrowAway", boost::program_options::value<int>()->default_value(5),
            "Maximal number of stations to throw away between greedy iterations.")

        // Genetic Parameters:
        ("geneticPopulatioSize", boost::program_options::value<int>()->default_value(10),
            "Size of the population for the genetic algorithm optimalization.")
        ("geneticNumGenerations", boost::program_options::value<int>()->default_value(10),
            "Number of generations of the genetic algorithm optimalization.")
        ("geneticNumBestSelection", boost::program_options::value<int>()->default_value(2),
            "Number of best individuals to immediately copy to next generation \
            for the genetic algorithm optimalization.")
        ("geneticTournamentSelectionTreshold", boost::program_options::value<double>()->default_value(0.8),
            "Probability of choosing the better individual during the tournament selection\
             in the genetic algorithm optimalization.")
        ("geneticMutationTreshold", boost::program_options::value<double>()->default_value(0.01),
            "Probability of the mutation of one station in the individual in \
            the genetic algorithm optimalization.")
        ("geneticMemberSizeVariance", boost::program_options::value<double>()->default_value(2),
            "Variance of the member size for normal distribution (mean is current population size).")

        // K-Means Parameters:
        ("kMeansNumIterationsOneRun", boost::program_options::value<int>()->default_value(50),
            "Number of iterations of one run of K-Means algorithm \
            (after how many steps we expect convergence).")
        ("kMeansNumGenerations", boost::program_options::value<int>()->default_value(5),
            "Number of generations of the K-Means algorithm \
            (how many times we aply K-Means centroids as station aspirants)")

        // Algorithm choice:
        ("randomModels", boost::program_options::bool_switch()->default_value(false),
            "Perform simulation on randomly generated stations.")
        ("greedy", boost::program_options::bool_switch()->default_value(false),
            "Perform greedy optimalization algorithm.")
        ("genetic", boost::program_options::bool_switch()->default_value(false),
            "Perform genetic algorithm optimization.")
        ("kMeans", boost::program_options::bool_switch()->default_value(false),
            "Perform K-Means algorithm for the optimization.")
        ;
}


void parseFilenames(boost::program_options::variables_map& vm,
    SimulationParameters& simulationParameters,
    std::string& edgesFile,
    std::string& nodeCityFile,
    std::string& citiesFile,
    std::string& addedEdgesFile,
    std::string& preparedNodesFile,
    std::string& preparedEdgesFile)
{
    simulationParameters.EdgesFile = edgesFile;
    simulationParameters.NodeCityFile = nodeCityFile;
    simulationParameters.CitiesFile = citiesFile;
    simulationParameters.AddedEdgesFile = addedEdgesFile;
    simulationParameters.PreparedNodesFile = preparedNodesFile;
    simulationParameters.PreparedEdgesFile = preparedEdgesFile;

    if (vm.count("savePrepared"))
    {
        simulationParameters.SavePrepared = vm["savePrepared"].as<bool>();
        std::cout << "Save prepared is: " << simulationParameters.SavePrepared << std::endl;
    }

    if (vm.count("edgesFile"))
    {
        simulationParameters.EdgesFile = vm["edgesFile"].as<std::string>();
        std::cout << "Edges: " << simulationParameters.EdgesFile << std::endl;
    }

    if (vm.count("nodeCityFile"))
    {
        simulationParameters.NodeCityFile = vm["nodeCityFile"].as<std::string>();
        std::cout << "Node City: " << simulationParameters.NodeCityFile << std::endl;
    }

    if (vm.count("citiesFile"))
    {
        simulationParameters.CitiesFile = vm["citiesFile"].as<std::string>();
        std::cout << "Cities: " << simulationParameters.CitiesFile << std::endl;
    }

    if (vm.count("addedEdgesFile"))
    {
        simulationParameters.AddedEdgesFile = vm["addedEdgesFile"].as<std::string>();
        std::cout << "Added edges: " << simulationParameters.AddedEdgesFile << std::endl;
    }

    if (vm.count("preparedNodesFile"))
    {
        simulationParameters.PreparedNodesFile = vm["preparedNodesFile"].as<std::string>();
        std::cout << "Prepared nodes: " << simulationParameters.PreparedNodesFile << std::endl;
    }

    if (vm.count("modelFile"))
    {
        std::cout << "Best model: " << vm.count("modelFile") << std::endl;
    }

    if (vm.count("preparedEdgesFile"))
    {
        simulationParameters.PreparedEdgesFile = vm["preparedEdgesFile"].as<std::string>();
        std::cout << "Prepared edges: " << simulationParameters.PreparedEdgesFile << std::endl;
    }
}

void parseSimulatorArgumets(boost::program_options::variables_map& vm, SimulationParameters& simulationParameters)
{
    if (vm.count("simulationTime"))
    {
        simulationParameters.SimulationTime = vm["simulationTime"].as<int>();
        std::cout << "Simulation time: " << simulationParameters.SimulationTime << std::endl;
    }

    if (vm.count("segmentLength"))
    {
        simulationParameters.SegmentLength = vm["segmentLength"].as<double>();
        std::cout << "Length of the edge segment: " << simulationParameters.SegmentLength << std::endl;
    }

    if (vm.count("numClosestStations"))
    {
        simulationParameters.NumClosestStations = vm["numClosestStations"].as<int>();
        std::cout << "Number of closest stations: " << simulationParameters.NumClosestStations << std::endl;
    }
    if (vm.count("carConsumption"))
    {
        simulationParameters.CarConsumption = vm["carConsumption"].as<double>();
        std::cout << "Car consumtion: " << simulationParameters.CarConsumption << std::endl;
    }
    if (vm.count("numStations"))
    {
        simulationParameters.NumStations = vm["numStations"].as<int>();
        std::cout << "Number of charging stations: " << simulationParameters.NumStations << std::endl;
    }
    if (vm.count("stationCapacity"))
    {
        simulationParameters.StationCapacity = vm["stationCapacity"].as<int>();
        std::cout << "Capacity of the charging station: " << simulationParameters.StationCapacity << std::endl;
    }
    if (vm.count("exponentialLambdaCities"))
    {
        simulationParameters.ExponentialLambdaCities = vm["exponentialLambdaCities"].as<double>();
        std::cout << "Exponential parameter for city distances: " << simulationParameters.ExponentialLambdaCities << std::endl;
    }
    if (vm.count("exponentialLambdaDepartures"))
    {
        simulationParameters.ExponentialLambdaDepartures = vm["exponentialLambdaDepartures"].as<double>();
        std::cout << "Exponential parameter for car departures: " << simulationParameters.ExponentialLambdaDepartures << std::endl;
    }
    if (vm.count("endCityRatio"))
    {
        simulationParameters.EndCityRatio = vm["endCityRatio"].as<double>();
        std::cout << "End city distance versus population ratio: " << simulationParameters.EndCityRatio << std::endl;
    }
    if (vm.count("batteryTresholdLambda"))
    {
        simulationParameters.BatteryTresholdLambda = vm["batteryTresholdLambda"].as<double>();
        std::cout << "Exponential parameter for battery treshold: " << simulationParameters.BatteryTresholdLambda << std::endl;
    }
    if (vm.count("carBatteryMean"))
    {
        simulationParameters.CarBatteryMean = vm["carBatteryMean"].as<double>();
        std::cout << "Car battery level mean: " << simulationParameters.CarBatteryMean << std::endl;
    }
    if (vm.count("carBatteryDeviation"))
    {
        simulationParameters.CarBatteryDeviation = vm["carBatteryDeviation"].as<double>();
        std::cout << "Car battery level deviation: " << simulationParameters.CarBatteryDeviation << std::endl;
    }
    if (vm.count("carStartBatteryBottomLimit"))
    {
        simulationParameters.CarStartBatteryBottomLimit = vm["carStartBatteryBottomLimit"].as<double>();
        std::cout << "Car start battery level bottom limit: " << simulationParameters.CarStartBatteryBottomLimit << std::endl;
    }
    if (vm.count("chargingTreshold"))
    {
        simulationParameters.ChargingTreshold = vm["chargingTreshold"].as<double>();
        std::cout << "Battery treshold to always charge: " << simulationParameters.ChargingTreshold << std::endl;
    }
    if (vm.count("notChargingTreshold"))
    {
        simulationParameters.NotChargingTreshold = vm["notChargingTreshold"].as<double>();
        std::cout << "Battery treshold to never go charging: " << simulationParameters.NotChargingTreshold << std::endl;
    }
    if (vm.count("batteryTolerance"))
    {
        simulationParameters.BatteryTolerance = vm["batteryTolerance"].as<double>();
        std::cout << "Battery tolerance level in the end of path: " << simulationParameters.BatteryTolerance << std::endl;
    }
    if (vm.count("carVelocity"))
    {
        simulationParameters.CarVelocity = vm["carVelocity"].as<double>();
        std::cout << "Car velocity: " << simulationParameters.CarVelocity << std::endl;
    }
    if (vm.count("chargingWaitingTime"))
    {
        simulationParameters.ChargingWaitingTime = vm["chargingWaitingTime"].as<double>();
        std::cout << "Waiting time of fully charge: " << simulationParameters.ChargingWaitingTime << std::endl;
    }
    if (vm.count("meanChargingLevel"))
    {
        simulationParameters.MeanChargingLevel = vm["meanChargingLevel"].as<double>();
        std::cout << "Estimated mean battery level to charge: " << simulationParameters.MeanChargingLevel << std::endl;
    }
}


void parseOptimizerArgumets(boost::program_options::variables_map& vm, OptimizerParameters& optimizerParameters)
{
    if (vm.count("lossDifferenceTreshold"))
    {
        optimizerParameters.ScoreDifferenceTreshold_ = vm["lossDifferenceTreshold"].as<double>();
        std::cout << "Loss difference treshold: " << optimizerParameters.ScoreDifferenceTreshold_ << std::endl;
    }
    if (vm.count("stationNumberParameter"))
    {
        optimizerParameters.StationNumberParameter_ = vm["stationNumberParameter"].as<double>();
        std::cout << "Station number loss parameter: " << optimizerParameters.StationNumberParameter_ << std::endl;
    }
    if (vm.count("runDownParameter"))
    {
        optimizerParameters.RunDownParameter_ = vm["runDownParameter"].as<double>();
        std::cout << "Battery run down loss parameter: " << optimizerParameters.RunDownParameter_ << std::endl;
    }
    if (vm.count("durationParameter"))
    {
        optimizerParameters.DurationParameter_ = vm["durationParameter"].as<double>();
        std::cout << "Travel duration loss parameter: " << optimizerParameters.DurationParameter_ << std::endl;
    }
    if (vm.count("batteryDifferenceParameter"))
    {
        optimizerParameters.BatteryDifferenceParameter_ = vm["batteryDifferenceParameter"].as<double>();
        std::cout << "Battery difference loss paramter: " << optimizerParameters.BatteryDifferenceParameter_ << std::endl;
    }
    if (vm.count("waitingTimesParameter"))
    {
        optimizerParameters.WaitingTimesParameter_ = vm["waitingTimesParameter"].as<double>();
        std::cout << "Station waiting time paramter: " << optimizerParameters.WaitingTimesParameter_ << std::endl;
    }
    if (vm.count("greedyMaxIterations"))
    {
        optimizerParameters.GreedyParameters_.MaxIterations_ = vm["greedyMaxIterations"].as<int>();
        std::cout << "Maximal number of iteration for greedy optimization: " 
            << optimizerParameters.GreedyParameters_.MaxIterations_ << std::endl;
    }
    if (vm.count("greedyNumThrowAway"))
    {
        optimizerParameters.GreedyParameters_.NumThrowAway_ = vm["greedyNumThrowAway"].as<int>();
        std::cout << "Number of station to throw away in greedy optimization: " << optimizerParameters.GreedyParameters_.NumThrowAway_ << std::endl;
    }
    if (vm.count("geneticPopulatioSize"))
    {
        optimizerParameters.GeneticParameters_.PopulationSize_ = vm["geneticPopulatioSize"].as<int>();
        std::cout << "Size of the population in genetic algorithm: " 
            << optimizerParameters.GeneticParameters_.PopulationSize_ << std::endl;
    }
    if (vm.count("geneticNumGenerations"))
    {
        optimizerParameters.GeneticParameters_.NumGenerations_ = vm["geneticNumGenerations"].as<int>();
        std::cout << "Number of generations in genetic algorithm: " 
            << optimizerParameters.GeneticParameters_.NumGenerations_ << std::endl;
    }
    if (vm.count("geneticNumBestSelection"))
    {
        optimizerParameters.GeneticParameters_.NumBestSelection_ = vm["geneticNumBestSelection"].as<int>();
        std::cout << "Number of best to choose automaticaly during selection in genetic algorithm: " 
            << optimizerParameters.GeneticParameters_.NumBestSelection_ << std::endl;
    }
    if (vm.count("geneticTournamentSelectionTreshold"))
    {
        optimizerParameters.GeneticParameters_.TournamentSelectionTreshold_ = vm["geneticTournamentSelectionTreshold"].as<double>();
        std::cout << "Tournament selection probability of the better model choice in genetic algorithm: "
            << optimizerParameters.GeneticParameters_.TournamentSelectionTreshold_ << std::endl;
    }
    if (vm.count("geneticMutationTreshold"))
    {
        optimizerParameters.GeneticParameters_.MutationTreshold_ = vm["geneticMutationTreshold"].as<double>();
        std::cout << "Probability of mutation of one station in genetic algorithm: " 
            << optimizerParameters.GeneticParameters_.MutationTreshold_ << std::endl;
    }
    if (vm.count("geneticMemberSizeVariance"))
    {
        optimizerParameters.GeneticParameters_.MemberSizeVariance_ = vm["geneticMemberSizeVariance"].as<double>();
        std::cout << "Variance of the population size in genetic algorithm: "
            << optimizerParameters.GeneticParameters_.MemberSizeVariance_ << std::endl;
    }
    if (vm.count("kMeansNumIterationsOneRun"))
    {
        optimizerParameters.KMeansParameters_.NumIterationsOneRun_ = vm["kMeansNumIterationsOneRun"].as<int>();
        std::cout << "Number of iteration during one K-Means run: " << optimizerParameters.KMeansParameters_.NumIterationsOneRun_ << std::endl;
    }
    if (vm.count("kMeansNumGenerations"))
    {
        optimizerParameters.KMeansParameters_.NumGenerations_ = vm["kMeansNumGenerations"].as<int>();
        std::cout << "Number of generations of K-Means stations: " << optimizerParameters.KMeansParameters_.NumGenerations_ << std::endl;
    }
}

void parseAlgorithmChoice(boost::program_options::variables_map& vm, Optimizer& optimizer, bool logs)
{
    if (vm.count("randomModels") && vm["randomModels"].as<bool>())
    {
        std::cout << "Running random models algorithm: " << std::endl;
        optimizer.RunMultipleSimulations(10, logs);   
    }
    if (vm.count("greedy") && vm["greedy"].as<bool>())
    {
        std::cout << "Running greedy algorith." << std::endl;
        optimizer.GreedyAlgorithm(logs);        
    }
    if (vm.count("genetic") && vm["genetic"].as<bool>())
    {
        std::cout << "Running genetic algorithm." << std::endl;
        optimizer.GeneticAlgorithm(logs);        
    }
    if (vm.count("kMeans") && vm["kMeans"].as<bool>())
    {
        std::cout << "Running K-Means algorithm." << std::endl;
        optimizer.KMeansAlgorithm(logs);      
    }
}



int main(int argc, char** argv)
{
    // Init default filenames:
    // std::string nodeCityFile = "prepared_graph/combined_nodes.txt";
    // std::string edgesFile = "prepared_graph/edges.txt";
    std::string nodeCityFile = "prepared_graph/prepared_combined_nodes.txt";
    std::string edgesFile = "prepared_graph/prepared_edges.txt";
    std::string citiesFile = "prepared_graph/cities.txt";
    std::string addedEdgesFile = "prepared_graph/added_edges.txt";
    std::string preparedNodesFile = "prepared_graph/prepared_combined_nodes.txt";
    std::string preparedEdgesFile = "prepared_graph/prepared_edges.txt";
    std::string bestModelFile = "models/best_model.txt";

    // Declare the supported options.
    boost::program_options::options_description desc("Allowed options");
    declareOptions(desc);

    
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help")) 
    {
        std::cout << desc << "\n";
        return 0;
    }

    bool logs;
    if (vm.count("logs"))
    {
        logs = vm["logs"].as<bool>();
    }

    SimulationParameters simulationParameters;
    OptimizerParameters optimizerParameters;

    parseFilenames(vm, simulationParameters, edgesFile, nodeCityFile, 
        citiesFile, addedEdgesFile, preparedNodesFile, preparedEdgesFile);
    parseSimulatorArgumets(vm, simulationParameters);
    parseOptimizerArgumets(vm, optimizerParameters);


    Optimizer optimizer(simulationParameters, optimizerParameters);

    parseAlgorithmChoice(vm, optimizer, logs);

    if (vm.count("saveBestModel"))
    {
        if (vm.count("modelFile"))
        {
            bestModelFile = vm["modelFile"].as<std::string>();
        }
        std::ofstream bestModelStream(bestModelFile);
        optimizer.SaveBestModel(bestModelStream);
    }


    //End of application
    return 0;
}
