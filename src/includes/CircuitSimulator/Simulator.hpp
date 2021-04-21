#ifndef _SIMULATOR_HPP_INC_
#define _SIMULATOR_HPP_INC_

#ifdef WITH_PYTHON
#include "matplotlibcpp.h"
#endif

#include "CircuitElements/CircuitElements.hpp"
#include "Maths/DynamicMatrix.hpp"

#ifdef WITH_MATLAB
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#endif


#include <string>
#include <regex>
#include <iostream>
#include <fstream>

/// @brief The first character of each line for each component type
enum class LineType {
    Resistor = 'R',
    Capacitor = 'C',
    Inductor = 'L',
    CurrentSource = 'I',
    VoltageSource = 'V',
    SParameterBlock = 'S',
    Transistor = 'Q',
    Diode = 'D',
    Comment = '%',
    Directive = '.',
};

/// @brief The type of (voltage) source
enum class SourceType {
    TimeSeries = 'T',
    Sinusoidal = 'S',
};

/// @brief The main class to hold all of the relevant simulation data
///
/// @tparam VT The type used for values. e.g. double, float, etc
template<typename VT>
class SimulationEnvironment {
public:
    /// @brief Setup the simulation environment from a netlist
    ///
    /// @param netlistPath Path to the netlist that is going to be used in the
    /// simulation
    SimulationEnvironment(std::string netlistPath)
        : netlistPath(netlistPath), luPair(0), scratchSpace(0, 0),
          solutionMat(0, 0) {
#ifdef WITH_MATLAB
        matlabDesktop = matlab::engine::findMATLAB().size() > 0;
        matlabEngine = matlab::engine::connectMATLAB();
        elements.matlabEngine = matlabEngine;
#endif
        std::ifstream netlist(netlistPath);

        std::string line;

        std::regex transientRegex(R"(^\.transient\((.+?),(.+?),(.+?)\)\s?$)");
        std::regex graphRegex(R"(^\.graph\((.+?)\)\s?$)");
        std::regex noDCRegex(R"(^\.nodc\s?$)");
        std::regex outputFileRegex(R"(^\.outputFile\(\s*['"](.+?)['"]\s*\)\s?$)");
        std::smatch matches;

        while (!netlist.eof()) {
            std::getline(netlist, line);

            LineType lType = static_cast<LineType>(line[0]);

            switch (lType) {
                case LineType::Resistor:
                    Resistor<VT>::addToElements(line, elements, numNodes,
                                                numCurrents, numDCCurrents);
                    break;
                case LineType::Capacitor:
                    if (line[1] == 'N') {
                        NLCapacitor<VT>::addToElements(line, elements, numNodes,
                                                       numCurrents, numDCCurrents);
                    } else {
                        Capacitor<VT>::addToElements(line, elements, numNodes,
                                                     numCurrents, numDCCurrents);
                    }
                    break;
                case LineType::Inductor:
                    Inductor<VT>::addToElements(line, elements, numNodes,
                                                numCurrents, numDCCurrents);
                    break;
                case LineType::CurrentSource:
                    if (line[1] == 'N') {
                        NLCurrentSource<VT>::addToElements(line, elements, numNodes,
                                                           numCurrents,
                                                           numDCCurrents);
                    } else {
                        CurrentSource<VT>::addToElements(line, elements, numNodes,
                                                         numCurrents, numDCCurrents);
                    }
                    break;
                case LineType::VoltageSource:
                    switch (static_cast<SourceType>(line[1])) {
                        case SourceType::TimeSeries:
                            TimeSeriesVoltageSource<
                                VT>::addToElements(line, elements, numNodes,
                                                   numCurrents, numDCCurrents);
                            break;
                        case SourceType::Sinusoidal:
                            SinusoidalVoltageSource<
                                VT>::addToElements(line, elements, numNodes,
                                                   numCurrents, numDCCurrents);
                            break;
                        default:
                            VoltageSource<VT>::addToElements(line, elements,
                                                             numNodes, numCurrents,
                                                             numDCCurrents);
                            break;
                    }
                    break;
                case LineType::SParameterBlock:
                    if (line[1] == 'V' && (line[2] == 'P' || line[2] == 'F')) {
#ifndef WITH_MATLAB
                        if (line[2] == 'F') {
                            std::throw(
                                "ERROR: Matlab not available at compile time");
                        }
#endif
                        SParameterBlockVF<VT>::addToElements(line, elements,
                                                             numNodes, numCurrents,
                                                             numDCCurrents);
                    } else {
                        SParameterBlock<VT>::addToElements(line, elements, numNodes,
                                                           numCurrents,
                                                           numDCCurrents);
                    }
                    break;
                case LineType::Transistor:
                    if (line[1] == 'N') {
                        BJTN<VT>::addToElements(line, elements, numNodes,
                                                numCurrents, numDCCurrents);
                    } else if (line[1] == 'P') {
                        BJTP<VT>::addToElements(line, elements, numNodes,
                                                numCurrents, numDCCurrents);
                    } else if (line[1] == 'M') {
                        if (line[2] == 'N') {
                            NLNMOS<VT>::addToElements(line, elements, numNodes,
                                                      numCurrents, numDCCurrents);
                        } else {
                            std::cout << "Other Transistors not implemented yet"
                                      << std::endl;
                        }
                    } else {
                        std::cout << "Other Transistors not implemented yet"
                                  << std::endl;
                    }
                    break;
                case LineType::Diode:
                    Diode<VT>::addToElements(line, elements, numNodes, numCurrents,
                                             numDCCurrents);
                    break;
                case LineType::Comment:
                    break;
                case LineType::Directive:
                    std::regex_match(line, matches, transientRegex);

                    if (matches.size()) {
                        if constexpr (std::is_same_v<VT, double> ||
                                      std::is_same_v<VT, float>) {
                            initialTime = std::stod(matches.str(1));
                            timestep = std::stod(matches.str(3));
                            finalTime = std::stod(matches.str(2));
                            steps = (finalTime - initialTime) / timestep;
                        } else {
                            static_assert("Unsupported Type");
                        }
                        break;
                    }

                    std::regex_match(line, matches, graphRegex);
                    if (matches.size()) {
                        parseGraph(matches.str(1));
                        break;
                    }

                    std::regex_match(line, matches, noDCRegex);
                    if (matches.size()) {
                        performDCAnalysis = false;
                        break;
                    }

                    std::regex_match(line, matches, outputFileRegex);
                    if (matches.size()) {
                        outputFilePath = matches.str(1);
                        break;
                    }

                    std::cout << "Unsupported Directive" << std::endl;

                    break;
            }
        }

        elements.setNewStampSize(numNodes, numCurrents, numDCCurrents);

        size_t sizeMat = elements.staticStamp.G.M;

        solutionMat = Matrix<VT>(sizeMat, steps, 0);

        luPair = LUPair<VT>(sizeMat);
        scratchSpace = Matrix<VT>(sizeMat, 1);

        for (auto & comp : elements.staticElements) {
            comp->setTimestep(timestep);
        }
        for (auto & comp : elements.dynamicElements) {
            comp->setTimestep(timestep);
        }
        for (auto & comp : elements.nonLinearElements) {
            comp->setTimestep(timestep);
        }

        if (performDCAnalysis) {
            setDCOpPoint();
        }
    }

    /// @brief a function to determine and set the DC operating point
    void setDCOpPoint() {
        Matrix<VT> dcSoln = Matrix<VT>(solutionMat.M + numDCCurrents, 1);
        Matrix<VT> scratchSpace = Matrix<VT>(solutionMat.M + numDCCurrents, 1);
        LUPair<VT> luPair = LUPair<VT>(solutionMat.M + numDCCurrents);

        auto simStartTime = std::chrono::high_resolution_clock::now();
        for (size_t nr = 0; nr < 35; nr++) {
            auto & stamp = elements.generateDCStamp(dcSoln, numCurrents);
            stamp.G.luPair(luPair);
            stamp.G.leftDivide(stamp.s, luPair, scratchSpace, dcSoln);
        }

        for (size_t k = 0; k < solutionMat.M; k++) {
            solutionMat(k, 0) = dcSoln(k, 0);
        }

        elements.updateDCStoredState(dcSoln, numCurrents);

        auto simEndTime = std::chrono::high_resolution_clock::now();
        auto timeTaken = std::chrono::duration_cast<std::chrono::nanoseconds>(
                             simEndTime - simStartTime)
                             .count();


        std::cout << "DC OP in: " << timeTaken * 1e-6 << " ms (" << timeTaken
                  << " ns)" << std::endl;
    }

    /// @brief Where a lot of the magic starts. This is what runs the simulation
    ///
    /// @details This function is a simple newton-raphson solver, it still has room
    ///          for improvement, such as early termination when the loop converges.
    ///          \n\n
    ///          After the simulation has run to completion, the raw data is dumped,
    ///          and any graphs that were due to be generated are created.
    void simulate() {
        constexpr VT convergedThreshold = 1e-12;
        constexpr size_t maxNR = 32;
        Matrix<VT> tempSoln = Matrix<VT>(solutionMat.M, 1);
        VT maxDiff;
        VT singleVarDiff;
        auto simStartTime = std::chrono::high_resolution_clock::now();
        for (size_t n = 1; n < steps; n++) {
            size_t nr;
            for (nr = 0; nr < maxNR; nr++) {
                auto & stamp = elements.generateNonLinearStamp(solutionMat, n,
                                                               timestep);
                stamp.G.luPair(luPair);
                stamp.G.leftDivide(stamp.s, luPair, scratchSpace, tempSoln);

                maxDiff = 0;
                for (size_t k = 0; k < solutionMat.M; k++) {
                    singleVarDiff = std::abs(solutionMat(k, n) - tempSoln(k, 0));
                    maxDiff = std::max(maxDiff, singleVarDiff);
                }

                for (size_t k = 0; k < solutionMat.M; k++) {
#ifdef _DEBUG
                    if (std::isnan(tempSoln(k, 0))) {
                        std::cout << "Simulation Error: NaN found in solution"
                                  << std::endl;
                    }
#endif
                    solutionMat(k, n) = tempSoln(k, 0);
                }
                if (maxDiff < convergedThreshold) {
                    break;
                }
                elements.nonLinearStampIsFresh = false;
            }

#ifdef _DEBUG
            if (nr < maxNR) {
                std::cout << "NR terminated at: " << nr << " steps" << std::endl;
            }
#endif

            elements.updateTimeStep(solutionMat, n, timestep);
            if (n == 1) {
                elements.staticStampIsFresh = false; // for VF s-param model update
            }
        }
        auto simEndTime = std::chrono::high_resolution_clock::now();
        auto timeTaken = std::chrono::duration_cast<std::chrono::nanoseconds>(
                             simEndTime - simStartTime)
                             .count();
        // std::cout << "Time taken for simulation: " << timeTaken << std::endl;
        std::cout << timeTaken * 1e-6 << " ms (" << timeTaken << " ns)" << std::endl;
        std::ofstream runtimeFile("RunTimes.txt", std::ofstream::app);
        runtimeFile << netlistPath << " " << timeTaken << std::endl;

        dataDump();

        size_t graphNum = 1;
        for (auto nodes : nodesToGraph) {
            printMultipleOnGraph(nodes, std::to_string(graphNum++));
        }
    }

    /// @brief Outputs a single node's (or current's) time series to a graph, saving
    ///        as both eps and png.
    ///
    /// @param node The index of the node to plot
    void printGraph(size_t node) {
#ifdef WITH_PYTHON
        namespace plt = matplotlibcpp;
        assert(node > 0);
        std::vector<double> timeVector(steps);
        std::vector<double> voltageNode(steps);
        for (size_t n = 0; n < steps; n++) {
            voltageNode[n] = solutionMat(node - 1, n);
            timeVector[n] = n * timestep;
        }
        // Set the size of output image to 1200x780 pixels
        plt::figure_size(1200, 780);
        // Plot line from given x and y data. Color is selected automatically.
        plt::plot(timeVector, voltageNode);

        std::string filename = "Node " + std::to_string(node) + ".png";
        std::string filenameEps = "Node " + std::to_string(node) + ".eps";
        // std::cout << "Saving result to " << filename << std::endl;
        plt::save(filename.c_str());
        plt::save(filenameEps.c_str());
        plt::close();
#endif
    }

    /// @brief Similar to printGraph, but instead can plot multiple series on the
    ///        same graph.
    ///
    /// @param nodeVec A vector of node indices
    /// @param suffix What's appended to "Graph" to make the name of the output file.
    void printMultipleOnGraph(std::vector<size_t> nodeVec, std::string suffix = "") {
#ifdef WITH_PYTHON
        namespace plt = matplotlibcpp;
        plt::figure_size(1200, 780);
        for (auto node : nodeVec) {
            assert(node > 0);
            std::vector<double> timeVector(steps);
            std::vector<double> voltageNode(steps);
            for (size_t n = 0; n < steps; n++) {
                voltageNode[n] = solutionMat(node - 1, n);
                timeVector[n] = n * timestep;
            }
            // Set the size of output image to 1200x780 pixels
            // Plot line from given x and y data. Color is selected automatically.
            std::map<std::string, std::string> kwArgs;
            kwArgs["label"] = "Node " + std::to_string(node);
            plt::plot(timeVector, voltageNode, kwArgs);
        }
        plt::legend();
        std::string name = "Graph";
        name += suffix;

        plt::save(name + ".png");
        plt::save(name + ".eps");
        plt::close();
#endif
    }

    /// @brief Simple function to dump the output of the simulator in a matlab
    ///        table readable format
    void dataDump() {
#ifdef WITH_MATLAB
        // Create matlab data array factory
        matlab::data::ArrayFactory factory;

        std::vector<std::string> varNames = {"t"};
#endif

        std::ofstream outputFile(outputFilePath);
        outputFile << "time";
        for (int i = 1; i <= numNodes; i++) {
            outputFile << "\t"
                       << "n" << i;
#ifdef WITH_MATLAB
            if (matlabDesktop) {
                varNames.emplace_back(std::string("n") + std::to_string(i));
            }
#endif
        }

        for (int i = 1; i <= numCurrents; i++) {
            outputFile << "\t"
                       << "i" << i;
#ifdef WITH_MATLAB
            if (matlabDesktop) {
                varNames.emplace_back(std::string("i") + std::to_string(i));
            }
#endif
        }

#ifdef WITH_MATLAB
        auto sArray = factory.createStructArray({solutionMat.N, 1}, varNames);
#endif
        for (size_t n = 0; n < solutionMat.N; n++) {
            outputFile << std::endl;
            outputFile << std::setprecision(9) << n * timestep;
#ifdef WITH_MATLAB
            if (matlabDesktop) {
                sArray[n][varNames[0]] = factory.createArray({1, 1}, {n * timestep});
            }
#endif
            for (int i = 0; i < numNodes + numCurrents; i++) {
                outputFile << "\t" << std::setprecision(9) << solutionMat(i, n);
#ifdef WITH_MATLAB
                if (matlabDesktop) {
                    sArray[n][varNames[i + 1]] = factory
                                                     .createArray({1, 1},
                                                                  {solutionMat(i,
                                                                               n)});
                }
#endif
            }
        }
        outputFile.close();

#ifdef WITH_MATLAB
        if (matlabDesktop) {
            matlabEngine->setVariable(u"solutionData", sArray,
                                      matlab::engine::WorkspaceType::GLOBAL);
        }
#endif
    }


private:
#ifdef WITH_MATLAB
    std::shared_ptr<matlab::engine::MATLABEngine> matlabEngine;
    bool matlabDesktop = false;
#endif


    /// @brief Helper function to pull indices from graph netlist directive
    ///
    /// @param line The line to parse
    void parseGraph(std::string line) {
        std::stringstream sstr(line);
        nodesToGraph.emplace_back(std::vector<size_t>());
        std::vector<size_t> & toGraph = nodesToGraph.back();
        size_t nodeNum = 0;
        sstr >> nodeNum;

        while (!sstr.fail()) {
            toGraph.emplace_back(nodeNum);
            sstr.ignore(1);
            sstr >> nodeNum;
        }
    }

    std::string outputFilePath = "datadump.txt";
    std::string netlistPath = "";

    double initialTime;
    double timestep;
    double finalTime;
    size_t steps;

    size_t numNodes = 1;
    size_t numCurrents = 0;
    size_t numDCCurrents = 0;
    bool performDCAnalysis = true;
    /// @brief A collection of all the circuit elements
    CircuitElements<VT> elements;

    /// @brief Preallocated space to prevent repeated allocations and deallocations
    LUPair<VT> luPair;
    /// @brief Preallocated space to prevent repeated allocations and deallocations.
    ///        using during the leftDivide to solve the MNA system
    Matrix<VT> scratchSpace;

    /// @brief Keeps track of the nodes to be graphed after simulation
    std::vector<std::vector<size_t> > nodesToGraph;

    /// @brief Preallocated space to store the results in
    Matrix<VT> solutionMat;
};

#endif
