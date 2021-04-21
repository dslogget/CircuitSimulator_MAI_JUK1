#ifndef _SPARAMETERBLOCKVF_HPP_INC_
#define _SPARAMETERBLOCKVF_HPP_INC_
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"
#include "Maths/dft.hpp"
#include "Maths/ForceCausal.hpp"

#ifdef WITH_MATLAB
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#endif

template<typename T>
struct SParamVFDataFrom {
    size_t numPoles = 0;

    std::vector<std::complex<T> > pole;
    std::vector<std::complex<T> > residue;
    std::complex<T> remainder = 0;

    /// @brief the per-pole contribution of the current awave
    std::vector<std::complex<T> > lambda_p;
    /// @brief the per-pole contribution of the previous awave
    std::vector<std::complex<T> > mu_p;
    /// @brief the per-pole contribution of the 2nd previous awave
    std::vector<std::complex<T> > nu_p;

    std::vector<std::complex<T> > exp_alpha;

    /// @brief the contribution of the current awave
    std::complex<T> lambda = 0;
    /// @brief the contribution of the previous awave
    std::complex<T> mu = 0;
    /// @brief the contribution of the 2nd previous awave
    std::complex<T> nu = 0;

    /// @brief The previous x values
    std::vector<std::complex<T> > x;
};

/// @brief a helper struct to store the information for the ports of an
///        S-Parameter Block
///
/// @tparam T The value type
template<typename T>
struct SParameterPortVF {
    /// @brief The positive index
    size_t positive = 0;
    /// @brief The negative index
    size_t negative = 0;
    /// @brief The current index
    size_t current = 0;


    /// @brief a constant factor equal to 1 / ( 1 - lambda )
    std::complex<T> beta = 0;

    /// @brief The equivalent scalar for controlled sources. Index is the second
    ///        port
    std::vector<std::complex<T> > alpha;

    /// @brief The equivalent resistance
    std::complex<T> R = 0;

    std::vector<SParamVFDataFrom<T> > from;
};


/// @brief A vectorfitting based model of an s-parameter block.
///
/// @tparam T
template<typename T>
struct SParameterBlockVF : public Component<T> {
    std::vector<SParameterPortVF<T> > port;
    size_t numPorts = 0;
    bool firstOrder = true;

    T z_ref = 0;

    std::complex<T> history_p(const size_t p, const Matrix<T> & solutionMatrix,
                              const size_t currentSolutionIndex, T timestep,
                              const size_t sizeG_A) const {
        std::complex<T> toRet = 0;
        for (size_t c = 0; c < numPorts; c++) {
            for (size_t rho = 0; rho < port[p].from[c].numPoles; rho++) {
                toRet += port[p].from[c].x[rho] * port[p].from[c].exp_alpha[rho];
            }

            toRet += port[p].from[c].mu *
                     awave_p(c, solutionMatrix, currentSolutionIndex - 1, sizeG_A);
            if (currentSolutionIndex > 1) {
                toRet += port[p].from[c].nu * awave_p(c, solutionMatrix,
                                                      currentSolutionIndex - 2,
                                                      sizeG_A);
            }
        }
        return 2.0 * toRet * std::sqrt(z_ref);
    }

    T
    V_p(const size_t p, const Matrix<T> & solutionMatrix,
        const size_t currentSolutionIndex, T timestep, const size_t sizeG_A) const {
        return std::real(
            history_p(p, solutionMatrix, currentSolutionIndex, timestep, sizeG_A) *
            port[p].beta);
    }

    T awave_p(const size_t p, const Matrix<T> & solutionMatrix,
              const size_t currentSolutionIndex, const size_t sizeG_A) const {
        size_t np = port[p].positive - 1;
        size_t nn = port[p].negative - 1;
        size_t curr = port[p].current - 1;
        size_t n = currentSolutionIndex;
        T toRet = 0;

        if (port[p].positive) {
            toRet += solutionMatrix(np, n);
        }

        if (port[p].negative) {
            toRet += -solutionMatrix(nn, n);
        }

        return (toRet + solutionMatrix(sizeG_A + curr, n) * z_ref) /
               (2 * std::sqrt(z_ref));
    }

    T bwave_p(const size_t p, const Matrix<T> & solutionMatrix,
              const size_t currentSolutionIndex, const size_t sizeG_A) const {
        size_t np = port[p].positive - 1;
        size_t nn = port[p].negative - 1;
        size_t curr = port[p].current - 1;
        size_t n = currentSolutionIndex;
        T toRet = 0;

        if (port[p].positive) {
            toRet += solutionMatrix(np, n);
        }

        if (port[p].negative) {
            toRet += -solutionMatrix(nn, n);
        }

        return (toRet - solutionMatrix(sizeG_A + curr, n) * z_ref) /
               (2 * std::sqrt(z_ref));
    }

    void addStaticStampTo(Stamp<T> & stamp) const {
        for (size_t p = 0; p < port.size(); p++) {
            size_t np = port[p].positive - 1;
            size_t nn = port[p].negative - 1;
            size_t curr = port[p].current - 1;
            // Voltage source and resistance
            //
            stamp.G(stamp.sizeG_A + curr,
                    stamp.sizeG_A + curr) += std::real(-port[p].R);

            if (port[p].positive != 0) {
                stamp.G(stamp.sizeG_A + curr, np) += 1;
                stamp.G(np, stamp.sizeG_A + curr) += 1;
            }

            if (port[p].negative != 0) {
                stamp.G(stamp.sizeG_A + curr, nn) += -1;
                stamp.G(nn, stamp.sizeG_A + curr) += -1;
            }
            // controlled sources
            for (size_t c = 0; c < port.size(); c++) {
                if (c != p) {
                    if (port[c].positive != 0) {
                        stamp.G(stamp.sizeG_A + curr,
                                port[c].positive -
                                    1) += std::real(-port[p].alpha[c]);
                    }
                    if (port[c].negative != 0) {
                        stamp.G(stamp.sizeG_A + curr,
                                port[c].negative - 1) += std::real(port[p].alpha[c]);
                    }
                    stamp.G(stamp.sizeG_A + curr,
                            stamp.sizeG_A + port[c].current -
                                1) += std::real(-z_ref * port[p].alpha[c]);
                }
            }
        }
    }

    void addDynamicStampTo(Stamp<T> & stamp, const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex,
                           T simulationTimestep) const {
        for (size_t p = 0; p < port.size(); p++) {
            size_t curr = port[p].current - 1;
            // V_p
            stamp.s(stamp.sizeG_A + curr,
                    0) += V_p(p, solutionMatrix, currentSolutionIndex,
                              simulationTimestep, stamp.sizeG_A);
        }
    }

    void updateStoredState(const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A) {
        for (size_t p = 0; p < numPorts; p++) {
            for (size_t c = 0; c < numPorts; c++) {
                for (size_t rho = 0; rho < port[p].from[c].numPoles; rho++) {
                    port[p].from[c].x[rho] = port[p].from[c].x[rho] *
                                                 port[p].from[c].exp_alpha[rho] +
                                             port[p].from[c].lambda_p[rho] *
                                                 awave_p(c, solutionMatrix,
                                                         currentSolutionIndex,
                                                         sizeG_A) +
                                             port[p].from[c].mu_p[rho] *
                                                 awave_p(c, solutionMatrix,
                                                         currentSolutionIndex - 1,
                                                         sizeG_A);
                    if (!firstOrder) {
                        port[p].from[c].x[rho] += port[p].from[c].nu_p[rho] *
                                                  awave_p(c, solutionMatrix,
                                                          currentSolutionIndex - 2,
                                                          sizeG_A);
                    }
                }
            }
        }

        if (firstOrder && currentSolutionIndex >= 1) {
            setSecondOrder(timestep);
        }
    }

    void setConstants(T timestep) {
        for (size_t p = 0; p < numPorts; p++) {
            port[p].beta = 1.0 - port[p].from[p].lambda - port[p].from[p].remainder;
            port[p].beta = 1.0 / port[p].beta;

            port[p].R = 1.0 + port[p].from[p].lambda + port[p].from[p].remainder;

            port[p].R = z_ref * port[p].R * port[p].beta;

            for (size_t c = 0; c < numPorts; c++) {
                if (c == p) {
                    port[p].alpha[c] = 0;
                } else {
                    port[p].alpha[c] = port[p].from[c].lambda +
                                       port[p].from[c].remainder;
                    port[p].alpha[c] = port[p].alpha[c] * port[p].beta;
                }
            }
        }
    }

    void setFirstOrder(T timestep) {
        firstOrder = true;

        for (size_t p = 0; p < numPorts; p++) {
            for (size_t c = 0; c < numPorts; c++) {
                port[p].from[c].lambda = 0;
                port[p].from[c].mu = 0;
                port[p].from[c].nu = 0;
                for (size_t rho = 0; rho < port[p].from[c].numPoles; rho++) {
                    const auto & pole = port[p].from[c].pole[rho];
                    const auto & residue = port[p].from[c].residue[rho];
                    const auto a = pole * timestep;
                    const auto ea = std::exp(a);
                    port[p].from[c].lambda_p[rho] = -(residue / pole) *
                                                    (1.0 + (1.0 - ea) / (a));
                    port[p].from[c].lambda += port[p].from[c].lambda_p[rho];

                    port[p].from[c].mu_p[rho] = -(residue / pole) *
                                                ((ea - 1.0) / a - ea);
                    port[p].from[c].mu += port[p].from[c].mu_p[rho];

                    port[p].from[c].nu_p[rho] = 0;
                }
            }
        }

        setConstants(timestep);
    }

    void setSecondOrder(T timestep) {
        firstOrder = false;

        for (size_t p = 0; p < numPorts; p++) {
            for (size_t c = 0; c < numPorts; c++) {
                port[p].from[c].lambda = 0;
                port[p].from[c].mu = 0;
                port[p].from[c].nu = 0;
                for (size_t rho = 0; rho < port[p].from[c].numPoles; rho++) {
                    const auto & pole = port[p].from[c].pole[rho];
                    const auto & residue = port[p].from[c].residue[rho];
                    const auto a = pole * timestep;
                    const auto ea = std::exp(a);
                    port[p].from[c].lambda_p[rho] = -(residue / pole) *
                                                    ((1.0 - ea) / (a * a) +
                                                     (3.0 - ea) / (2.0 * a) + 1.0);
                    port[p].from[c].lambda += port[p].from[c].lambda_p[rho];

                    port[p].from[c].mu_p[rho] = -(residue / pole) *
                                                (-2.0 * (1.0 - ea) / (a * a) -
                                                 (2.0 / a) - ea);
                    port[p].from[c].mu += port[p].from[c].mu_p[rho];

                    port[p].from[c].nu_p[rho] = -(residue / pole) *
                                                ((1.0 - ea) / (a * a) +
                                                 (1.0 + ea) / (2.0 * a));
                    port[p].from[c].nu += port[p].from[c].nu_p[rho];
                }
            }
        }

        setConstants(timestep);
    }

    void setTimestep(T timestep) {
        for (size_t p = 0; p < numPorts; p++) {
            for (size_t c = 0; c < numPorts; c++) {
                for (size_t rho = 0; rho < port[p].from[c].numPoles; rho++) {
                    port[p].from[c].exp_alpha[rho] = std::exp(
                        port[p].from[c].pole[rho] * timestep);
                }
            }

            if (firstOrder) {
                setFirstOrder(timestep);
            } else {
                setSecondOrder(timestep);
            }
        }
    }
#ifdef WITH_MATLAB
    void
    performVectorFit(std::string filePath, size_t numPorts,
                     std::shared_ptr<matlab::engine::MATLABEngine> matlabEngine) {
        matlab::data::ArrayFactory factory;

        matlabEngine->eval(u"addpath('./Matlab');");
        std::vector<matlab::data::Array> args;
        args.push_back(factory.createCharArray(filePath));
        auto result = matlabEngine->feval(u"CPPVectFitAdaptor", 3, args);

        assert(result[0].getDimensions()[0] == numPorts);
        assert(result[0].getDimensions()[1] == numPorts);

        z_ref = static_cast<T>(result[1][0]);
        for (size_t a = 0; a < numPorts; a++) {
            port[a].alpha.resize(numPorts);
            port[a].from.resize(numPorts);
            for (size_t b = 0; b < numPorts; b++) {
                auto structArrayResult = static_cast<
                    matlab::data::TypedArray<matlab::data::Struct> >(result[0]);
                auto structResult = static_cast<matlab::data::Struct>(
                    structArrayResult[a][b]);

                port[a].from[b].numPoles = structResult["poles"]
                                               .getNumberOfElements();
                port[a].from[b].lambda_p.resize(port[a].from[b].numPoles);
                port[a].from[b].mu_p.resize(port[a].from[b].numPoles);
                port[a].from[b].nu_p.resize(port[a].from[b].numPoles);
                port[a].from[b].exp_alpha.resize(port[a].from[b].numPoles);
                port[a].from[b].x.resize(port[a].from[b].numPoles);
                port[a].from[b].pole.resize(port[a].from[b].numPoles);
                port[a].from[b].residue.resize(port[a].from[b].numPoles);

                for (size_t p = 0; p < port[a].from[b].numPoles; p++) {
                    port[a].from[b].pole[p] = static_cast<std::complex<T> >(
                        static_cast<matlab::data::TypedArray<std::complex<T> > >(
                            structResult["poles"])[p]);
                    port[a].from[b].residue[p] = static_cast<std::complex<T> >(
                        static_cast<matlab::data::TypedArray<std::complex<T> > >(
                            structResult["residues"])[p]);
                }
                port[a].from[b].remainder = static_cast<std::complex<T> >(
                    static_cast<matlab::data::TypedArray<std::complex<T> > >(
                        structResult["remainder"])[0]);
            }
        }
    }
#endif

    void readInPRR(std::string filePath, size_t numPorts) {
        std::ifstream file(filePath);


        std::string line;
        while (file.peek() == '#' || file.peek() == '!') {
            std::getline(file, line);
        }

        T rval;
        T cval;

        std::stringstream polesLine;
        std::stringstream residuesLine;

        file >> z_ref;

        while (!file.eof()) {
            for (size_t a = 0; a < numPorts; a++) {
                port[a].alpha.resize(numPorts);
                port[a].from.resize(numPorts);
                for (size_t c = 0; c < numPorts; c++) {
                    file >> rval >> cval;
                    port[a].from[c].remainder = std::complex(rval, cval);

                    if (file.fail()) {
                        break;
                    }

                    std::
                        getline(file,
                                line); // have to clear the end of the remainder line
                    std::getline(file, line);
                    polesLine = std::stringstream(line);
                    std::getline(file, line);
                    residuesLine = std::stringstream(line);

                    while (!polesLine.eof() && !residuesLine.eof()) {
                        polesLine >> rval >> cval;
                        port[a].from[c].pole.emplace_back(rval, cval);

                        residuesLine >> rval >> cval;
                        port[a].from[c].residue.emplace_back(rval, cval);
                    }

                    port[a].from[c].numPoles = port[a].from[c].pole.size();
                    port[a].from[c].lambda_p.resize(port[a].from[c].numPoles);
                    port[a].from[c].mu_p.resize(port[a].from[c].numPoles);
                    port[a].from[c].nu_p.resize(port[a].from[c].numPoles);
                    port[a].from[c].exp_alpha.resize(port[a].from[c].numPoles);
                    port[a].from[c].x.resize(port[a].from[c].numPoles);
                }
            }

            if (file.fail()) {
                break;
            }
        }
    }

    void addDCAnalysisStampTo(Stamp<T> & stamp, const Matrix<T> & solutionVector,
                              size_t numCurrents) const {
        size_t numPorts = port.size();
        std::vector<std::complex<T> > xSum(port.size() * port.size());
        for (size_t p = 0; p < port.size(); p++) {
            for (size_t c = 0; c < port.size(); c++) {
                for (size_t rho = 0; rho < port[p].from[c].numPoles; rho++) {
                    xSum[p * numPorts + c] += -(port[p].from[c].lambda_p[rho] +
                                                port[p].from[c].mu_p[rho]) /
                                              (port[p].from[c].exp_alpha[rho] - 1.0);
                }
                xSum[p * numPorts + c] += port[p].from[c].remainder;
            }
        }

        for (size_t p = 0; p < port.size(); p++) {
            size_t np = port[p].positive - 1;
            size_t nn = port[p].negative - 1;
            size_t curr = port[p].current - 1;

            std::complex<T> beta = 1.0 / (1.0 - xSum[p * numPorts + p]);

            stamp.G(stamp.sizeG_A + curr, stamp.sizeG_A + curr) += std::real(
                -z_ref * (1.0 + xSum[p * numPorts + p]) * beta);

            if (port[p].positive != 0) {
                stamp.G(stamp.sizeG_A + curr, np) += 1;
                stamp.G(np, stamp.sizeG_A + curr) += 1;
            }

            if (port[p].negative != 0) {
                stamp.G(stamp.sizeG_A + curr, nn) += -1;
                stamp.G(nn, stamp.sizeG_A + curr) += -1;
            }
            // controlled sources
            for (size_t c = 0; c < port.size(); c++) {
                if (c != p) {
                    if (port[c].positive != 0) {
                        stamp.G(stamp.sizeG_A + curr,
                                port[c].positive -
                                    1) += std::real(-xSum[p * numPorts + c] * beta);
                    }
                    if (port[c].negative != 0) {
                        stamp.G(stamp.sizeG_A + curr,
                                port[c].negative -
                                    1) += std::real(xSum[p * numPorts + c] * beta);
                    }
                    stamp.G(stamp.sizeG_A + curr,
                            stamp.sizeG_A + port[c].current -
                                1) += std::real(-z_ref * xSum[p * numPorts + c] *
                                                beta);
                }
            }
        }
    }

    static void
    addToElements(const std::string & line, CircuitElements<T> & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents) {
        // std::regex SParameterBlockInitialRegex( R"(^S(.*?)\s(\d+?)\s)" );
        std::regex SParameterBlockInitialRegex = generateRegex("SV", "n s", true,
                                                               false);
        std::smatch matches;
        std::regex_search(line, matches, SParameterBlockInitialRegex);

        SParameterBlockVF<T> block;

        block.designator = "SV";
        block.designator += matches.str(1);

        size_t numPorts = std::stoul(matches.str(2));
        block.numPorts = numPorts;
        block.port = std::vector<SParameterPortVF<T> >(numPorts);

        auto strIter = matches.suffix().first;
        std::regex SParameterBlockPortRegex(R"(^(\d+?)\s(\d+?)\s)");
        for (size_t p = 0; p < numPorts; p++) {
            std::regex_search(strIter, line.cend(), matches,
                              SParameterBlockPortRegex);
            block.port[p].positive = std::stoi(matches.str(1));
            block.port[p].negative = std::stoi(matches.str(2));

            numNodes = std::max(numNodes, std::stoull(matches.str(1)));
            numNodes = std::max(numNodes, std::stoull(matches.str(2)));

            block.port[p].current = ++numCurrents;

            strIter = matches.suffix().first;
        }
        std::regex endOfLineRegex(R"(^(.*)$)");
        std::regex_search(strIter, line.cend(), matches, endOfLineRegex);
        if (line[2] == 'F') {
#ifdef WITH_MATLAB
            block.performVectorFit(matches.str(1), numPorts, elements.matlabEngine);
#endif
        } else {
            block.readInPRR(matches.str(1), numPorts);
        }

        elements.dynamicElements.emplace_back(
            std::make_shared<SParameterBlockVF<T> >(block));
        for (size_t p = 0; p < numPorts; p++) {
            elements.nodeComponentMap.insert(
                {{block.port[p].positive, elements.dynamicElements.back()},
                 {block.port[p].negative, elements.dynamicElements.back()}});
        }
    }
};

#endif
