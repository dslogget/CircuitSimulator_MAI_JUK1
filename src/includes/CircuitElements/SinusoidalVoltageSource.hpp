#ifndef _SINUSOIDALVOLTAGESOURCE_INC_
#define _SINUSOIDALVOLTAGESOURCE_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"
#include <math.h>

/// @brief A sinusoidal voltage source model
///
/// @tparam T the value type
template<typename T>
struct SinusoidalVoltageSource : public Component<T> {
    size_t n1 = 0;
    size_t n2 = 0;
    size_t currentIndex = 0;

    T V = 1;
    T phase = 0;
    T frequency = 1;
    T offset = 0;
    bool degrees = true;

    void addDynamicStampTo(Stamp<T> & stamp, const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep) const {
        size_t n1p = n1 - 1;
        size_t n2p = n2 - 1;
        size_t currentIndexp = currentIndex - 1;

        if (n1) {
            stamp.G(n1p, stamp.sizeG_A + currentIndexp) += 1;
            stamp.G(stamp.sizeG_A + currentIndexp, n1p) += 1;
        }

        if (n2) {
            stamp.G(n2p, stamp.sizeG_A + currentIndexp) += -1;
            stamp.G(stamp.sizeG_A + currentIndexp, n2p) += -1;
        }

        if (degrees) {
            stamp.s(stamp.sizeG_A + currentIndexp,
                    0) += offset + V * std::sin(2 * std::numbers::pi * frequency *
                                                    currentSolutionIndex * timestep +
                                                std::numbers::pi * phase / 180);
        } else {
            stamp.s(stamp.sizeG_A + currentIndexp,
                    0) += offset + V * std::sin(2 * std::numbers::pi * frequency *
                                                    currentSolutionIndex * timestep +
                                                phase);
        }
    }

    void updateStoredState(const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A) {
    }

    void addDCAnalysisStampTo(Stamp<T> & stamp, const Matrix<T> & solutionVector,
                              size_t numCurrents) const {
        addDynamicStampTo(stamp, solutionVector, 0, 0);
    }

    static void
    addToElements(const std::string & line, CircuitElements<T> & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents) {
        // n1 n2 amp freq off phase
        std::regex sinusoidalVoltageSourceRegex = generateRegex("VS",
                                                                "n n ? w w w w");

        SinusoidalVoltageSource<T> voltageSource;
        std::smatch matches;

        std::regex_match(line, matches, sinusoidalVoltageSourceRegex);

        voltageSource.designator = "VS";
        voltageSource.designator += matches.str(1);

        voltageSource.n1 = std::stoi(matches.str(2));
        voltageSource.n2 = std::stoi(matches.str(3));

        numNodes = std::max(numNodes, std::stoull(matches.str(2)));
        numNodes = std::max(numNodes, std::stoull(matches.str(3)));

        if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
            voltageSource.V = std::stod(matches.str(4));
        } else {
            static_assert("Unsupported Type");
        }

        if (matches.size() > 4 && matches.str(5) != "") {
            voltageSource.frequency = std::stod(matches.str(5));
        }

        if (matches.size() > 5 && matches.str(6) != "") {
            voltageSource.offset = std::stod(matches.str(6));
        }

        if (matches.size() > 6 && matches.str(7) != "") {
            voltageSource.phase = std::stod(matches.str(7));
        }

        voltageSource.currentIndex = ++numCurrents;

        elements.dynamicElements.emplace_back(
            std::make_shared<SinusoidalVoltageSource<T> >(voltageSource));
        elements.nodeComponentMap.insert(
            {{voltageSource.n1, elements.dynamicElements.back()},
             {voltageSource.n2, elements.dynamicElements.back()}});
    }
};

#endif
