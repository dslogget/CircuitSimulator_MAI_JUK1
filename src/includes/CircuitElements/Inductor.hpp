#ifndef _INDUCTOR_HPP_INC_
#define _INDUCTOR_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"


/// @brief An ideal inductor model
///
/// @tparam T the value type
template<typename T>
struct Inductor : public Component<T> {
public:
    T value = 0;

    size_t n1 = 0;
    size_t n2 = 0;
    T lastCurrent = 0;

    size_t dcCurrentIndex = 0;

    bool trapezoidalRule = true;
    void addDynamicStampTo(Stamp<T> & stamp, const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep) const {
        size_t n1p = n1 - 1;
        size_t n2p = n2 - 1;

        T u0 = 0;
        if (n1) {
            u0 = solutionMatrix(n1p, currentSolutionIndex - 1);
        }

        if (n2) {
            u0 -= solutionMatrix(n2p, currentSolutionIndex - 1);
        }


        T G_eq = 0;
        T I_eq = 0;

        if (trapezoidalRule) {
            G_eq = timestep / (2 * value);
            I_eq = lastCurrent + G_eq * u0;
        } else {
            G_eq = timestep / timestep;
            I_eq = lastCurrent;
        }

        if (n1) {
            stamp.G(n1p, n1p) += G_eq;
            stamp.s(n1p, 0) += -I_eq;
        }

        if (n2) {
            stamp.G(n2p, n2p) += G_eq;
            stamp.s(n2p, 0) += I_eq;
        }

        if (n1 && n2) {
            stamp.G(n1p, n2p) += -G_eq;
            stamp.G(n2p, n1p) += -G_eq;
        }
    }

    void updateStoredState(const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A) {
        size_t n1p = n1 - 1;
        size_t n2p = n2 - 1;
        T u0 = 0;
        T u1 = 0;
        if (n1) {
            u0 = solutionMatrix(n1p, currentSolutionIndex - 1);
            u1 = solutionMatrix(n1p, currentSolutionIndex);
        }

        if (n2) {
            u0 -= solutionMatrix(n2p, currentSolutionIndex - 1);
            u1 -= solutionMatrix(n2p, currentSolutionIndex);
        }

        if (trapezoidalRule) {
            T G_eq = timestep / (2 * value);
            lastCurrent = G_eq * u1 + (lastCurrent + G_eq * u0);
        } else {
            T G_eq = timestep / value;
            lastCurrent = G_eq * u1 + lastCurrent;
        }
    }

    void updateDCStoredState(const Matrix<T> & solutionVector, size_t sizeG_A,
                             size_t numCurrents) {
        lastCurrent = solutionVector(sizeG_A + numCurrents + dcCurrentIndex - 1, 0);
    }

    void addDCAnalysisStampTo(Stamp<T> & stamp, const Matrix<T> & solutionVector,
                              size_t numCurrents) const {
        // Short circuit
        size_t n1p = n1 - 1;
        size_t n2p = n2 - 1;
        size_t dcCurrentIndexp = dcCurrentIndex - 1;

        if (n1 > 0) {
            stamp.G(n1p, stamp.sizeG_A + numCurrents + dcCurrentIndexp) += 1;
            stamp.G(stamp.sizeG_A + numCurrents + dcCurrentIndexp, n1p) += 1;
        }

        if (n2 > 0) {
            stamp.G(n2p, stamp.sizeG_A + numCurrents + dcCurrentIndexp) += -1;
            stamp.G(stamp.sizeG_A + numCurrents + dcCurrentIndexp, n2p) += -1;
        }
    }

    static void
    addToElements(const std::string & line, CircuitElements<T> & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents) {
        // std::regex inductorRegex( R"(^L(.*?)\s(\d+?)\s(\d+?)\s(.+?)\s?$)" );
        std::regex inductorRegex = generateRegex("L", "n n w");
        Inductor<T> inductor;
        std::smatch matches;

        std::regex_match(line, matches, inductorRegex);

        inductor.designator = "L";
        inductor.designator += matches.str(1);

        inductor.n1 = std::stoi(matches.str(2));
        inductor.n2 = std::stoi(matches.str(3));
        inductor.trapezoidalRule = true;
        inductor.dcCurrentIndex = ++numDCCurrents;

        numNodes = std::max(numNodes, std::stoull(matches.str(2)));
        numNodes = std::max(numNodes, std::stoull(matches.str(3)));

        if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
            inductor.value = std::stod(matches.str(4));
        } else {
            static_assert("Unsupported Type");
        }

        elements.dynamicElements.emplace_back(
            std::make_shared<Inductor<T> >(inductor));

        elements.nodeComponentMap.insert(
            {{inductor.n1, elements.dynamicElements.back()},
             {inductor.n2, elements.dynamicElements.back()}});
    }
};

#endif
