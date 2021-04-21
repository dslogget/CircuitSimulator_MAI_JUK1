#ifndef _RESISTOR_HPP_INC_
#define _RESISTOR_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"

/// @brief An ideal resistor model
///
/// @tparam T the value type
template<typename T>
struct Resistor : public Component<T> {
public:
    T value = 0;

    size_t n1 = 0;
    size_t n2 = 0;
    size_t currentIndex = 0;

    bool group1 = true;

    void addStaticStampTo(Stamp<T> & stamp) const {
        // the p means prime ( ' ) and is used for the index - 1
        size_t n1p = n1 - 1;
        size_t n2p = n2 - 1;
        size_t currentIndexp = currentIndex - 1;

        if (group1) {
            T conductance = 1 / value;

            if (n1) {
                stamp.G(n1p, n1p) += conductance;
            }
            if (n2) {
                stamp.G(n2p, n2p) += conductance;
            }
            if (n1 && n2) {
                stamp.G(n1p, n2p) += -conductance;
                stamp.G(n2p, n1p) += -conductance;
            }
        } else {
            if (n1) {
                stamp.G(n1p, stamp.sizeG_A + currentIndexp) += 1;
                stamp.G(stamp.sizeG_A + currentIndexp, n1p) += 1;
            }
            if (n2) {
                stamp.G(n2p, stamp.sizeG_A + currentIndexp) += -1;
                stamp.G(stamp.sizeG_A + currentIndexp, n2p) += -1;
            }

            stamp.G(stamp.sizeG_A + currentIndexp,
                    stamp.sizeG_A + currentIndexp) += -value;
        }
    }

    void addDCAnalysisStampTo(Stamp<T> & stamp, const Matrix<T> & solutionVector,
                              size_t numCurrents) const {
        addStaticStampTo(stamp);
    }

    static void
    addToElements(const std::string & line, CircuitElements<T> & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents) {
        // std::regex resistorRegex(
        // R"(^R(.*?)\s(\d+?)\s(\d+?)\s(.+?)(?:\s(.))?\s?$)"
        // );
        std::regex resistorRegex = generateRegex("R", "n n w ? c");
        Resistor<T> resistor;
        std::smatch matches;

        std::regex_match(line, matches, resistorRegex);

        resistor.designator = "R";
        resistor.designator += matches.str(1);

        resistor.n1 = std::stoi(matches.str(2));
        resistor.n2 = std::stoi(matches.str(3));

        numNodes = std::max(numNodes, std::stoull(matches.str(2)));
        numNodes = std::max(numNodes, std::stoull(matches.str(3)));

        if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
            resistor.value = std::stod(matches.str(4));
        } else {
            static_assert("Unsupported Type");
        }

        if (matches.size() < 6 && matches.str(5) != "") {
            resistor.group1 = false;
            resistor.currentIndex = ++numCurrents;
        } else {
            resistor.group1 = true;
        }

        elements.staticElements.emplace_back(
            std::make_shared<Resistor<T> >(resistor));

        elements.nodeComponentMap.insert(
            {{resistor.n1, elements.staticElements.back()},
             {resistor.n2, elements.staticElements.back()}});
    }
};

#endif
