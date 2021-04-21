#ifndef _BJT_HPP_INC_
#define _BJT_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"

/// @brief A simple NPN BJT model
///
/// @tparam T the value type
template<typename T>
struct BJTN : public Component<T> {
public:
    size_t c = 0;
    size_t b = 0;
    size_t e = 0;

    T alpha_f = 0.99;
    T alpha_r = 0.02;

    const T I_es = 2e-14;
    const T V_Te = 26e-3;
    const T I_cs = 99e-14;
    const T V_Tc = 26e-3;

    T V_bc_crit = V_Tc * std::log(V_Tc / (I_cs * std::sqrt(2)));
    T V_be_crit = V_Te * std::log(V_Te / (I_es * std::sqrt(2)));


    void
    addNonLinearStampTo(Stamp<T> & stamp, const Matrix<T> & solutionMatrix,
                        const size_t currentSolutionIndex, T timestep = 0) const {
        const size_t bp = b - 1;
        const size_t cp = c - 1;
        const size_t ep = e - 1;

        T v_be = 0;
        T v_bc = 0;

        if (b > 0) {
            v_be = solutionMatrix(bp, currentSolutionIndex);
            v_bc = solutionMatrix(bp, currentSolutionIndex);
        }

        if (e > 0) {
            v_be -= solutionMatrix(ep, currentSolutionIndex);
        }

        if (c > 0) {
            v_bc -= solutionMatrix(cp, currentSolutionIndex);
        }

        v_be = std::min(V_be_crit, v_be);
        v_bc = std::min(V_bc_crit, v_bc);

        T i_e = -I_es * (std::exp(v_be / V_Te) - 1) +
                alpha_r * I_cs * (std::exp(v_bc / V_Tc) - 1);
        T i_c = alpha_f * I_es * (std::exp(v_be / V_Te) - 1) -
                I_cs * (std::exp(v_bc / V_Tc) - 1);
        // T i_b = - ( i_e + i_c ); This is unused

        T g_ee = (I_es / V_Te) * std::exp(v_be / V_Te);
        T g_ec = alpha_r * (I_cs / V_Tc) * std::exp(v_bc / V_Tc);
        T g_ce = alpha_f * (I_es / V_Te) * std::exp(v_be / V_Te);
        T g_cc = (I_cs / V_Tc) * std::exp(v_bc / V_Tc);

        T I_e = i_e + g_ee * v_be - g_ec * v_bc;
        T I_c = i_c - g_ce * v_be + g_cc * v_bc;

        if (e > 0) {
            stamp.G(ep, ep) += g_ee;
            stamp.s(ep, 0) += -I_e;

            if (c > 0) {
                stamp.G(ep, cp) += -g_ec;
            }

            if (b > 0) {
                stamp.G(ep, bp) += (g_ec - g_ee);
            }
        }

        if (c > 0) {
            stamp.G(cp, cp) += g_cc;
            stamp.s(cp, 0) += -I_c;

            if (e > 0) {
                stamp.G(cp, ep) += -g_ce;
            }

            if (b > 0) {
                stamp.G(cp, bp) += (g_ce - g_cc);
            }
        }

        if (b > 0) {
            stamp.G(bp, bp) += g_cc + g_ee - g_ce - g_ec;
            stamp.s(bp, 0) += I_e + I_c;

            if (e > 0) {
                stamp.G(bp, ep) += g_ce - g_ee;
            }

            if (c > 0) {
                stamp.G(bp, cp) += (g_ec - g_cc);
            }
        }
    }

    void updateStoredState(const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A) {
    }

    void addDCAnalysisStampTo(Stamp<T> & stamp, const Matrix<T> & solutionVector,
                              size_t numCurrents) const {
        addNonLinearStampTo(stamp, solutionVector, 0, 0);
    }

    static void
    addToElements(const std::string & line, CircuitElements<T> & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents) {
        std::regex BJTRegex = generateRegex("QN", "n n n");
        BJTN<T> bjt;
        std::smatch matches;

        std::regex_match(line, matches, BJTRegex);

        bjt.designator = "QN";
        bjt.designator += matches.str(1);

        bjt.c = std::stoi(matches.str(2));
        bjt.b = std::stoi(matches.str(3));
        bjt.e = std::stoi(matches.str(4));

        numNodes = std::max(numNodes, std::stoull(matches.str(2)));
        numNodes = std::max(numNodes, std::stoull(matches.str(3)));
        numNodes = std::max(numNodes, std::stoull(matches.str(4)));


        elements.nonLinearElements.emplace_back(std::make_shared<BJTN<T> >(bjt));
        elements.nodeComponentMap.insert(
            {{bjt.b, elements.nonLinearElements.back()},
             {bjt.c, elements.nonLinearElements.back()},
             {bjt.e, elements.nonLinearElements.back()}});
    }
};

/// @brief A simple PNP BJT model
///
/// @tparam T the value type
template<typename T>
struct BJTP : public Component<T> {
public:
    size_t c = 0;
    size_t b = 0;
    size_t e = 0;

    T alpha_f = 0.99;
    T alpha_r = 0.02;

    const T I_es = 2e-14;
    const T V_Te = 26e-3;
    const T I_cs = 99e-14;
    const T V_Tc = 26e-3;

    T V_bc_crit = V_Tc * std::log(V_Tc / (I_cs * std::sqrt(2)));
    T V_be_crit = V_Te * std::log(V_Te / (I_es * std::sqrt(2)));


    void
    addNonLinearStampTo(Stamp<T> & stamp, const Matrix<T> & solutionMatrix,
                        const size_t currentSolutionIndex, T timestep = 0) const {
        const size_t bp = b - 1;
        const size_t cp = c - 1;
        const size_t ep = e - 1;

        T v_be = 0;
        T v_bc = 0;

        if (b > 0) {
            v_be = solutionMatrix(bp, currentSolutionIndex);
            v_bc = solutionMatrix(bp, currentSolutionIndex);
        }

        if (e > 0) {
            v_be -= solutionMatrix(ep, currentSolutionIndex);
        }

        if (c > 0) {
            v_bc -= solutionMatrix(cp, currentSolutionIndex);
        }

        v_be = std::max(-V_be_crit, v_be);
        v_bc = std::max(-V_bc_crit, v_bc);

        T i_F = I_cs * (std::exp(-v_bc / V_Tc) - 1);
        T i_R = I_es * (std::exp(-v_be / V_Te) - 1);
        T di_F = -(I_cs / V_Tc) * std::exp(-v_bc / V_Tc);
        T di_R = -(I_es / V_Te) * std::exp(-v_be / V_Te);

        T i_e = i_R - alpha_f * i_F;
        T i_b = (alpha_f - 1) * i_F + (alpha_r - 1) * i_R;
        T i_c = i_F - alpha_r * i_R;

        T g_ee = di_F;
        T g_ec = -alpha_r * di_R;
        T g_ce = -alpha_f * di_F;
        T g_cc = di_R;
        T g_be = (alpha_r - 1) * di_R;
        T g_bc = (alpha_f - 1) * di_F;

        T I_e = i_e - g_ee * v_be - g_ec * v_bc;
        T I_c = i_c - g_ce * v_be - g_cc * v_bc;
        T I_b = i_b - g_be * v_be - g_bc * v_bc;

        if (e > 0) {
            stamp.G(ep, ep) += -g_ee;
            stamp.s(ep, 0) += -I_e;

            if (c > 0) {
                stamp.G(ep, cp) += -g_ec;
            }

            if (b > 0) {
                stamp.G(ep, bp) += g_ec + g_ee;
            }
        }

        if (c > 0) {
            stamp.G(cp, cp) += -g_cc;
            stamp.s(cp, 0) += -I_c;

            if (e > 0) {
                stamp.G(cp, ep) += -g_ce;
            }

            if (b > 0) {
                stamp.G(cp, bp) += g_ce + g_cc;
            }
        }

        if (b > 0) {
            stamp.G(bp, bp) += g_be + g_bc;
            stamp.s(bp, 0) += -I_b;

            if (e > 0) {
                stamp.G(bp, ep) += -g_be;
            }

            if (c > 0) {
                stamp.G(bp, cp) += -g_bc;
            }
        }
    }

    void updateStoredState(const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A) {
    }

    void addDCAnalysisStampTo(Stamp<T> & stamp, const Matrix<T> & solutionVector,
                              size_t numCurrents) const {
        addNonLinearStampTo(stamp, solutionVector, 0, 0);
    }

    static void
    addToElements(const std::string & line, CircuitElements<T> & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents) {
        std::regex BJTRegex = generateRegex("QP", "n n n");
        BJTP<T> bjt;
        std::smatch matches;

        std::regex_match(line, matches, BJTRegex);

        bjt.designator = "QP";
        bjt.designator += matches.str(1);

        bjt.c = std::stoi(matches.str(2));
        bjt.b = std::stoi(matches.str(3));
        bjt.e = std::stoi(matches.str(4));

        numNodes = std::max(numNodes, std::stoull(matches.str(2)));
        numNodes = std::max(numNodes, std::stoull(matches.str(3)));
        numNodes = std::max(numNodes, std::stoull(matches.str(4)));


        elements.nonLinearElements.emplace_back(std::make_shared<BJTP<T> >(bjt));
        elements.nodeComponentMap.insert(
            {{bjt.b, elements.nonLinearElements.back()},
             {bjt.c, elements.nonLinearElements.back()},
             {bjt.e, elements.nonLinearElements.back()}});
    }
};
#endif
