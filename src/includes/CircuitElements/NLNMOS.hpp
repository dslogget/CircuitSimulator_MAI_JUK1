#ifndef _NLNMOS_HPP_INC_
#define _NLNMOS_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"
#include <math.h>

/// @brief a non-linear FET model
/// @tparam T the value type
template<typename T>
struct NLNMOS : public Component<T> {
public:
    size_t d = 0;
    size_t g = 0;
    size_t s = 0;

    // constant params of the model
    const T C_GSp = 0.01;
    const T C_GSo = 0.5;
    const T P_S10 = 0;
    const T P_S11 = 0.5;
    const T C_GDp = 0.5;
    const T C_GDo = 1;
    const T P_D10 = -1;
    const T P_D11 = 0.4;

    const T beta_DS = 1.3;
    const T alpha_DS = 0.42;

    T u_gd_last = 0;
    T u_gs_last = 0;


    T i_gd_last = 0;
    T i_gs_last = 0;

    T C_GD_last = C_GDp + C_GDo * (1.0 + std::tanh(P_D10 + P_D11 * u_gd_last));
    T C_GS_last = C_GSp + C_GSo * (1.0 + std::tanh(P_S10 + P_S11 * u_gs_last));

    void
    addNonLinearStampTo(Stamp<T> & stamp, const Matrix<T> & solutionMatrix,
                        const size_t currentSolutionIndex, T timestep = 0) const {
        const size_t gp = g - 1;
        const size_t dp = d - 1;
        const size_t sp = s - 1;

        T u_gs = 0;
        T u_gd = 0;

        if (g > 0) {
            u_gs = solutionMatrix(gp, currentSolutionIndex);
            u_gd = solutionMatrix(gp, currentSolutionIndex);
        }

        if (s > 0) {
            u_gs -= solutionMatrix(sp, currentSolutionIndex);
        }

        if (d > 0) {
            u_gd -= solutionMatrix(dp, currentSolutionIndex);
        }

        T C_GD = C_GDp + C_GDo * (1.0 + std::tanh(P_D10 + P_D11 * u_gd));
        T C_GS = C_GSp + C_GSo * (1.0 + std::tanh(P_S10 + P_S11 * u_gs));

        T dC_GD = C_GDo * P_D11 / std::pow(std::cosh(P_D10 + P_D11 * u_gd), 2);
        T dC_GS = C_GSo * P_S11 / std::pow(std::cosh(P_S10 + P_S11 * u_gs), 2);

        T i_ds = beta_DS * std::tanh(alpha_DS * (u_gs - u_gd));
        T di_ds_d = -beta_DS * alpha_DS /
                    std::pow(std::cosh(alpha_DS * (u_gs - u_gd)), 2);
        T di_ds_s = beta_DS * alpha_DS /
                    std::pow(std::cosh(alpha_DS * (u_gs - u_gd)), 2);

        T i_gd = C_GD *
                 (2.0 * (u_gd - u_gd_last) / timestep - i_gd_last / C_GD_last);
        T i_gs = C_GS *
                 (2.0 * (u_gs - u_gs_last) / timestep - i_gs_last / C_GS_last);

        T i_d = -i_gd + i_ds;
        T i_s = -i_gs - i_ds;
        T i_g = i_gs + i_gd;

        T di_gd = dC_GD *
                      (2.0 * (u_gd - u_gd_last) / timestep - i_gd_last / C_GD_last) +
                  2.0 * C_GD / timestep;
        T di_gs = dC_GS *
                      (2.0 * (u_gs - u_gs_last) / timestep - i_gs_last / C_GS_last) +
                  2.0 * C_GS / timestep;

        T g_dd = -di_gd + di_ds_d;
        T g_sd = -di_ds_d;
        T g_gd = di_gd;

        T g_ds = di_ds_d;
        T g_ss = -di_gs - di_ds_s;
        T g_gs = di_gs;

        T I_d = i_d - g_dd * u_gd - g_ds * u_gs;
        T I_s = i_s - g_sd * u_gd - g_ss * u_gs;
        T I_g = i_g - g_gd * u_gd - g_gs * u_gs;

        if (d > 0) {
            stamp.G(dp, dp) += -g_dd;
            stamp.s(dp, 0) += -I_d;

            if (s > 0) {
                stamp.G(dp, sp) += -g_ds;
            }

            if (g > 0) {
                stamp.G(dp, gp) += g_dd + g_ds;
            }
        }

        if (s > 0) {
            stamp.G(sp, sp) += -g_ss;
            stamp.s(sp, 0) += -I_s;

            if (d > 0) {
                stamp.G(sp, dp) += -g_sd;
            }

            if (g > 0) {
                stamp.G(sp, gp) += g_sd + g_ss;
            }
        }

        if (g > 0) {
            stamp.G(gp, gp) += g_gd + g_gs;
            stamp.s(gp, 0) += -I_g;

            if (d > 0) {
                stamp.G(gp, dp) += -g_gd;
            }

            if (s > 0) {
                stamp.G(gp, sp) += -g_gs;
            }
        }
    }

    void updateStoredState(const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A) {
        const size_t gp = g - 1;
        const size_t dp = d - 1;
        const size_t sp = s - 1;

        T u_gs = 0;
        T u_gd = 0;

        if (g > 0) {
            u_gs = solutionMatrix(gp, currentSolutionIndex);
            u_gd = solutionMatrix(gp, currentSolutionIndex);
        }

        if (s > 0) {
            u_gs -= solutionMatrix(sp, currentSolutionIndex);
        }

        if (d > 0) {
            u_gd -= solutionMatrix(dp, currentSolutionIndex);
        }

        T C_GD = C_GDp + C_GDo * (1.0 + std::tanh(P_D10 + P_D11 * u_gd));
        T C_GS = C_GSp + C_GSo * (1.0 + std::tanh(P_S10 + P_S11 * u_gs));

        i_gd_last = C_GD *
                    (2.0 * (u_gd - u_gd_last) / timestep - i_gd_last / C_GD_last);
        i_gs_last = C_GS *
                    (2.0 * (u_gs - u_gs_last) / timestep - i_gs_last / C_GS_last);

        C_GD_last = C_GD;
        C_GS_last = C_GS;

        u_gd_last = u_gd;
        u_gs_last = u_gs;
    }

    static void
    addToElements(const std::string & line, CircuitElements<T> & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents) {
        std::regex BJTRegex = generateRegex("QMN", "n n n");
        NLNMOS<T> nmos;
        std::smatch matches;

        std::regex_match(line, matches, BJTRegex);

        nmos.designator = "QMN";
        nmos.designator += matches.str(1);

        nmos.d = std::stoi(matches.str(2));
        nmos.g = std::stoi(matches.str(3));
        nmos.s = std::stoi(matches.str(4));

        numNodes = std::max(numNodes, std::stoull(matches.str(2)));
        numNodes = std::max(numNodes, std::stoull(matches.str(3)));
        numNodes = std::max(numNodes, std::stoull(matches.str(4)));


        elements.nonLinearElements.emplace_back(std::make_shared<NLNMOS<T> >(nmos));
        elements.nodeComponentMap.insert(
            {{nmos.d, elements.nonLinearElements.back()},
             {nmos.g, elements.nonLinearElements.back()},
             {nmos.s, elements.nonLinearElements.back()}});
    }
};
#endif
