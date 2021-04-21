#ifndef _CIRCUITELEMENTS_HPP_INC_
#define _CIRCUITELEMENTS_HPP_INC_
#include "CircuitElements/Resistor.hpp"
#include "CircuitElements/Capacitor.hpp"
#include "CircuitElements/NLCapacitor.hpp"
#include "CircuitElements/Inductor.hpp"
#include "CircuitElements/VoltageSource.hpp"
#include "CircuitElements/SinusoidalVoltageSource.hpp"
#include "CircuitElements/TimeSeriesVoltageSource.hpp"
#include "CircuitElements/CurrentSource.hpp"
#include "CircuitElements/BJT.hpp"
#include "CircuitElements/Diode.hpp"
#include "CircuitElements/SParameterBlock.hpp"
#include "CircuitElements/SParameterBlockVF.hpp"
#include "CircuitElements/NLNMOS.hpp"
#include "CircuitElements/NLCurrentSource.hpp"

#ifdef WITH_MATLAB
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#endif

/// @brief An enum to track how far we want to go in a solution.
enum class SolutionStage {
    StaticSolution = 0,
    DynamicSolution = 1,
    NonLinearSolution = 2
};

/// @brief An enum for component types. Current unused
enum class ComponentType {
    Resistor,
    Capacitor,
    Inductor,
    VoltageSource,
    CurrentSource,
    BJT,
    Diode
};


/// @brief a glorified container for the different types of components.
///
/// @tparam T
template<typename T>
struct CircuitElements {
    /// @brief Preallocated stamp. Used for caching between loop iterations. Static
    ///        stamps will only be generated once as a result.
    Stamp<T> staticStamp;
    /// @brief Preallocated stamp. Used for caching between loop iterations. Dynamic
    ///        stamps need to be updated on every timestep.
    Stamp<T> dynamicStamp;
    /// @brief Preallocated stamp. Used for caching between loop iterations.
    /// Non-Linear
    ///        stamps must be updated on every newton-raphson iteration.
    Stamp<T> nonLinearStamp;
    /// @brief Preallocated stamp. Used for caching between loop iterations.
    /// DC
    ///        stamps must be updated on every newton-raphson iteration.
    Stamp<T> dcStamp;

#ifdef WITH_MATLAB
    std::shared_ptr<matlab::engine::MATLABEngine> matlabEngine;
#endif

    /// @brief A container to store the static components
    std::vector<std::shared_ptr<Component<T> > > staticElements;
    /// @brief A container to store the dynamic components
    std::vector<std::shared_ptr<Component<T> > > dynamicElements;
    /// @brief A container to store the Non-Linear components
    std::vector<std::shared_ptr<Component<T> > > nonLinearElements;

    /// @brief A variable used to track if the cached stamp is current.
    bool staticStampIsFresh = false;
    /// @brief A variable used to track if the cached stamp is current.
    bool dynamicStampIsFresh = false;
    /// @brief A variable used to track if the cached stamp is current.
    bool nonLinearStampIsFresh = false;

    /// @brief A map to pair nodes with the components connected to them.
    std::multimap<size_t, std::shared_ptr<Component<T> > > nodeComponentMap;

    /// @brief Initialisation.
    ///
    /// @param numNodes The size of the stamps voltage dependants
    /// @param numCurrents The size of the stamps voltage dependants
    /// @param numDCCurrents The size of the stamps voltage dependants
    CircuitElements(size_t numNodes = 0, size_t numCurrents = 0,
                    size_t numDCCurrents = 0)
        : staticStamp(numNodes, numCurrents), dynamicStamp(numNodes, numCurrents),
          nonLinearStamp(numNodes, numCurrents),
          dcStamp(numNodes, numCurrents + numDCCurrents) {
    }

    /// @brief Updates the size of all stamps
    ///
    /// @param numNodes The size of the stamps voltage dependants
    /// @param numCurrents The size of the stamps voltage dependants
    /// @param numDCCurrents The size of the stamps voltage dependants
    void
    setNewStampSize(size_t numNodes, size_t numCurrents, size_t numDCCurrents = 0) {
        staticStamp = Stamp<T>(numNodes, numCurrents);
        dynamicStamp = Stamp<T>(numNodes, numCurrents);
        nonLinearStamp = Stamp<T>(numNodes, numCurrents);
        dcStamp = Stamp<T>(numNodes, numCurrents + numDCCurrents);
        staticStampIsFresh = false;
        dynamicStampIsFresh = false;
        nonLinearStampIsFresh = false;
    }

    /// @brief Forces a clear of the static stamp, and generates a new one.
    ///
    /// @return A reference to the cached stamp.
    Stamp<T> & generateStaticStamp() {
        staticStamp.clear();

        for (const auto & component : staticElements) {
            staticStamp.addStaticStamp(component);
        }

        for (const auto & component : dynamicElements) {
            staticStamp.addStaticStamp(component);
        }

        for (const auto & component : nonLinearElements) {
            staticStamp.addStaticStamp(component);
        }

        staticStampIsFresh = true;
        return staticStamp;
    }

    /// @brief Obtains the static stamp, then adds dynamic components to it.
    ///
    /// @param solutionMatrix The solution matrix to use for the dynamic and
    /// non-linear stamp.
    /// @param currentSolutionIndex The current index we are at.
    /// @param timestep The time step being used.
    ///
    /// @return A reference to the cached stamp.
    Stamp<T> & generateDynamicStamp(const Matrix<T> & solutionMatrix,
                                    const size_t currentSolutionIndex, T timestep) {
        if (!staticStampIsFresh) {
            generateStaticStamp();
        }
        dynamicStamp = staticStamp;

        for (const auto & component : dynamicElements) {
            dynamicStamp.addDynamicStamp(component, solutionMatrix,
                                         currentSolutionIndex, timestep);
        }

        for (const auto & component : nonLinearElements) {
            dynamicStamp.addDynamicStamp(component, solutionMatrix,
                                         currentSolutionIndex, timestep);
        }

        dynamicStampIsFresh = true;
        return dynamicStamp;
    }

    /// @brief Obtains the dynamic stamp, then adds dynamic components to it.
    ///
    /// @param solutionMatrix The solution matrix to use for the dynamic and
    /// non-linear stamp.
    /// @param currentSolutionIndex The current index we are at.
    /// @param timestep The time step being used.
    ///
    /// @return A reference to the cached stamp.
    Stamp<T> &
    generateNonLinearStamp(const Matrix<T> & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep) {
        if (!dynamicStampIsFresh) {
            generateDynamicStamp(solutionMatrix, currentSolutionIndex, timestep);
        }
        nonLinearStamp = dynamicStamp;

        for (const auto & component : nonLinearElements) {
            nonLinearStamp.addNonLinearStamp(component, solutionMatrix,
                                             currentSolutionIndex, timestep);
        }

        nonLinearStampIsFresh = true;
        return nonLinearStamp;
    }

    /// @brief Generates the complete stamp up to a certain point.
    ///
    /// @param solutionMatrix The solution matrix to use for the dynamic and
    /// non-linear stamp.
    /// @param currentSolutionIndex The current index we are at.
    /// @param timestep The time step being used.
    ///
    /// @return The complete stamp.
    Stamp<T> &
    generateCompleteStamp(SolutionStage stage, const Matrix<T> & solutionMatrix,
                          const size_t currentSolutionIndex, T timestep) {
        switch (stage) {
            case SolutionStage::StaticSolution:
                if (!staticStampIsFresh) {
                    generateStaticStamp()
                }
                return staticStamp;
                break;

            case SolutionStage::DynamicSolution:
                if (!dynamicStampIsFresh) {
                    generateDynamicStamp(solutionMatrix, currentSolutionIndex,
                                         timestep)
                }
                return dynamicStamp;
                break;

            case SolutionStage::NonLinearSolution:
                if (!nonLinearStampIsFresh) {
                    generateNonLinearStamp(solutionMatrix, currentSolutionIndex,
                                           timestep)
                }
                return nonLinearStamp;
                break;
        }
    }

    /// @brief Generates the complete DC stamp
    ///
    /// @param solutionVector The solution vector to use for the dynamic and
    /// non-linear stamp.
    ///
    /// @param currentSolutionIndex The current index we are at.
    /// @param timestep The time step being used.
    ///
    /// @return The complete stamp.
    Stamp<T> &
    generateDCStamp(const Matrix<T> & solutionVector, size_t numCurrents) {
        dcStamp.clear();

        for (const auto & component : staticElements) {
            dcStamp.addDCAnalysisStamp(component, solutionVector, numCurrents);
        }

        for (const auto & component : dynamicElements) {
            dcStamp.addDCAnalysisStamp(component, solutionVector, numCurrents);
        }

        for (const auto & component : nonLinearElements) {
            dcStamp.addDCAnalysisStamp(component, solutionVector, numCurrents);
        }

        return dcStamp;
    }

    /// @brief Updates the components at the end of each time step. Applies to
    ///        dynamic and non-linear components.
    ///
    /// @param solutionMatrix The solution matrix to use for the dynamic and
    /// non-linear stamp.
    /// @param currentSolutionIndex The current index we are at.
    /// @param timestep The time step being used.
    void updateTimeStep(const Matrix<T> & solutionMatrix,
                        const size_t currentSolutionIndex, T timestep) {
        dynamicStampIsFresh = false;
        nonLinearStampIsFresh = false;
        for (const auto & component : dynamicElements) {
            component->updateStoredState(solutionMatrix, currentSolutionIndex,
                                         timestep, staticStamp.sizeG_A);
        }
        for (const auto & component : nonLinearElements) {
            component->updateStoredState(solutionMatrix, currentSolutionIndex,
                                         timestep, staticStamp.sizeG_A);
        }
    }

    /// @brief Updates the components based on their DC value. Applies to
    ///        dynamic and non-linear components.
    ///
    /// @param solutionVector The solution vector to use for DC value
    /// @param numCurrents The number of currents used by the transient simulation.
    void updateDCStoredState(const Matrix<T> & solutionVector, size_t numCurrents) {
        for (const auto & component : staticElements) {
            component->updateDCStoredState(solutionVector, dcStamp.sizeG_A,
                                           numCurrents);
        }
        for (const auto & component : dynamicElements) {
            component->updateDCStoredState(solutionVector, dcStamp.sizeG_A,
                                           numCurrents);
        }
        for (const auto & component : nonLinearElements) {
            component->updateDCStoredState(solutionVector, dcStamp.sizeG_A,
                                           numCurrents);
        }
    }
};

#endif
