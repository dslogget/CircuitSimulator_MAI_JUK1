#ifndef _COMPONENT_HPP_INC_
#define _COMPONENT_HPP_INC_
#include "Maths/DynamicMatrix.hpp"
#include "CircuitElements/ElementsRegexBuilder.h"
#include <regex>
#include <algorithm>
#include <map>

template< typename T >
struct Component;

/// @brief A helper struct to store the preallocated stamps for MNA
///
/// @tparam T The value type
///
/// @details The matrix has the following structure:\n
/// \verbatim
///       | G_A | G_B |
///   G = ------|------
///       | G_C | G_D |\endverbatim
///
template< typename T >
struct Stamp {
   size_t sizeG_A;
   size_t sizeG_D;

   Matrix< T > G;
   Matrix< T > s;

   /// @brief Sets the initial size of the stamp pair.
   ///
   /// @param _sizeG_A Size of the voltage dependence portion of the stamps (Group I)
   /// @param _sizeG_D Size of the current dependence portion of the stamps (Group
   /// II)
   Stamp( size_t _sizeG_A, size_t _sizeG_D )
       : sizeG_A( _sizeG_A ), sizeG_D( _sizeG_D ),
         G( _sizeG_A + _sizeG_D, _sizeG_A + _sizeG_D, 0 ),
         s( _sizeG_A + _sizeG_D, 1, 0 ) {
   }

   /// @brief Clears the stamps to 0s
   void clear() {
      G.fill( 0 );
      s.fill( 0 );
   }

   /// @brief Combines two stamps together into the current stamp.
   ///        This is not done via the operator, as it has side-effects
   ///
   /// @param rhs The stamps to add to the current stamp
   void add( const Stamp< T > & rhs ) {
      G.add( rhs.G, G );
      s.add( rhs.s, s );
   }

   /// @brief A helper function to add a static component to the stamp.
   ///
   /// @param rhs Component to be added
   void addStaticStamp( const std::shared_ptr< Component< T > > & rhs ) {
      rhs->addStaticStampTo( *this );
   }

   /// @brief A helper function to add a dynamic component to the stamp.
   ///
   /// @param rhs Component to be added
   void addDynamicStamp( const std::shared_ptr< Component< T > > & rhs,
                         const Matrix< T > & solutionMatrix,
                         const size_t currentSolutionIndex, T timestep ) {
      rhs->addDynamicStampTo( *this, solutionMatrix, currentSolutionIndex,
                              timestep );
   }

   /// @brief A helper function to add a non-linear component to the stamp.
   ///
   /// @param rhs Component to be added
   void addNonLinearStamp( const std::shared_ptr< Component< T > > & rhs,
                           const Matrix< T > & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep = 0 ) {
      rhs->addNonLinearStampTo( *this, solutionMatrix, currentSolutionIndex,
                                timestep );
   }

   /// @brief A helper function to add a DC component to the stamp.
   ///
   /// @param rhs Component to be added
   void addDCAnalysisStamp( const std::shared_ptr< Component< T > > & rhs,
                            const Matrix< T > & solutionMatrix,
                            const size_t numCurrents ) {
      rhs->addDCAnalysisStampTo( *this, solutionMatrix, numCurrents );
   }

   /// @brief An alias for left dividing G by s.
   ///
   /// @return a solution vector the same dimension as s.
   Matrix< T > solve() {
      return G.leftDivide( s );
   }
};

template< typename T >
struct CircuitElements;
/// @brief A template base class to define the fundamental things a component
///        should define.
///
/// @tparam T The value type
template< typename T >
struct Component {
   /// @brief The designator as in the netlist for e.g
   std::string designator = "";

   /// @brief Adds this component's static stamp to the target stamp.
   ///
   /// @param destination The stamp to be added to.
   virtual void addStaticStampTo( Stamp< T > & destination ) const {
   }

   /// @brief Adds this component's dynamic stamp to the target stamp.
   ///
   /// @param destination The stamp to be added to.
   /// @param solutionMatrix A vector containing all past solutions to the circuit
   /// @param currentSolutionIndex The current timeStep index
   /// @param timestep The length of each time step
   virtual void
   addDynamicStampTo( Stamp< T > & destination,
                      const Matrix< T > & solutionMatrix,
                      const size_t currentSolutionIndex, T timestep ) const {
   }

   /// @brief adds this component's non-linear stamp to the target stamp.
   ///
   /// @param destination The stamp to be added to.
   /// @param solutionMatrix A vector containing all past solutions to the circuit
   /// @param currentSolutionIndex The current timeStep index
   /// @param timestep The length of each time step
   virtual void
   addNonLinearStampTo( Stamp< T > & destination,
                        const Matrix< T > & solutionMatrix,
                        const size_t currentSolutionIndex, T timestep = 0 ) const {
      throw std::exception( "not implemented" );
   }

   /// @brief Updates any stored state based on the current solution index
   ///
   /// @param solutionMatrix A vector containing all past solutions to the circuit
   /// @param currentSolutionIndex The current timeStep index
   /// @param timestep The length of each time step
   /// @param sizeG_A the size of the A portion of G, marks the end of the equiv
   /// currents
   virtual void updateStoredState( const Matrix< T > & solutionMatrix,
                                   const size_t currentSolutionIndex, T timestep,
                                   size_t numCurrents ) {
   }

   /// @brief adds this component's DC stamp to the target stamp.
   ///
   /// @param destination The stamp to be added to.
   /// @param solutionMatrix A vector containing all past solutions to the circuit
   /// @param numCurrents The number of currents used by the transient simulation
   virtual void addDCAnalysisStampTo( Stamp< T > & destination,
                                    const Matrix< T > & solutionVector,
                                    size_t numCurrents ) const {
      throw std::exception( "not implemented" );
   }

   /// @brief a function to update the stored state of a component based on a DC value
   ///
   /// @param solutionVector The DC vector
   /// @param sizeG_A the number of voltages
   /// @param numCurrents the number of transient currents
   virtual void updateDCStoredState( const Matrix< T > & solutionVector,
                                     size_t sizeG_A, size_t numCurrents ) {

   }

   /// @brief initialises the component
   ///
   /// @param timestep The length of each time step
   virtual void setTimestep( T timestep ) {
   }

   /// @brief Called as a helper to add the component to the elements class.
   ///
   /// @param line The line to be parsed.
   /// @param elements A reference to the elements object the component is being
   ///                 added to.
   /// @param numNodes A reference to the current max number of nodes.
   /// @param numCurrents A reference to the current max number of currents.
   /// @param numDCCurrents A reference to the current max number of DC currents.
   static void
   addToElements( const std::string & line, CircuitElements< T > & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents ) {
      throw std::exception( "not implemented" );
   }

   virtual ~Component() {};
};

#endif
