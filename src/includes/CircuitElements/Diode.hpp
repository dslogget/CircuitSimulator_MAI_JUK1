#ifndef _DIODE_HPP_INC_
#define _DIODE_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"


/// @brief An ebbers moll diode model
///
/// @tparam T Value type
template< typename T >
struct Diode : public Component< T > {
public:
   size_t n1 = 0;
   size_t n2 = 0;

   const T I_sat = 2.52e-9;
   const T V_T = 25.8563e-3;
   const T eta = 2;

   T V_crit = eta * V_T * std::log( eta * V_T / ( I_sat * std::sqrt( 2 ) ) );


   void
   addNonLinearStampTo( Stamp< T > & stamp,
                        const Matrix< T > & solutionMatrix,
                        const size_t currentSolutionIndex, T timestep = 0 ) const {
      const size_t n1p = n1 - 1;
      const size_t n2p = n2 - 1;

      T v = 0;

      if ( n1 > 0 ) {
         v = solutionMatrix( n1p , currentSolutionIndex  );
      }

      if ( n2 > 0 ) {
         v -= solutionMatrix( n2p , currentSolutionIndex  );
      }

      v = std::min( V_crit, v );

      T G_eq = ( I_sat / ( eta * V_T ) ) * std::exp( v / ( eta * V_T ) );
      T I_eq = I_sat * ( std::exp( v / ( eta * V_T ) ) - 1 ) - G_eq * v;

      if ( n1 > 0 ) {
         stamp.G( n1p, n1p ) += G_eq;
         stamp.s( n1p, 0 ) += -I_eq;
      }

      if ( n2 > 0 ) {
         stamp.G( n2p, n2p ) += G_eq;
         stamp.s( n2p, 0 ) += +I_eq;
      }

      if ( n1 > 0 && n2 > 0 ) {
         stamp.G( n1p, n2p ) += -G_eq;
         stamp.G( n2p, n1p ) += -G_eq;
      }
   }

   void updateStoredState( const Matrix< T > & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A ) {
   }

   void addDCAnalysisStampTo( Stamp< T > & stamp,
                              const Matrix< T > & solutionVector,
                              size_t numCurrents ) const {
      addNonLinearStampTo( stamp, solutionVector, 0, 0 );
   }

   static void
   addToElements( const std::string & line, CircuitElements< T > & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents ) {
      std::regex diodeRegex = generateRegex( "D", "n n" );
      Diode< T > diode;
      std::smatch matches;

      std::regex_match( line, matches, diodeRegex );

      diode.designator = "D";
      diode.designator += matches.str( 1 );

      diode.n1 = std::stoi( matches.str( 2 ) );
      diode.n2 = std::stoi( matches.str( 3 ) );

      numNodes = std::max( numNodes, std::stoull( matches.str( 2 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 3 ) ) );

      elements.nonLinearElements.emplace_back(
          std::make_shared< Diode< T > >( diode ) );
      elements.nodeComponentMap.insert(
          { { diode.n1, elements.nonLinearElements.back() },
            { diode.n2, elements.nonLinearElements.back() } } );
   }
};

#endif
