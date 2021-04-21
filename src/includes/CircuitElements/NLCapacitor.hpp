#ifndef _NLCAPACITOR_HPP_INC_
#define _NLCAPACITOR_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"
#include <math.h>

/// @brief a non-linear capacitor model of the form
/// C = C_p + C_o * ( 1.0 + tanh( P_10 + P_11 * u ) )
/// @tparam T the value type
template< typename T >
struct NLCapacitor : public Component< T > {
public:
   size_t n1 = 0;
   size_t n2 = 0;

   T C_p = 0;
   T C_o = 0;
   T P_10 = 0;
   T P_11 = 0;

   T u_last = 0;
   T i_last = 0;

   T C_last = C_p + C_o * ( 1.0 + std::tanh( P_10 + P_11 * u_last ) );

   void
   addNonLinearStampTo( Stamp< T > & stamp,
                        const Matrix< T > & solutionMatrix,
                        const size_t currentSolutionIndex, T timestep = 0 ) const {
      const size_t n1p = n1 - 1;
      const size_t n2p = n2 - 1;

      T u = 0;

      if ( n1 > 0 ) {
         u = solutionMatrix( n1p, currentSolutionIndex );
      }

      if ( n2 > 0 ) {
         u -= solutionMatrix( n2p, currentSolutionIndex );
      }

      T C = C_p + C_o * ( 1.0 + std::tanh( P_10 + P_11 * u ) );

      T dC = C_o*P_11 / std::pow( std::cosh( P_10 + P_11 * u ), 2 );

      T i = C * ( 2.0 * ( u - u_last ) / timestep - i_last / C_last );

      T di = dC * ( 2.0 * ( u - u_last ) / timestep - i_last / C_last ) +
             2.0 * C / timestep;

      T G_eq = di;

      T I_eq = -G_eq * u + i;

      if ( n1 > 0 ) {
         stamp.G( n1p, n1p ) += G_eq;
         stamp.s( n1p, 0 ) +=  -I_eq;

         if ( n2 > 0 ) {
            stamp.G( n1p, n2p ) += -G_eq;
         }
      }

      if ( n2 > 0 ) {
         stamp.G( n2p, n2p ) += G_eq;
         stamp.s( n2p, 0 ) +=  +I_eq;

         if ( n1 > 0 ) {
            stamp.G( n2p, n1p ) += -G_eq;
         }
      }
   }

   void updateStoredState( const Matrix< T > & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A ) {
      const size_t n1p = n1 - 1;
      const size_t n2p = n2 - 1;

      T u = 0;

      if ( n1 > 0 ) {
         u = solutionMatrix( n1p , currentSolutionIndex  );
      }

      if ( n2 > 0 ) {
         u -= solutionMatrix( n2p , currentSolutionIndex  );
      }

      T C = C_p + C_o * ( 1.0 + std::tanh( P_10 + P_11 * u ) );

      i_last = C * ( 2.0 * ( u - u_last ) / timestep - i_last / C_last );

      C_last = C;

      u_last = u;
   }

   void updateDCStoredState( const Matrix< T > & solutionVector,
                             size_t sizeG_A,
                             size_t numCurrents ) {
      const size_t n1p = n1 - 1;
      const size_t n2p = n2 - 1;

      T u = 0;

      if ( n1 > 0 ) {
         u = solutionVector( n1p , 0 );
      }

      if ( n2 > 0 ) {
         u -= solutionVector( n2p , 0 );
      }

      T C = C_p + C_o * ( 1.0 + std::tanh( P_10 + P_11 * u ) );

      i_last = 0;

      C_last = C;

      u_last = u;
   }

   void addDCAnalysisStampTo( Stamp< T > & stamp,
                              const Matrix< T > & solutionVector,
                              size_t numCurrents ) const {
      // open circuit
      if ( n1 > 0 ) {
         stamp.G( n1 - 1, n1 - 1 ) += 1e-9;
      }
      if ( n2 > 0 ) {
         stamp.G( n2 - 1, n2 - 1 ) += 1e-9;
      }
   }

   static void
   addToElements( const std::string & line, CircuitElements< T > & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents ) {
      std::regex capacitorRegex = generateRegex( "CN", "n n w w w w" );
      NLCapacitor< T > cap;
      std::smatch matches;

      std::regex_match( line, matches, capacitorRegex );

      cap.designator = "CN";
      cap.designator += matches.str( 1 );

      cap.n1 = std::stoi( matches.str( 2 ) );
      cap.n2 = std::stoi( matches.str( 3 ) );

      numNodes = std::max( numNodes, std::stoull( matches.str( 2 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 3 ) ) );

      if constexpr ( std::is_same_v< T, double > || std::is_same_v< T, float > ) {
         cap.C_p =  std::stod( matches.str( 4 ) );
         cap.C_o = std::stod( matches.str( 5 ) );
         cap.P_10 = std::stod( matches.str( 6 ) );
         cap.P_11 = std::stod( matches.str( 7 ) );
         cap.C_last = cap.C_p + cap.C_o * ( 1.0 + std::tanh( cap.P_10 + cap.P_11 * cap.u_last ) );
      } else {
         static_assert( "Unsupported Type" );
      }


      elements.nonLinearElements.emplace_back(
          std::make_shared< NLCapacitor< T > >( cap ) );
      elements.nodeComponentMap.insert(
          { { cap.n1, elements.nonLinearElements.back() },
            { cap.n2, elements.nonLinearElements.back() } } );
   }
};
#endif
