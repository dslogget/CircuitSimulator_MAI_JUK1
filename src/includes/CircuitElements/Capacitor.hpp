#ifndef _CAPACITOR_HPP_INC_
#define _CAPACITOR_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"


/// @brief An ideal capacitor model
///
/// @tparam T the value type
template< typename T >
struct Capacitor : public Component< T > {
public:
   T value = 0;

   size_t n1 = 0;
   size_t n2 = 0;
   T lastCurrent = 0;

   bool trapezoidalRule = true;
   void addDynamicStampTo( Stamp< T > & stamp,
                           const Matrix< T > & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep ) const {
      size_t n1p = n1 - 1;
      size_t n2p = n2 - 1;

      T u0 = 0;
      if ( n1 ) {
         u0 = solutionMatrix( n1p , currentSolutionIndex - 1 );
      }

      if ( n2 ) {
         u0 -= solutionMatrix( n2p , currentSolutionIndex - 1 );
      }


      T G_eq = 0;
      T I_eq = 0;

      if ( trapezoidalRule ) {
         G_eq = 2 * value / timestep;
         I_eq = lastCurrent + G_eq * u0;
      } else {
         G_eq = value / timestep;
         I_eq = value * u0 / timestep;
      }

      if ( n1 ) {
         stamp.G( n1p, n1p ) += G_eq;
         stamp.s( n1p, 0 ) += I_eq;
      }

      if ( n2 ) {
         stamp.G( n2p, n2p ) += G_eq;
         stamp.s( n2p, 0 ) += -I_eq;
      }

      if ( n1 && n2 ) {
         stamp.G( n1p, n2p ) += -G_eq;
         stamp.G( n2p, n1p ) += -G_eq;
      }
   }

   void updateStoredState( const Matrix< T > & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A ) {
      if ( trapezoidalRule ) {
         size_t n1p = n1 - 1;
         size_t n2p = n2 - 1;
         T u0 = 0;
         T u1 = 0;
         if ( n1 ) {
            u0 = solutionMatrix( n1p , currentSolutionIndex - 1 );
            u1 = solutionMatrix( n1p , currentSolutionIndex );
         }

         if ( n2 ) {
            u0 -= solutionMatrix( n2p , currentSolutionIndex - 1 );
            u1 -= solutionMatrix( n2p , currentSolutionIndex );
         }

         T G_eq = 2 * value / timestep;
         lastCurrent = G_eq * u1 - ( lastCurrent + G_eq * u0 );
      }
   }

   void addDCAnalysisStampTo( Stamp< T > & stamp,
                              const Matrix< T > & solutionVector,
                              size_t numCurrents ) const {
      // open circuit. Maybe add large connection to ref if unstable
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
      // std::regex capacitorRegex( R"(^C(.*?)\s(\d+?)\s(\d+?)\s(.+?)\s?$)" );
      std::regex capacitorRegex = generateRegex( "C", "n n w" );
      Capacitor< T > capacitor;
      std::smatch matches;

      std::regex_match( line, matches, capacitorRegex );

      capacitor.designator = "C";
      capacitor.designator += matches.str( 1 );

      capacitor.n1 = std::stoi( matches.str( 2 ) );
      capacitor.n2 = std::stoi( matches.str( 3 ) );
      capacitor.trapezoidalRule = true;

      numNodes = std::max( numNodes, std::stoull( matches.str( 2 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 3 ) ) );

      if constexpr ( std::is_same_v< T, double > || std::is_same_v< T, float > ) {
         capacitor.value = std::stod( matches.str( 4 ) );
      } else {
         static_assert( "Unsupported Type" );
      }

      elements.dynamicElements.emplace_back(
          std::make_shared< Capacitor< T > >( capacitor ) );

      elements.nodeComponentMap.insert(
          { { capacitor.n1, elements.dynamicElements.back() },
            { capacitor.n2, elements.dynamicElements.back() } } );
   }
};


#endif
