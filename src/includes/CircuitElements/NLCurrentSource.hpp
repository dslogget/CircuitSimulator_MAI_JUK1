#ifndef _NLCURRENTSOURCE_HPP_INC_
#define _NLCURRENTSOURCE_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"
#include "Maths/AutoDifferentiation.hpp"

/// @brief An non-linear current source for the COBRA transistor model
///
/// @tparam T The value type
template< typename T >
struct NLCurrentSource : public Component< T > {
public:
   T value = 0;

   size_t n1 = 0;
   size_t n2 = 0;

   size_t r1_pos = 0;
   size_t r1_neg = 0;
   size_t r2_pos = 0;
   size_t r2_neg = 0;

   void addNonLinearStampTo( Stamp< T > & stamp,
                             const Matrix< T > & solutionMatrix,
                             const size_t currentSolutionIndex, T timestep = 0 ) const {
      constexpr T alpha = 1.3;
      constexpr T beta0 = 0.42;
      constexpr T gamma = 0.0005;
      constexpr T delta = 0.3;
      constexpr T xi = 0.06;
      constexpr T lambda = 1.5;
      constexpr T mu = 0.0;
      constexpr T zeta = 0.18;
      constexpr T Vto = -2.4;

      size_t n1p = n1 - 1;
      size_t n2p = n2 - 1;
      size_t r1p_pos = r1_pos - 1;
      size_t r1p_neg = r1_neg - 1;
      size_t r2p_pos = r2_pos - 1;
      size_t r2p_neg = r2_neg - 1;

      T u = 0;
      T r1 = 0;
      T r2 = 0;

      if ( n1 > 0 ) {
         u = solutionMatrix( n1p, currentSolutionIndex );
      }

      if ( n2 > 0 ) {
         u -= solutionMatrix( n2p, currentSolutionIndex );
      }

      if ( r1_pos > 0 ) {
         r1 = solutionMatrix( r1p_pos, currentSolutionIndex );
      }

      if ( r1_neg > 0 ) {
         r1 -= solutionMatrix( r1p_neg, currentSolutionIndex );
      }

      if ( r2_pos > 0 ) {
         r2 = solutionMatrix( r2p_pos, currentSolutionIndex );
      }

      if ( r2_neg > 0 ) {
         r2 -= solutionMatrix( r2p_neg, currentSolutionIndex );
      }

      namespace AD = AutoDifferentiation;

      using ADT = AD::DiffVar< T, 2 >;
      ADT V_gs( r1, 1, 0 );
      ADT V_ds( r2, 0, 1 );

      auto beta = beta0;
      auto Vgst = V_gs - ( 1 + beta*beta ) * Vto + gamma * V_ds;
      auto Veff = 0.5 * ( Vgst + AD::sqrt( AD::pow( Vgst, 2 ) + delta*delta ) );
      auto power = lambda / ( 1 + mu * AD::pow( V_ds, 2 ) + xi * Veff );
      auto area = alpha * V_ds * ( 1 + zeta * Veff );
      auto f1 = AD::tanh( area );
      auto Ids_lim = beta * AD::pow( Veff, power );
      auto Idrain = Ids_lim * f1;
      auto I_ds = Idrain[ 0 ] - Idrain[ 1 ] * r1 - Idrain[ 2 ] * r2;

      if ( n1 > 0 ) {
         stamp.s( n1p, 0 ) += - I_ds;
         if ( r1_pos > 0 ) {
            stamp.G( n1p, r1p_pos ) += Idrain[ 1 ];
         }
         if ( r1_neg > 0 ) {
            stamp.G( n1p, r1p_neg ) += - Idrain[ 1 ];
         }
         if ( r2_pos > 0 ) {
            stamp.G( n1p, r2p_pos ) += Idrain[ 2 ];
         }
         if ( r2_neg > 0 ) {
            stamp.G( n1p, r2p_neg ) += - Idrain[ 2 ];
         }
      }

      if ( n2 > 0 ) {
         stamp.s( n2p, 0 ) += + I_ds;
         if ( r1_pos > 0 ) {
            stamp.G( n2p, r1p_pos ) += - Idrain[ 1 ];
         }
         if ( r1_neg > 0 ) {
            stamp.G( n2p, r1p_neg ) += Idrain[ 1 ];
         }
         if ( r2_pos > 0 ) {
            stamp.G( n2p, r2p_pos ) += - Idrain[ 2 ];
         }
         if ( r2_neg > 0 ) {
            stamp.G( n2p, r2p_neg ) += Idrain[ 2 ];
         }
      }

   }

   void addDCAnalysisStampTo( Stamp< T > & stamp,
                              const Matrix< T > & solutionVector,
                              size_t numCurrents ) const {
      addNonLinearStampTo( stamp, solutionVector, 0, 0 );
   }


   static void
   addToElements( const std::string & line, CircuitElements< T > & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents ) {
      std::regex currentSourceRegex = generateRegex( "I", "n n n n n n" );
      NLCurrentSource< T > currentSource;
      std::smatch matches;

      std::regex_match( line, matches, currentSourceRegex );
      currentSource.n1 = std::stoi( matches.str( 2 ) );
      currentSource.n2 = std::stoi( matches.str( 3 ) );
      currentSource.r1_pos = std::stoi( matches.str( 4 ) );
      currentSource.r1_neg = std::stoi( matches.str( 5 ) );
      currentSource.r2_pos = std::stoi( matches.str( 6 ) );
      currentSource.r2_neg = std::stoi( matches.str( 7 ) );

      numNodes = std::max( numNodes, std::stoull( matches.str( 2 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 3 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 4 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 5 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 6 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 7 ) ) );

      elements.nonLinearElements.emplace_back(
          std::make_shared< NLCurrentSource< T > >( currentSource ) );
   }
};

#endif
