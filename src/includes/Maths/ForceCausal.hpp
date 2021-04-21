#ifndef _FORCECAUSAL_HPP_INC_
#define _FORCECAUSAL_HPP_INC_
#include "Maths/dft.hpp"

namespace ForceCausal {

// TODO: Maybe this could be made more efficient by worrying about vectorisation
// of the calculations. I.E storing precomputed F values etc

template< typename T >
std::complex< T >
F( const std::vector< T > & freq, const std::vector< std::complex< T > > & data,
   T tau, T k, size_t n ) {
   return ( data[ n ] - k ) *
          exp( std::complex< T >( 0, -2 * std::numbers::pi * freq[ n ] * tau ) );
}

template< typename T >
T
K( const std::vector< T > & freq, const std::vector< std::complex< T > > & data,
   T tau ) {
   return std::real( data.back() ) -
          std::imag( data.back() ) /
              std::tan( 2 * std::numbers::pi * freq.back() * tau );
}

template< typename T >
T
f0( const std::vector< T > & freq, const std::vector< std::complex< T > > & data,
    T tau ) {
   std::complex< T > toRet = 0;
   T k = K( freq, data, tau );
   for ( size_t i = 1; i < freq.size() - 1; i++ ) {
      toRet += 2 * ( std::real( F( freq, data, tau, k, i ) ) );
   }
   toRet += F( freq, data, tau, k, 0 );
   toRet += std::real( F( freq, data, tau, k, freq.size() - 1 ) );
   toRet *= 1e3 / ( 2 * freq.size() - 2 );
   return std::real( toRet );
}

template< typename T >
T
f0derivative( const std::vector< T > & freq,
              const std::vector< std::complex< T > > & data, T tau, T step ) {
   return ( f0( freq, data, tau + step ) - f0( freq, data, tau ) ) / step;
}


template< typename T >
T
getTau( const std::vector< T > & freq, const std::vector< std::complex< T > > & data,
        T tol = 1e-7, size_t maxIter = 30, T step = 1e-8 ) {
   T currentGuess = 1e-8;
   for ( size_t i = 0; i < maxIter; i++ ) {
      T f0Curr = f0( freq, data, currentGuess );
      if ( f0Curr * f0Curr < tol ) {
         break;
      }
      T diff = ( f0( freq, data, currentGuess + step ) - f0Curr ) / step;
      currentGuess = currentGuess - f0Curr / diff;
   }
   return currentGuess;
}

/// @brief a helper struct to return the result of the forced causal IDFT
///
/// @tparam T
template< typename T >
struct CausalData {
   T tau;
   T Ts;
   std::vector< T > data;
};

} // namespace ForceCausal

template< typename T >
ForceCausal::CausalData< T >
forceCausal( const std::vector< T > & freq,
             const std::vector< std::complex< T > > & data ) {
   ForceCausal::CausalData< T > toRet;
   toRet.data = std::vector< T >( 2 * freq.size() - 2 );
   toRet.Ts = 1.0 / ( toRet.data.size() * ( freq[ 1 ] - freq[ 0 ] ) );
   // Add in conjugate Symmetric Data
   std::vector< std::complex< T > > hermitianData( 2 * freq.size() - 2 );
   T k = 0;

   if ( abs( std::imag( data.back() ) ) < 1e-5 ) {
      toRet.tau = 0;

      for ( size_t i = 0; i < freq.size() - 1; i++ ) {
         hermitianData[ i ] = data[ i ];
      }

      for ( size_t i = 1; i < freq.size(); i++ ) {
         hermitianData[ hermitianData.size() - i ] = std::conj( data[ i ] );
      }

   } else {
      toRet.tau = ForceCausal::getTau( freq, data );
      k = ForceCausal::K( freq, data, toRet.tau );

      for ( size_t i = 0; i < freq.size() - 1; i++ ) {
         hermitianData[ i ] = ForceCausal::F( freq, data, toRet.tau, k, i );
      }

      for ( size_t i = 1; i < freq.size(); i++ ) {
         hermitianData[ hermitianData.size() - i ] = std::conj(
             ForceCausal::F( freq, data, toRet.tau, k, i ) );
      }
   }

   auto idftVal = idft( hermitianData );
   for ( size_t i = 0; i < idftVal.size(); i++ ) {
      toRet.data[ i ] = std::real( idftVal[ i ] );
   }

   if ( abs( std::imag( data.back() ) ) >= 1e-5 ) {
      toRet.data[ 0 ] = k;
   }
   return toRet;
}


#endif
