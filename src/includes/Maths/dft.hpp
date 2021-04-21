#ifndef _FFT_H_INC_
#define _FFT_H_INC_
#include <vector>
#include <complex>
#include <math.h>
#include <numbers>

template< typename T >
std::complex< T >
nthRootOfUnity( T fraction ) {
   return std::exp( std::complex< T >( 0, -2 * std::numbers::pi * fraction ) );
}

template< typename T >
std::complex< T >
nthRootOfUnity( int numerator, size_t denominator ) {
   return std::exp(
       std::complex< T >( 0, -2 * std::numbers::pi * numerator /
                                 static_cast< double >( denominator ) ) );
}


template< typename T >
std::vector< std::complex< T > >
dft( const std::vector< T > & inputData ) {
   const size_t length = inputData.size();
   std::vector< std::complex< T > > toRet( length, std::complex< T >( 0.0, 0.0 ) );
   for ( size_t k = 0; k < length; k++ ) {
      for ( size_t n = 0; n < length; n++ ) {
         toRet[ k ] += inputData[ n ] * nthRootOfUnity< T >( k * n, length );
      }
   }
   return toRet;
}

template< typename T >
std::vector< std::complex< T > >
idft( const std::vector< std::complex< T > > & inputData ) {
   std::vector< std::complex< T > > toRet;
   const size_t length = inputData.size();
   toRet.resize( length );
   for ( size_t n = 0; n < length; n++ ) {
      for ( size_t k = 0; k < length; k++ ) {
         toRet[ n ] += inputData[ k ] *
                       nthRootOfUnity< T >( -static_cast< int >( k * n ), length );
      }
      toRet[ n ] /= length;
   }
   return toRet;
}


template< typename T >
std::vector< std::complex< T > >
fft( const std::vector< T > & inputData ) {
   std::vector< std::complex< T > > result( inputData.size(),
                                            std::complex< T >( 0, 0 ) );
   std::vector< std::complex< T > > scratch( inputData.size(),
                                             std::complex< T >( 0, 0 ) );
   _fftHelperRadix2( inputData, result.begin(), scratch.begin(), 0, 0 );
   return result;
}

template< typename T >
std::vector< std::complex< T > >
ifft( const std::vector< std::complex< T > > & inputData ) {
   std::vector< std::complex< T > > result( inputData.size(),
                                            std::complex< T >( 0, 0 ) );
   std::vector< std::complex< T > > scratch( inputData.size(),
                                             std::complex< T >( 0, 0 ) );
   _fftHelperRadix2< std::complex< T >, decltype( result.begin() ), T,
                     -1 >( inputData, result.begin(), scratch.begin(), 0, 0 );
   for ( auto & num : result ) {
      num /= result.size();
   }
   return result;
}

template< typename T, typename Iter, typename U = T, int dir = 1 >
void
_fftHelperRadix2( const std::vector< T > & inputData, Iter result, Iter scratch,
                  size_t offset, size_t stride ) {
   const size_t len = ( inputData.size() >> stride );
   if ( len > 2 ) {
      _fftHelperRadix2< T, Iter, U, dir >( inputData, scratch, result, offset,
                                           stride + 1 );

      _fftHelperRadix2< T, Iter, U, dir >( inputData, scratch + ( len / 2 ),
                                           result + ( len / 2 ),
                                           offset + ( 1 << stride ), stride + 1 );

      for ( size_t i = 0; i < len / 2; i++ ) {
         result[ i ] = scratch[ i ] +
                       nthRootOfUnity< U >( dir * i, len ) * scratch[ i + len / 2 ];
         result[ i + len / 2 ] = scratch[ i ] - nthRootOfUnity< U >( dir * i, len ) *
                                                    scratch[ i + len / 2 ];
      }
   } else {
      result[ 0 ] = inputData[ offset ] + inputData[ offset + ( 1 << stride ) ];
      result[ 1 ] = inputData[ offset ] - inputData[ offset + ( 1 << stride ) ];
   }
}

#endif
