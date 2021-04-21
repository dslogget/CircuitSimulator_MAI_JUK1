#ifndef _MATRIX_HPP_INC_
#define _MATRIX_HPP_INC_
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <concepts>
#include <array>
#include <assert.h>
#include <complex.h>

template< typename T, size_t N = 1 >
using storageType = typename std::conditional<
    N * sizeof( T ) < 20000, std::array< T, N >, std::vector< T > >::type;

template< typename T >
std::complex< double >
operator>( const std::complex< T > & a, const std::complex< T > & b ) {
   return std::norm( a ) > std::norm( b );
}

template< typename T >
std::complex< double >
operator<( const std::complex< T > & a, const std::complex< T > & b ) {
   return std::norm( a ) < std::norm( b );
}

template< typename T >
std::complex< double >
operator>=( const std::complex< T > & a, const std::complex< T > & b ) {
   return std::norm( a ) >= std::norm( b );
}

template< typename T >
std::complex< double >
operator<=( const std::complex< T > & a, const std::complex< T > & b ) {
   return std::norm( a ) <= std::norm( b );
}

template< typename T >
concept arithmetic = requires( T a, T b ) {
   { a == b };
   { a != b };
   { a >= b };
   { a <= b };
   { a * b };
   { a / b };
   { a + b };
   { a - b };
};

/// @brief A compile-time sized matrix row
///
/// @tparam T The value type
/// @tparam N The row length
/// @tparam ST the type of storage used (vector vs array)
template< typename T, size_t N, typename ST = storageType< T, N > >
requires arithmetic< T > struct StaticRow {
   ST columns;

   static constexpr size_t sizeN = N;

   StaticRow() {
      if constexpr ( std::is_same< ST, std::vector< T > >::value ) {
         columns.resize( N );
      }
   }

   StaticRow( T initialValue ) {
      if constexpr ( std::is_same< ST, std::vector< T > >::value ) {
         columns.resize( N );
      }
      columns.fill( initialValue );
   }

   void fill( T value ) {
      columns.fill( value );
   }

   T & operator[]( size_t index ) {
      return columns[ index ];
   }

   const T & operator[]( size_t index ) const {
      return columns[ index ];
   }

   T dot( StaticRow< T, N > other ) {
      T toRet = 0;
      for ( size_t i = 0; i < sizeN; i++ ) {
         toRet += columns[ i ] * other.columns[ i ];
      }
      return toRet;
   }
};

template< typename T, size_t M >
struct StaticLUPair;

/// @brief A compile-time sized matrix
///
/// @tparam T The value type
/// @tparam M The matrix column length
/// @tparam N The row length
/// @tparam ST the type of storage used (vector vs array)
template< typename T, size_t M, size_t N = M,
          typename ST = storageType< StaticRow< T, N >, M > >
requires arithmetic< T > struct StaticMatrix {
   ST rows;

   StaticMatrix() {
      if constexpr ( std::is_same< ST, std::vector< StaticRow< T, N > > >::value ) {
         rows.resize( M );
      }
   }

   StaticMatrix( T initialValue ) {
      if constexpr ( std::is_same< ST, std::vector< T > >::value ) {
         rows.resize( M );
      }
      fill( initialValue );
   }

   /*
   StaticMatrix( const StaticMatrix< T, M, N, ST > & other ) {
      if constexpr ( std::is_same< ST, std::vector< T > >::value ) {
         rows.resize( M );
      }
      (*this) = other;
   }

   StaticMatrix< T, M, N, ST > & operator=( const StaticMatrix< T, M, N, ST > & other
   ) { for ( size_t m = 0; m < M; m++ ) { for ( size_t n = 0; n < N; n++ ) { rows[ m
   ][ n ] = other[ m ][ n ];
         }
      }
      return *this;
   }
   */

   static constexpr size_t sizeM = M;
   static constexpr size_t sizeN = N;

   void fill( T fillVal ) {
      for ( auto & row : rows ) {
         row.fill( fillVal );
      }
   }

   StaticRow< T, N > & operator[]( size_t index ) {
      return rows[ index ];
   }

   StaticRow< T, N > operator[]( size_t index ) const {
      return rows[ index ];
   }

   void rowAddition( size_t destinationRow, size_t sourceRow, T scalingFactor ) {
      assert( 0 <= destinationRow && destinationRow <= N );
      assert( 0 <= sourceRow && sourceRow <= M );
      for ( size_t i = 0; i < N; i++ ) {
         rows[ destinationRow ][ i ] += scalingFactor * rows[ sourceRow ][ i ];
      }
   }

   void swapRows( size_t row1, size_t row2 ) {
      std::swap( rows[ row1 ], rows[ row2 ] );
   }

   StaticMatrix< T, N, M > transpose() {
      StaticMatrix< T, N, M > toRet;
      for ( size_t m = 0; m < M; m++ ) {
         for ( size_t n = 0; n < N; n++ ) {
            toRet[ n ][ m ] = rows[ m ][ n ];
         }
      }
      return toRet;
   }

   template< size_t N2 >
   StaticMatrix< T, M, N2 > multiply( const StaticMatrix< T, N, N2 > & rhs ) const {
      StaticMatrix< T, M, N2 > toRet( 0.0 );
      multiply( rhs, toRet );
      return toRet;
   }

   template< size_t N2 >
   void multiply( const StaticMatrix< T, N, N2 > & rhs,
                  StaticMatrix< T, M, N2 > & dest ) const {
      // for simd reasons, the order of operations is slightly weird. This is to
      // minimise row changes which may cause cache misses. This is due to the fact
      // that in this model rows are represented contiguously in memory,
      // so streaming instructions and local caching can be used to our advantage
      //
      // It is also worth stating that there is strong potential for multithreading
      // here, and especially if paired with the Strassen Method for matrix
      // multiplication However for small sizes, naive multiplication can be better
      // due to less memory allocations

      for ( size_t m = 0; m < M; m++ ) {
         for ( size_t k = 0; k < N; k++ ) {
            for ( size_t n = 0; n2 < N2; n2++ ) {
               destination[ m ][ n ] += rows[ m ][ k ] * rows[ k ][ n ];
            }
         }
      }
   }

   StaticMatrix< T, M, N > add( const StaticMatrix< T, M, N > & rhs ) const {
      StaticMatrix< T, M, N > toRet();
      add( rhs, toRet );
      return toRet;
   }

   void
   add( const StaticMatrix< T, M, N > & rhs, StaticMatrix< T, M, N > & dest ) const {
      for ( size_t m = 0; m < M; m++ ) {
         for ( size_t n = 0; n < N; n++ ) {
            dest[ m ][ n ] = rows[ m ][ n ] + rhs[ m ][ n ];
         }
      }
   }

   StaticMatrix< T, M, N > subtract( const StaticMatrix< T, N, N > & rhs ) const {
      StaticMatrix< T, M, N > toRet();
      add( rhs, toRet );
      return toRet;
   }

   void subtract( const StaticMatrix< T, N, N > & rhs,
                  StaticMatrix< T, M, N > & dest ) const {
      for ( size_t m = 0; m < M; m++ ) {
         for ( size_t n = 0; n < N; n++ ) {
            dest[ m ][ n ] = rows[ m ][ n ] - dest[ m ][ n ]
         }
      }
   }

   std::string toString() const {
      std::stringstream toRet;
      for ( const auto & row : rows ) {
         for ( const auto & val : row.columns ) {
            toRet << std::setw( 5 ) << std::setprecision( 2 ) << val << " ";
         }

         toRet << std::endl;
      }

      return toRet.str();
   }

   StaticLUPair< T, M > luPair() const {
      static_assert( N == M, "StaticMatrix must be a square" );

      StaticLUPair< T, M > toRet;
      luPair( toRet );
      return toRet;
   }

   void luPair( StaticLUPair< T, M > & dest ) const {
      static_assert( N == M, "StaticMatrix must be a square" );

      dest.u = *this;
      dest.l.fill( 0.0 );
      for ( size_t n = 0; n < N; n++ ) {
         dest.l[ n ][ n ] = 1.0;
         dest.p[ n ] = n;
      }

      for ( size_t r = 0; r < M - 1; r++ ) {
         // find largest in column
         size_t largestRow = r;
         auto maxV = abs( dest.u[ r ][ r ] );
         for ( size_t r2 = r + 1; r2 < M; r2++ ) {
            if ( abs( dest.u[ r2 ][ r ] ) > maxV ) {
               maxV = abs( dest.u[ r2 ][ r ] );
               largestRow = r2;
            }
         }

         // swap rows in U and indices in p
         dest.u.swapRows( r, largestRow );
         std::swap( dest.p[ r ], dest.p[ largestRow ] );
         // swap subdiagonal entries in L
         for ( size_t n = 0; n < r; n++ ) {
            std::swap( dest.l[ r ][ n ], dest.l[ largestRow ][ n ] );
         }

         // Gaussian elimination
         for ( size_t m = r + 1; m < M; m++ ) {
            // TODO: potential for multithreading here
            // need to take into account how many processors there are
            // https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency
            T multiplier = dest.u[ m ][ r ] / dest.u[ r ][ r ];
            dest.u.rowAddition( m, r, -multiplier );
            dest.l[ m ][ r ] = multiplier;
         }
         // std::cout << "After Gaussian\n " << dest.toString();
      }
   }

   StaticMatrix< T, M, 1 > leftDivide( const StaticMatrix< T, M, 1 > & rhs ) const {
      auto lu = luPair();
      StaticMatrix< T, M, 1 > scratchSpace;
      StaticMatrix< T, M, 1 > toRet;
      leftDivide( rhs, lu, scratchSpace, toRet );
      return toRet;
   }

   void
   leftDivide( const StaticMatrix< T, M, 1 > & rhs, const StaticLUPair< T, M > & lu,
               StaticMatrix< T, M, 1 > & scratchSpace,
               StaticMatrix< T, M, 1 > & dest ) const {
      for ( size_t i = 0; i < M; i++ ) {
         dest[ i ] = rhs[ lu.p[ i ] ];
      }

      // scratchSpace: solve LY = Pb for y using substitution
      for ( size_t m = 0; m < M; m++ ) {
         T val = dest[ m ][ 0 ];
         for ( size_t n = 0; n < m; n++ ) {
            val -= scratchSpace[ n ][ 0 ] * lu.l[ m ][ n ];
         }
         scratchSpace[ m ][ 0 ] = val / lu.l[ m ][ m ];
      }

      // stage2: solve Ux = Y for x using substitution
      for ( size_t m = 0; m < M; m++ ) {
         T val = scratchSpace[ M - m - 1 ][ 0 ];
         for ( size_t n = 0; n < m; n++ ) {
            val -= dest[ M - n - 1 ][ 0 ] * lu.u[ M - m - 1 ][ N - n - 1 ];
         }
         dest[ M - m - 1 ][ 0 ] = val / lu.u[ M - m - 1 ][ M - m - 1 ];
      }
   }
};

/// @brief A compile-time sized L U and pivot grouping
///
/// @tparam T The value type
/// @tparam M The L and U sizes
template< typename T, size_t M >
struct StaticLUPair {
   StaticMatrix< T, M, M > l;
   StaticMatrix< T, M, M > u;
   std::array< size_t, M > p;

   std::string toString() {
      std::stringstream toRet;
      toRet << " U\n" << u.toString();
      toRet << " L\n" << l.toString();
      toRet << " p\n";
      for ( size_t i = 0; i < p.size(); i++ ) {
         toRet << std::setw( 5 ) << p[ i ] << " ";
      }
      toRet << std::endl;
      return toRet.str();
   }
};
// --------------------------------------------------------------------------
//                   Implementation
// --------------------------------------------------------------------------


#endif
