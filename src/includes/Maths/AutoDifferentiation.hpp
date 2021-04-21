#ifndef _AUTODIFFERENTIATION_HPP_INC_
#define _AUTODIFFERENTIATION_HPP_INC_
#include <array>
#include <concepts>
#include <tuple>
#include <cmath>
#include <functional>
#include <iostream>

/// @brief a namespace to hold the messiness of my auto-differentiator
///
/// In essence the idea of this header is to allow easy determination of the jacobian
/// of a non-linear element without having to manually differentiate a function.
/// Sometimes this can also lead to faster implementations than an analytical derivative.
/// The main advantage though, is being able to express the equation once in a natural
/// format.
namespace AutoDifferentiation {

template< typename RT, typename LT, typename OT >
concept MultipliableResult =
requires ( LT lhs, OT rhs ) { std::is_same< RT, decltype( lhs *= rhs ) >::value; };

template< typename RT, typename LT, typename OT >
concept AddableResult =
requires ( LT lhs, OT rhs ) { std::is_same< RT, decltype( lhs += rhs ) >::value; };

template< typename RT, typename LT, typename OT >
concept SubtractableResult =
requires ( LT lhs, OT rhs ) { std::is_same< RT, decltype( lhs -= rhs ) >::value; };

template< typename RT, typename LT, typename OT >
concept DivisableResult =
requires ( LT lhs, OT rhs ) { std::is_same< RT, decltype( lhs /= rhs ) >::value; };


template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
struct DiffVar {

   using ThisDiff = DiffVar< ValType, NumVars >;

   ValType var = 0;
   std::array< ValType, NumVars > diffVars = { 0 };

   template< typename ... VariadicType >
      requires ( ... && std::is_convertible< VariadicType, ValType >::value )
   constexpr
   explicit DiffVar( ValType _var, VariadicType ... _diffVars )
      : var( _var ), diffVars{{ static_cast< ValType >( _diffVars )... }} {
   }

   constexpr
   DiffVar( ValType _var, std::array< ValType, NumVars > _diffVars )
      : var( _var ), diffVars( _diffVars ) {
   }

   template< typename OtherValType >
      requires ( std::is_convertible< OtherValType, ValType >::value )
   constexpr
   explicit DiffVar( const OtherValType & other ) {
      var = static_cast< ValType >( other );
   }

   constexpr
   explicit operator ValType() const {
      return var;
   }

   constexpr
   ValType & operator[]( size_t index ) {
      if ( index == 0 ) {
         return var;
      } else {
         return diffVars[ index - 1 ];
      }
   }

   constexpr
   const ValType & operator[]( size_t index ) const {
      if ( index == 0 ) {
         return var;
      } else {
         return diffVars[ index - 1 ];
      }
   }

   template< typename OtherValType >
      requires ( std::is_convertible< OtherValType, ValType >::value )
   constexpr
   int operator<=>( OtherValType & rhs ) {
      return var <=> static_cast< ValType >( rhs );
   }

   constexpr
   ThisDiff operator+=( const ThisDiff & rhs ) {
      var += rhs.var;
      auto it1 = diffVars.begin();
      auto it2 = rhs.diffVars.begin();
      while ( it1 != diffVars.end() && it2 != rhs.diffVars.end() ) {
         *it1 += *it2;
         ++it1;
         ++it2;
      }
      return *this;
   }

   template< typename OtherValType >
      requires ( std::is_convertible< OtherValType, ValType >::value )
   constexpr
   ThisDiff operator+=( const OtherValType & rhs ) {
      var += static_cast< ValType >( rhs );
      return *this;
   }

   constexpr
   ThisDiff operator-=( const ThisDiff & rhs ) {
      var -= rhs.var;
      auto it1 = diffVars.begin();
      auto it2 = rhs.diffVars.begin();
      while ( it1 != diffVars.end() && it2 != rhs.diffVars.end() ) {
         *it1 -= *it2;
         ++it1;
         ++it2;
      }
      return *this;
   }

   template< typename OtherValType >
      requires ( std::is_convertible< OtherValType, ValType >::value )
   constexpr
   ThisDiff operator-=( const OtherValType & rhs ) {
      var -= static_cast< ValType >( rhs );
      return *this;
   }

   constexpr
   ThisDiff operator-() const {
      auto toRet = *this;
      toRet.var = -var;
      auto it1 = toRet.diffVars.begin();
      while ( it1 != toRet.diffVars.end() ) {
         *it1 = -*it1;
         ++it1;
      }
      return toRet;
   }


   constexpr
   ThisDiff operator*=( const ThisDiff & rhs ) {
      auto it1 = diffVars.begin();
      auto it2 = rhs.diffVars.begin();
      while ( it1 != diffVars.end() && it2 != rhs.diffVars.end() ) {
         *it1 = ( rhs.var * *it1 + var * *it2 );
         ++it1;
         ++it2;
      }
      var *= rhs.var;
      return *this;
   }

   template< typename OtherValType >
      requires ( std::is_convertible< OtherValType, ValType >::value )
   constexpr
   ThisDiff operator*=( const OtherValType & rhs ) {
      auto it1 = diffVars.begin();
      while ( it1 != diffVars.end() ) {
         *it1 = static_cast< ValType >( rhs ) * *it1;
         ++it1;
      }
      var *= static_cast< ValType >( rhs );
      return *this;
   }

   constexpr
   ThisDiff operator/=( const ThisDiff & rhs ) {
      auto it1 = diffVars.begin();
      auto it2 = rhs.diffVars.begin();
      while ( it1 != diffVars.end() && it2 != rhs.diffVars.end() ) {
         *it1 = ( rhs.var * *it1 - var * *it2 ) / ( rhs.var * rhs.var );
         ++it1;
         ++it2;
      }
      var /= rhs.var;
      return *this;
   }

   template< typename OtherValType >
      requires ( std::is_convertible< OtherValType, ValType >::value )
   constexpr
   ThisDiff operator/=( const OtherValType & rhs ) {
      auto it1 = diffVars.begin();
      while ( it1 != diffVars.end() ) {
         *it1 =  *it1 / static_cast< ValType >( rhs );
         ++it1;
      }
      var /= static_cast< ValType >( rhs );
      return *this;
   }

   // Friends
   //    Binary Operators
   //       Addition
   template< typename LT, typename RT >
      requires AddableResult< ThisDiff, LT, RT >
   constexpr
   friend auto operator+( LT lhs,
                          const RT & rhs ) {
      lhs += rhs;
      return lhs;
   }

   template< typename LT, typename RT >
      requires AddableResult< ThisDiff, RT, LT > &&
               ( !AddableResult< ThisDiff, LT, RT > )
   constexpr
   friend auto operator+( const LT & lhs,
                          RT rhs ) {
      rhs += lhs;
      return rhs;
   }

   //       Subtraction
   template< typename LT, typename RT >
      requires SubtractableResult< ThisDiff, LT, RT >
   constexpr
   friend auto operator-( LT lhs,
                          const RT & rhs ) {
      lhs -= rhs;
      return lhs;
   }

   template< typename LT, typename RT >
      requires SubtractableResult< ThisDiff, RT, LT > &&
               ( !SubtractableResult< ThisDiff, LT, RT > )
   constexpr
   friend auto operator-( const LT & lhs,
                          RT rhs ) {
      rhs -= lhs;
      return -rhs;
   }

   template< typename LT, typename RT >
      requires ( !SubtractableResult< ThisDiff, LT, RT > ) &&
               AddableResult< ThisDiff, LT, RT >
   constexpr
   friend auto operator-( LT lhs,
                          const RT & rhs ) {
      lhs += -rhs;
      return lhs;
   }

   template< typename LT, typename RT >
      requires ( !SubtractableResult< ThisDiff, RT, LT > ) &&
               AddableResult< ThisDiff, RT, LT > &&
               ( !AddableResult< ThisDiff, LT, RT > )
   constexpr
   friend auto operator-( const LT & lhs,
                          RT rhs ) {
      rhs += -lhs;
      return -rhs;
   }

   //       Multiplication
   template< typename LT, typename RT >
      requires MultipliableResult< ThisDiff, LT, RT >
   constexpr
   friend auto operator*( LT lhs,
                          const RT & rhs ) {
      lhs *= rhs;
      return lhs;
   }

   template< typename LT, typename RT >
      requires MultipliableResult< ThisDiff, RT, LT > &&
               ( !MultipliableResult< ThisDiff, LT, RT > )
   constexpr
   friend auto operator*( const LT & lhs,
                          RT rhs ) {
      rhs *= lhs;
      return rhs;
   }

   //       Division
   template< typename OT >
      requires ( !DivisableResult< ThisDiff,
                                   OT,
                                   ThisDiff > ) &&
               ( std::is_convertible< OT, ValType >::value )
   constexpr
   friend auto operator/( const OT & lhs,
                          const ThisDiff & rhs ) {
      auto lhsDiff = static_cast< ThisDiff >( lhs );
      lhsDiff /= rhs;
      return lhsDiff;
   }

   template< typename OT >
      requires ( std::is_convertible< OT, ValType >::value )
   constexpr
   friend auto operator/( ThisDiff lhs,
                          const OT & rhs ) {
      lhs /= rhs;
      return lhs;
   }

   constexpr
   friend std::ostream & operator<<( std::ostream & output, const ThisDiff & v ) {
      output << "Value: " << v.var;
      size_t n = 0;
      for ( auto diffVar : v.diffVars ) {
         output << std::endl << "Derivative " << n++ << ": " << diffVar;
      }

      return output;
   }

};


template< typename ValType, size_t NumVars, typename F1, typename F2 >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
constexpr
auto diffFunc( const DiffVar< ValType, NumVars > & arg,
               F1 func,
               F2 deriv ) {
   DiffVar< ValType, NumVars > toRet( func( arg.var ) );

   auto it1 = toRet.diffVars.begin();
   auto it2 = arg.diffVars.begin();
   auto derivEval = deriv( arg.var );
   while ( it1 != toRet.diffVars.end() && it2 != arg.diffVars.end() ) {
      *it1 = *it2 * derivEval;
      ++it1;
      ++it2;
   }
   return toRet;
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto sin( const DiffVar< ValType, NumVars > & arg ) {
   return diffFunc( arg,
                    []( ValType v ){ return std::sin( v ); },
                    []( ValType v ){ return std::cos( v ); } );
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto cos( const DiffVar< ValType, NumVars > & arg ) {
   return diffFunc( arg,
                    []( ValType v ){ return std::cos( v ); },
                    []( ValType v ){ return -std::sin( v ); } );
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto tan( const DiffVar< ValType, NumVars > & arg ) {
   return diffFunc( arg,
                    []( ValType v ){ return std::tan( v ); },
                    []( ValType v ){ return 1 / ( std::cos( v ) * std::cos( v ) ); } );
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto sinh( const DiffVar< ValType, NumVars > & arg ) {
   return diffFunc( arg,
                    []( ValType v ){ return std::sinh( v ); },
                    []( ValType v ){ return std::cosh( v ); } );
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto cosh( const DiffVar< ValType, NumVars > & arg ) {
   return diffFunc( arg,
                    []( ValType v ){ return std::cosh( v ); },
                    []( ValType v ){ return std::sinh( v ); } );
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto tanh( const DiffVar< ValType, NumVars > & arg ) {
   return diffFunc( arg,
                    []( ValType v ){ return std::tanh( v ); },
                    []( ValType v ){ return 1 / ( std::cosh( v ) * std::cosh( v ) ); } );
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto exp( const DiffVar< ValType, NumVars > & arg ) {
   auto derivEval = std::exp( arg.var );
   DiffVar< ValType, NumVars > toRet( derivEval);

   auto it1 = toRet.diffVars.begin();
   auto it2 = arg.diffVars.begin();
   while ( it1 != toRet.diffVars.end() && it2 != arg.diffVars.end() ) {
      *it1 = *it2 * derivEval;
      ++it1;
      ++it2;
   }
   toRet.var = derivEval;
   return toRet;
}

template< typename ValType, size_t NumVars, typename OT >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 ) &&
            ( !std::is_same< OT, DiffVar< ValType, NumVars > >::value )
auto pow( const DiffVar< ValType, NumVars > & arg, OT exponent ) {
   return diffFunc( arg,
                    [ &exponent ]( ValType v ){ return std::pow( v, exponent ); },
                    [ &exponent ]( ValType v ){ return exponent * std::pow( v, exponent - 1 ); } );
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto pow( const DiffVar< ValType, NumVars > & arg, DiffVar< ValType, NumVars > exponent ) {
   auto func = [] ( ValType v, ValType exponent ) {
      return std::pow( v, exponent );
   };

   auto deriv = [] ( ValType f, ValType fPrime, ValType g, ValType gPrime ) {
      return std::pow( f, g - 1 ) * ( g * fPrime + f * std::log( f ) * gPrime );
   };

   DiffVar< ValType, NumVars > toRet( func( arg.var, exponent.var ) );

   auto it1 = toRet.diffVars.begin();
   auto it2 = arg.diffVars.begin();
   auto it3 = exponent.diffVars.begin();
   while ( it1 != toRet.diffVars.end() && it2 != arg.diffVars.end() ) {
      *it1 = deriv( arg.var, *it2, exponent.var, *it3 );
      ++it1;
      ++it2;
      ++it3;
   }
   return toRet;
}

template< typename ValType, size_t NumVars >
   requires ( std::is_arithmetic< ValType >::value && NumVars >= 0 )
auto sqrt( const DiffVar< ValType, NumVars > & arg ) {
   auto funcResult = std::pow( arg.var, 0.5 );
   DiffVar< ValType, NumVars > toRet( funcResult );

   auto it1 = toRet.diffVars.begin();
   auto it2 = arg.diffVars.begin();
   auto derivEval = 0.5 / funcResult;
   while ( it1 != toRet.diffVars.end() && it2 != arg.diffVars.end() ) {
      *it1 = *it2 * derivEval;
      ++it1;
      ++it2;
   }
   return toRet;
}


}
#endif
