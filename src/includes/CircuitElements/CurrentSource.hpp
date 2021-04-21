#ifndef _CURRENTSOURCE_HPP_INC_
#define _CURRENTSOURCE_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"

/// @brief An ideal current source
///
/// @tparam T The value type
template< typename T >
struct CurrentSource : public Component< T > {
public:
   T value = 0;

   size_t n1 = 0;
   size_t n2 = 0;

   void addStaticStampTo( Stamp< T > & stamp ) const {
      // the p means prime ( ' ) and is used for the index - 1
      size_t n1p = n1 - 1;
      size_t n2p = n2 - 1;

      if ( n1 ) {
         stamp.s( n1p, 0 ) += -value;
      }

      if ( n2 ) {
         stamp.s( n2p, 0 ) += value;
      }
   }

   void addDCAnalysisStampTo( Stamp< T > & stamp,
                              const Matrix< T > & solutionVector,
                              size_t numCurrents ) const {
      addStaticStampTo( stamp );
   }

   static void
   addToElements( const std::string & line, CircuitElements< T > & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents ) {
      // std::regex currentSourceRegex( R"(^I(.*?)\s(\d+?)\s(\d+?)\s(.+?)\s?$)" );
      std::regex currentSourceRegex = generateRegex( "I", "n n w" );
      CurrentSource< T > currentSource;
      std::smatch matches;

      std::regex_match( line, matches, currentSourceRegex );
      currentSource.n1 = std::stoi( matches.str( 2 ) );
      currentSource.n2 = std::stoi( matches.str( 3 ) );

      numNodes = std::max( numNodes, std::stoull( matches.str( 2 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 3 ) ) );

      if constexpr ( std::is_same_v< T, double > || std::is_same_v< T, float > ) {
         currentSource.value = std::stod( matches.str( 4 ) );
      } else {
         static_assert( "Unsupported Type" );
      }

      elements.staticElements.emplace_back(
          std::make_shared< CurrentSource< T > >( currentSource ) );
   }
};

#endif
