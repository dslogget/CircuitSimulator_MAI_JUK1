#ifndef _VOLTAGESOURCE_HPP_INC_
#define _VOLTAGESOURCE_HPP_INC_
#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"

/// @brief An ideal voltageSource
///
/// @tparam T The value type
template< typename T >
struct VoltageSource : public Component< T > {
public:
   T value = 0;

   size_t n1 = 0;
   size_t n2 = 0;
   size_t currentIndex = 0;

   void addStaticStampTo( Stamp< T > & stamp ) const {
      // the p means prime ( ' ) and is used for the index - 1
      size_t n1p = n1 - 1;
      size_t n2p = n2 - 1;
      size_t currentIndexp = currentIndex - 1;

      if ( n1 ) {
         stamp.G( n1p, stamp.sizeG_A + currentIndexp ) += 1;
         stamp.G( stamp.sizeG_A + currentIndexp, n1p ) += 1;
      }
      if ( n2 ) {
         stamp.G( n2p, stamp.sizeG_A + currentIndexp ) += -1;
         stamp.G( stamp.sizeG_A + currentIndexp, n2p ) += -1;
      }

      stamp.s( stamp.sizeG_A + currentIndexp, 0 ) += value;
   }

   void addDCAnalysisStampTo( Stamp< T > & stamp,
                              const Matrix< T > & solutionVector,
                              size_t numCurrents ) const {
      addStaticStampTo( stamp );
   }

   static void
   addToElements( const std::string & line, CircuitElements< T > & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents ) {
      // std::regex voltageSourceRegex( R"(^V(.*?)\s(\d+?)\s(\d+?)\s(.+?)\s?$)" );
      std::regex voltageSourceRegex = generateRegex( "V", "n n w" );
      VoltageSource< T > voltageSource;
      std::smatch matches;

      std::regex_match( line, matches, voltageSourceRegex );

      voltageSource.designator = "V";
      voltageSource.designator += matches.str( 1 );

      voltageSource.n1 = std::stoi( matches.str( 2 ) );
      voltageSource.n2 = std::stoi( matches.str( 3 ) );


      numNodes = std::max( numNodes, std::stoull( matches.str( 2 ) ) );
      numNodes = std::max( numNodes, std::stoull( matches.str( 3 ) ) );

      if constexpr ( std::is_same_v< T, double > || std::is_same_v< T, float > ) {
         voltageSource.value = std::stod( matches.str( 4 ) );
      } else {
         static_assert( "Unsupported Type" );
      }

      voltageSource.currentIndex = ++numCurrents;

      elements.staticElements.emplace_back(
          std::make_shared< VoltageSource< T > >( voltageSource ) );
      elements.nodeComponentMap.insert(
          { { voltageSource.n1, elements.staticElements.back() },
            { voltageSource.n2, elements.staticElements.back() } } );
   }
};

#endif
