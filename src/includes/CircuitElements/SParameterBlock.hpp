#ifndef _SPARAMETERBLOCK_HPP_INC_
#define _SPARAMETERBLOCK_HPP_INC_
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "CircuitElements/Component.hpp"
#include "Maths/DynamicMatrix.hpp"
#include "Maths/dft.hpp"
#include "Maths/ForceCausal.hpp"
#include <immintrin.h>


/// @brief a helper struct to store the information for the ports of an
///        S-Parameter Block
///
/// @tparam T The value type
template< typename T >
struct SParameterPort {
   /// @brief The positive index
   size_t positive = 0;
   /// @brief The negative index
   size_t negative = 0;
   /// @brief The current index
   size_t current = 0;
   /// @brief The equivalent resistance
   T R = 0;
   /// @brief a constant factor
   T beta = 0;
   /// @brief The start of each SParameterSequence.
   std::vector< T > s0;
};

struct SParamLengthOffset {
   size_t length;
   size_t offset;
};

/// @brief a helper struct to store the DTIR sequence for the bloc
///
/// @tparam T The value type
template< typename T >
struct SParameterSequence {
   std::vector< T > _data;
   std::vector< T > _time;
   std::vector< SParamLengthOffset > sParamLengthOffset;
   size_t numPorts = 0;

   size_t & length( size_t a, size_t b ) {
      return sParamLengthOffset[ a * numPorts + b ].length;
   }

   const size_t length( size_t a, size_t b ) const {
      return sParamLengthOffset[ a * numPorts + b ].length;
   }

   size_t & offset( size_t a, size_t b ) {
      return sParamLengthOffset[ a * numPorts + b ].offset;
   }

   const size_t offset( size_t a, size_t b ) const {
      return sParamLengthOffset[ a * numPorts + b ].offset;
   }

   T & data( size_t a, size_t b, size_t n ) {
      return _data[ offset( a, b ) + n ];
   }

   const T & data( size_t a, size_t b, size_t n ) const {
      return _data[ offset( a, b ) + n ];
   }

   T & time( size_t a, size_t b, size_t n ) {
      return _time[ offset( a, b ) + n ];
   }

   const T & time( size_t a, size_t b, size_t n ) const {
      return _time[ offset( a, b ) + n ];
   }
};

/// @brief A DTIR based model of an s-parameter block.
///
/// @tparam T
template< typename T >
struct SParameterBlock : public Component< T > {
   std::string touchstoneFilePath = "";
   std::vector< SParameterPort< T > > port;
   SParameterSequence< T > s;

   T z_ref = 0;
   T fracMaxToKeep = 0;

   /// @brief performs a linear interpolation and returns the a wave value for
   ///        use in the convolution.
   ///
   /// @param portIndex The port being processed
   /// @param solutionMatrix The working solution vector
   /// @param n The current time step
   /// @param sTimePoint The time of DTIR
   /// @param simulationTimestep The timestep of the simulation
   /// @param sizeG_A The size of the voltage dependant portion of the stamp
   ///
   /// @return the value for use in the convolution
   T
   aWaveConvValue( size_t portIndex,
                   const Matrix< T > & solutionMatrix, const size_t n,
                   T sTimePoint, const T simulationTimestep, size_t sizeG_A ) const {
      T kprime = sTimePoint / simulationTimestep;
      // This check is incredibly slow because it occurs with all loops. A similar
      // check could occur in another place when S-Parameters are first formulated
      // maybe
      // if ( 0 <= kprime && kprime < 1.0 ) {
      //   std::cout << "WARNING: SIMULATION TIMESTEP TOO LARGE FOR CONVOLUTION" <<
      //   std::endl;
      //}
      T index = n - kprime;
      // TODO: Maybe look into changing the bounds of the convolution to avoid having
      // this checked too much. Maybe wrap in debug ifdef.
      if ( index <= 0 ) {
         return 0.0;
      }

      size_t floor = static_cast< size_t >( index );
      if ( floor <= 0 || floor + 1 >= n ) {
         return 0.0;
      }

      T mix = index - static_cast< T >( floor );

      T toRet;

      T uppArr[ 3 ] = { 0 };
      T lowArr[ 3 ] = { 0 };

      if ( port[ portIndex ].positive != 0 ) {
         uppArr[ 0 ] = solutionMatrix( port[ portIndex ].positive - 1 , floor + 1 );
         lowArr[ 0 ] = solutionMatrix( port[ portIndex ].positive - 1 , floor );
      }

      if ( port[ portIndex ].negative != 0 ) {
         uppArr[ 1 ] = solutionMatrix( port[ portIndex ].negative - 1 , floor + 1 );
         lowArr[ 1 ] = solutionMatrix( port[ portIndex ].negative - 1 , floor );
      }

      uppArr[ 2 ] = solutionMatrix( sizeG_A + port[ portIndex ].current - 1 , floor + 1 );
      lowArr[ 2 ] = solutionMatrix( sizeG_A + port[ portIndex ].current - 1 , floor );

      toRet = ( uppArr[ 0 ] - lowArr[ 0 ] ) * mix + lowArr[ 0 ] -
              ( uppArr[ 1 ] - lowArr[ 1 ] ) * mix - lowArr[ 1 ] +
              ( ( uppArr[ 2 ] - lowArr[ 2 ] ) * mix + lowArr[ 2 ] ) * z_ref;
      return toRet;
   }


   /// @brief Determines the equivalent port voltage source by convolving the
   ///        (historic) a wave values with the DTIR
   ///
   /// @param p The port index
   /// @param solutionMatrix The solution vector of the circuit
   /// @param n The current timestep
   /// @param simulationTimestep the simulations timestep
   /// @param sizeG_A the size of the voltage dependant part of the stamp
   ///
   /// @return Equivalent port voltage source voltage
   T V_p( size_t p, const Matrix< T > & solutionMatrix,
          const size_t n, T simulationTimestep, size_t sizeG_A ) const {
      // V_p = beta * sum of ports ( history of port )
      T toRet = 0;
      // TODO: This may finally be a good place for multithreading
      // We could have one thread per port perhaps. Mainly as the data overlap is
      // minimal This would require having a different "toRet" for each thread and
      // summing at the end before returning perhaps
      for ( size_t c = 0; c < port.size(); c++ ) {
         // TODO: Optimise the convolution here
         // size_t len = std::min( n, s.sParamLength );
         for ( size_t k = 1; k < s.length( p, c ); k++ ) {
            toRet += aWaveConvValue( c, solutionMatrix, n, s.time( p, c, k ),
                                     simulationTimestep, sizeG_A ) *
                     s.data( p, c, k );
         }
      }
      return port[ p ].beta * toRet;
   }

   T R_p( size_t p ) const {
      return port[ p ].R;
   }

   T beta_p( size_t p ) const {
      return port[ p ].beta
   }

   void addStaticStampTo( Stamp< T > & stamp ) const {
      for ( size_t p = 0; p < port.size(); p++ ) {
         size_t np = port[ p ].positive - 1;
         size_t nn = port[ p ].negative - 1;
         size_t curr = port[ p ].current - 1;
         // Voltage source and resistance
         stamp.G( stamp.sizeG_A + curr, stamp.sizeG_A + curr ) += -port[ p ].R;
         if ( port[ p ].positive != 0 ) {
            stamp.G( stamp.sizeG_A + curr, np ) += 1;
            stamp.G( np, stamp.sizeG_A + curr ) += 1;
         }

         if ( port[ p ].negative != 0 ) {
            stamp.G( stamp.sizeG_A + curr, nn ) += -1;
            stamp.G( nn, stamp.sizeG_A + curr ) += -1;
         }
         // controlled sources
         for ( size_t c = 0; c < port.size(); c++ ) {
            if ( c != p ) {
               T alpha = port[ p ].beta * port[ p ].s0[ c ];
               if ( port[ c ].positive != 0 ) {
                  stamp.G( stamp.sizeG_A + curr, port[ c ].positive - 1 ) += -alpha;
               }
               if ( port[ c ].negative != 0 ) {
                  stamp.G( stamp.sizeG_A + curr, port[ c ].negative - 1 ) += alpha;
               }
               stamp.G( stamp.sizeG_A + curr,
                        stamp.sizeG_A + port[ c ].current - 1 ) += -z_ref * alpha;
            }
         }
      }
   }

   void addDynamicStampTo( Stamp< T > & stamp,
                           const Matrix< T > & solutionMatrix,
                           const size_t currentSolutionIndex,
                           T simulationTimestep ) const {
      for ( size_t p = 0; p < port.size(); p++ ) {
         size_t curr = port[ p ].current - 1;
         // V_p
         stamp.s( stamp.sizeG_A + curr,
                  0 ) += V_p( p, solutionMatrix, currentSolutionIndex,
                              simulationTimestep, stamp.sizeG_A );
      }
   }

   void addDCAnalysisStampTo( Stamp< T > & stamp,
                              const Matrix< T > & solutionVector,
                              size_t numCurrents ) const {
      for ( size_t p = 0; p < port.size(); p++ ) {
         size_t np = port[ p ].positive - 1;
         size_t nn = port[ p ].negative - 1;
         size_t curr = port[ p ].current - 1;

         // Voltage source and resistance
         T sppSum = 0;
         for ( size_t k = 0; k < s.length( p, p ); k++ ) {
            sppSum += s.data( p, p, k );
         }
         T Rprime = port[ p ].beta * z_ref * ( 1 + sppSum ) / ( 1 - port[ p ].beta * sppSum );
         stamp.G( stamp.sizeG_A + curr, stamp.sizeG_A + curr ) += Rprime;
         if ( port[ p ].positive != 0 ) {
            stamp.G( stamp.sizeG_A + curr, np ) += 1;
            stamp.G( np, stamp.sizeG_A + curr ) += 1;
         }

         if ( port[ p ].negative != 0 ) {
            stamp.G( stamp.sizeG_A + curr, nn ) += -1;
            stamp.G( nn, stamp.sizeG_A + curr ) += -1;
         }

         // controlled sources
         for ( size_t c = 0; c < port.size(); c++ ) {
            if ( c != p ) {
               T alpha = port[ p ].beta * port[ p ].s0[ c ];
               T alphaPrime = 0;

               for ( size_t k = 0; k < s.length( p, c ); k++ ) {
                  alphaPrime += s.data( p, c, k );
               }

               alphaPrime = port[ p ].beta * alphaPrime;

               alphaPrime += alpha;

               alphaPrime = alphaPrime / ( 1 - port[ p ].beta * sppSum );


               if ( port[ c ].positive != 0 ) {
                  stamp.G( stamp.sizeG_A + curr, port[ c ].positive - 1 ) += -alphaPrime;
               }
               if ( port[ c ].negative != 0 ) {
                  stamp.G( stamp.sizeG_A + curr, port[ c ].negative - 1 ) += alphaPrime;
               }
               stamp.G( stamp.sizeG_A + curr,
                        stamp.sizeG_A + port[ c ].current - 1 ) += -z_ref * alphaPrime;
            }
         }

         // V_p = 0
      }
   }

   void updateStoredState( const Matrix< T > & solutionMatrix,
                           const size_t currentSolutionIndex, T timestep,
                           size_t sizeG_A ) {
   }


   /// @brief reads in the s-parameter data from a touchstone file. Currently it
   ///        ignores the units of frequency, and only works with real-imag
   ///        formatting
   void readInTouchstoneFile() {
      using CT = std::complex< T >;
      using FSP = std::vector< CT >;
      std::ifstream file( touchstoneFilePath );

      std::vector< T > freqs;
      std::vector< std::vector< FSP > > freqSParams( s.numPorts );
      z_ref = 50;

      std::string line;
      while ( file.peek() == '#' || file.peek() == '!' ) {
         std::getline( file, line );
      }

      T val1;
      T val2;
      while ( !file.eof() ) {
         file >> val1;
         if ( file.fail() ) {
            break;
         }
         freqs.emplace_back( val1 );
         for ( size_t a = 0; a < s.numPorts; a++ ) {
            freqSParams[ a ].resize( s.numPorts );
         }

         for ( size_t b = 0; b < s.numPorts; b++ ) {
            for ( size_t a = 0; a < s.numPorts; a++ ) {
               file >> val1 >> val2;
               freqSParams[ a ][ b ].emplace_back( val1, val2 );
            }
         }
      }

      // std::vector< T > symFreqs( 2 * freqs.size() - 1 );
      // T maxFreq = freqs.back();
      /*      for ( size_t k = 0; k < freqs.size(); k++ ) {
               symFreqs[ k ] = freqs[ k ];
               symFreqs[ freqs.size() + k ] = maxFreq + freqs[ k ];
            } */

      // s.sParamLength = ( 2 * freqs.size() - 2 );
      s.sParamLengthOffset.resize( s.numPorts * s.numPorts );
      for ( size_t a = 0; a < s.numPorts; a++ ) {
         port[ a ].s0.resize( s.numPorts );
         for ( size_t b = 0; b < s.numPorts; b++ ) {
            auto causal = forceCausal( freqs, freqSParams[ a ][ b ] );

            T thresholdToKeep = 1;
            for ( auto entry : causal.data ) {
               thresholdToKeep = std::max( std::abs( entry ), thresholdToKeep );
            }

            thresholdToKeep = thresholdToKeep * fracMaxToKeep;

            s.offset( a, b ) = s._data.size();
            for ( size_t n = 0; n < causal.data.size(); n++ ) {
               if ( n == 0 || std::abs( causal.data[ n ] ) > thresholdToKeep ) {
                  s._data.emplace_back( causal.data[ n ] );
                  s._time.emplace_back( n == 0 ? 0 : n * causal.Ts - causal.tau );
                  std::cout << s._time.back() << " " << s._data.back()<< std::endl;
               }

            }
            s.length( a, b ) = s._data.size() - s.offset( a, b );
            std::cout << "Pruned " << ( ( 2 * freqs.size() - 2 ) - s.length( a, b ) ) << " DTIR entries out of " << ( 2 * freqs.size() - 2 ) << " less than " << thresholdToKeep << " (" << fracMaxToKeep * 100 << "% of max val)" << std::endl;

            port[ a ].s0[ b ] = s.data( a, b, 0 );
         }
         port[ a ].beta = 1.0 / ( 1 - s.data( a, a, 0 ) );
         port[ a ].R = port[ a ].beta * 50 * ( 1 + s.data( a, a, 0 ) );
      }
   }

   static void
   addToElements( const std::string & line, CircuitElements< T > & elements,
                  size_t & numNodes, size_t & numCurrents, size_t & numDCCurrents ) {
      // std::regex SParameterBlockInitialRegex( R"(^S(.*?)\s(\d+?)\s)" );
      std::regex SParameterBlockInitialRegex = generateRegex( "S", "w n s", true,
                                                              false );
      std::smatch matches;
      std::regex_search( line, matches, SParameterBlockInitialRegex );

      SParameterBlock< T > block;

      block.designator = "S";
      block.designator += matches.str( 1 );

      block.fracMaxToKeep = std::stod(matches.str( 2 ));

      size_t numPorts = std::stoul( matches.str( 3 ) );
      block.s.numPorts = numPorts;
      block.port = std::vector< SParameterPort< T > >( numPorts );

      auto strIter = matches.suffix().first;
      std::regex SParameterBlockPortRegex( R"(^(\d+?)\s(\d+?)\s)" );
      for ( size_t p = 0; p < numPorts; p++ ) {
         std::regex_search( strIter, line.cend(), matches,
                            SParameterBlockPortRegex );
         block.port[ p ].positive = std::stoi( matches.str( 1 ) );
         block.port[ p ].negative = std::stoi( matches.str( 2 ) );

         numNodes = std::max( numNodes, std::stoull( matches.str( 1 ) ) );
         numNodes = std::max( numNodes, std::stoull( matches.str( 2 ) ) );

         block.port[ p ].current = ++numCurrents;

         strIter = matches.suffix().first;
      }
      std::regex endOfLineRegex( R"(^(.*)$)" );
      std::regex_search( strIter, line.cend(), matches, endOfLineRegex );
      block.touchstoneFilePath = matches.str( 1 );
      block.readInTouchstoneFile();

      elements.dynamicElements.emplace_back(
          std::make_shared< SParameterBlock< T > >( block ) );
      for ( size_t p = 0; p < numPorts; p++ ) {
         elements.nodeComponentMap.insert(
             { { block.port[ p ].positive, elements.dynamicElements.back() },
               { block.port[ p ].negative, elements.dynamicElements.back() } } );
      }
   }
};

#endif
