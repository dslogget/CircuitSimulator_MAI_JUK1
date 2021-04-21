#include <iostream>
#include <chrono>
#include <vector>
#include <complex.h>
#include <cmath>
#include <numbers>
#include <iomanip>
#include "Maths/dft.hpp"

double
myFunctionToSample( double t ) {
   constexpr double freqSine = 1;
   return std::sin( 2 * std::numbers::pi * freqSine * t );
}


[[clang::optnone]] int
main() {
   // generate data
   constexpr double T = 2;
   constexpr double f_s = 8;
   constexpr size_t numPoints = T * f_s;


   std::vector< double > inputData; // = { 0, 1, 0, -1, 0, 1, 0, -1 };
   for ( int i = 0; i < numPoints; i++ ) {
      inputData.push_back( myFunctionToSample( static_cast< double >( i ) / f_s ) );
   }

   // print data as csv
   for ( const auto & num : inputData ) {
      std::cout << std::setprecision( 5 ) << num << ", ";
   }
   std::cout << std::endl;
   std::cout << std::endl;
   std::cout << std::endl;

   auto dftres = dft( inputData );
   auto fftres = fft( inputData );

   for ( const auto & num : dftres ) {
      std::cout << std::setprecision( 5 ) << std::abs( num ) << ", ";
   }
   std::cout << std::endl;
   std::cout << std::endl;
   std::cout << std::endl;

   for ( const auto & num : fftres ) {
      std::cout << std::setprecision( 5 ) << std::abs( num ) << ", ";
   }
   std::cout << std::endl;
   std::cout << std::endl;
   std::cout << std::endl;


   auto idftres = idft( dftres );
   for ( const auto & num : idftres ) {
      std::cout << std::setprecision( 5 ) << std::real( num ) << ", ";
   }
   std::cout << std::endl;
   std::cout << std::endl;
   std::cout << std::endl;

   auto ifftres = ifft( dftres );
   for ( const auto & num : ifftres ) {
      std::cout << std::setprecision( 5 ) << std::real( num ) << ", ";
   }
   std::cout << std::endl;
   std::cout << std::endl;
   std::cout << std::endl;
   return 0;
}
