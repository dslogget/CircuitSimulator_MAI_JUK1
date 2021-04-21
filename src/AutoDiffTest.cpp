#include <iostream>
#include "Maths/AutoDifferentiation.hpp"
#include <assert.h>
#include <iomanip>
#include <chrono>

namespace AD = AutoDifferentiation;

void testBasicOutput_NoCheck() {
      constexpr AD::DiffVar< double, 2 > x( 1, 1, 0 );
      constexpr AD::DiffVar< double, 2 > y( 2, 0, 1 );

      constexpr auto f = x*y;
      auto sin = AD::sin( x );
      auto sin2 = AD::sin( 2*x );
      auto tanh = AD::tanh( 1 + 2*x );

      std::cout << "f\n" << f << std::endl;
      std::cout << "sin\n" << sin << std::endl;
      std::cout << "sin2\n" << sin2 << std::endl;
      std::cout << "tanh2\n" << tanh << std::endl;
}

template< typename T >
struct BJTResults {
   T g_ee;
   T g_ec;
   T g_ce;
   T g_cc;

   T I_e;
   T I_c;
};

template< typename T >
void testBJTModelControl( const T & base_v_be, const T & base_v_bc, BJTResults< T > & results ) {
   // The control

   constexpr T alpha_f = 0.99;
   constexpr T alpha_r = 0.02;

   constexpr T I_es = 2e-14;
   constexpr T V_Te = 26e-3;
   constexpr T I_cs = 99e-14;
   constexpr T V_Tc = 26e-3;

   T control_v_be = base_v_be;
   T control_v_bc = base_v_bc;

   T control_i_e = -I_es * ( std::exp( control_v_be / V_Te ) - 1 ) +
           alpha_r * I_cs * ( std::exp( control_v_bc / V_Tc ) - 1 );
   T control_i_c = alpha_f * I_es * ( std::exp( control_v_be / V_Te ) - 1 ) -
           I_cs * ( std::exp( control_v_bc / V_Tc ) - 1 );

   T control_g_ee = ( I_es / V_Te ) * std::exp( control_v_be / V_Te );
   T control_g_ec = alpha_r * ( I_cs / V_Tc ) * std::exp( control_v_bc / V_Tc );
   T control_g_ce = alpha_f * ( I_es / V_Te ) * std::exp( control_v_be / V_Te );
   T control_g_cc = ( I_cs / V_Tc ) * std::exp( control_v_bc / V_Tc );

   T control_I_e = control_i_e + control_g_ee * control_v_be - control_g_ec * control_v_bc;
   T control_I_c = control_i_c - control_g_ce * control_v_be + control_g_cc * control_v_bc;

   results.g_ee = control_g_ee;
   results.g_ec = control_g_ec;
   results.g_ce = control_g_ce;
   results.g_cc = control_g_cc;

   results.I_e = control_I_e;
   results.I_c = control_I_c;
}

template< typename T >
void testBJTModelAutoDiff( const T & base_v_be, const T & base_v_bc, BJTResults< T > & results ) {
   // AutoDifferentiation
   constexpr T alpha_f = 0.99;
   constexpr T alpha_r = 0.02;

   constexpr T I_es = 2e-14;
   constexpr T V_Te = 26e-3;
   constexpr T I_cs = 99e-14;
   constexpr T V_Tc = 26e-3;
   using ADT = AD::DiffVar< T, 2 >;

   ADT ad_v_be = ADT( base_v_be, 1, 0 );
   ADT ad_v_bc = ADT( base_v_bc, 0, 1 );

   ADT ad_i_e = -I_es * ( AD::exp( ad_v_be / V_Te ) - 1 ) +
           alpha_r * I_cs * ( AD::exp( ad_v_bc / V_Tc ) - 1 );
   ADT ad_i_c = alpha_f * I_es * ( AD::exp( ad_v_be / V_Te ) - 1 ) -
           I_cs * ( AD::exp( ad_v_bc / V_Tc ) - 1 );

   // Najm defines the g_ee and g_cc to be the negated partial
   T ad_g_ee = -ad_i_e[ 1 ];
   T ad_g_ec = ad_i_e[ 2 ];
   T ad_g_ce = ad_i_c[ 1 ];
   T ad_g_cc = -ad_i_c[ 2 ];

   T ad_I_e = ad_i_e[ 0 ] + ad_g_ee * ad_v_be[ 0 ] - ad_g_ec * ad_v_bc[ 0 ];
   T ad_I_c = ad_i_c[ 0 ] - ad_g_ce * ad_v_be[ 0 ] + ad_g_cc * ad_v_bc[ 0 ];

   results.g_ee = ad_g_ee;
   results.g_ec = ad_g_ec;
   results.g_ce = ad_g_ce;
   results.g_cc = ad_g_cc;

   results.I_e = ad_I_e;
   results.I_c = ad_I_c;
}

template< typename T >
[[clang::optnone]]
void testBJTModel() {
   // Model constants
   constexpr T I_es = 2e-14;
   constexpr T V_Te = 26e-3;
   constexpr T I_cs = 99e-14;
   constexpr T V_Tc = 26e-3;

   BJTResults< T > controlResults;
   BJTResults< T > autoDiffResults;

   T V_bc_crit = V_Tc * std::log( V_Tc / ( I_cs * std::sqrt( 2 ) ) );
   T V_be_crit = V_Te * std::log( V_Te / ( I_es * std::sqrt( 2 ) ) );

   std::array< T, 6 > base_v_be_vec = { 0, 1, 2, 3, 4, 5 };
   std::array< T, 6 > base_v_bc_vec = { 0, 1, 2, 3, 4, 5 };;

   size_t controlAccumulate = 0;
   size_t autoDiffAccumulate = 0;

   for ( auto base_v_be : base_v_be_vec ) {
      for ( auto base_v_bc : base_v_bc_vec ) {
         std::cout << "base_v_be=" << std::setw( 10 ) << base_v_be
                   << ", base_v_bc=" << std::setw( 10 ) << base_v_bc
                   << std::endl;

         base_v_be = std::min( V_be_crit, base_v_be );
         base_v_bc = std::min( V_bc_crit, base_v_bc );

         auto startTime = std::chrono::high_resolution_clock::now();
         for ( size_t i = 0; i < 10000; i++ ) {
            testBJTModelControl( base_v_be, base_v_bc, controlResults );
         }
         auto endTime = std::chrono::high_resolution_clock::now();
         auto timeTaken = std::chrono::duration_cast< std::chrono::microseconds >(
                        endTime - startTime )
                        .count();
         controlAccumulate += timeTaken;

         std::cout << std::setw( 15 ) << timeTaken << " us | ";

         startTime = std::chrono::high_resolution_clock::now();
         for ( size_t i = 0; i < 10000; i++ ) {
            testBJTModelAutoDiff( base_v_be, base_v_bc, autoDiffResults );
         }
         endTime = std::chrono::high_resolution_clock::now();
         timeTaken = std::chrono::duration_cast< std::chrono::microseconds >(
                        endTime - startTime )
                        .count();
         autoDiffAccumulate += timeTaken;

         std::cout << std::setw( 15 ) << timeTaken << " us" << std::endl;

         // Comparisons

         std::cout << "g_ee | control=" << std::setw( 15 ) << controlResults.g_ee
                       << " | autoDiff=" << std::setw( 15 ) << autoDiffResults.g_ee << std::endl;
         std::cout << "g_ec | control=" << std::setw( 15 ) << controlResults.g_ec
                       << " | autoDiff=" << std::setw( 15 ) << autoDiffResults.g_ec << std::endl;
         std::cout << "g_ce | control=" << std::setw( 15 ) << controlResults.g_ce
                       << " | autoDiff=" << std::setw( 15 ) << autoDiffResults.g_ce << std::endl;
         std::cout << "g_cc | control=" << std::setw( 15 ) << controlResults.g_cc
                       << " | autoDiff=" << std::setw( 15 ) << autoDiffResults.g_cc << std::endl;

         std::cout << "I_e  | control=" << std::setw( 15 ) << controlResults.I_e
                       << " | autoDiff=" << std::setw( 15 ) << autoDiffResults.I_e << std::endl;
         std::cout << "I_c  | control=" << std::setw( 15 ) << controlResults.I_c
                       << " | autoDiff=" << std::setw( 15 ) << autoDiffResults.I_c << std::endl;

         assert( std::abs( controlResults.I_e - autoDiffResults.I_e ) < 1e-12 );
         assert( std::abs( controlResults.I_c - autoDiffResults.I_c ) < 1e-12 );
      }
   }

   std::cout << std::setw( 15 ) << controlAccumulate << " us | ";
   std::cout << std::setw( 15 ) << autoDiffAccumulate << " us" << std::endl;

}

template< typename T >
struct TransistorTestResult {
   T var;
   T diff1;
   T diff2;
};

template< typename T >
void transistorTestControl( T V_gs_in, T V_ds_in, TransistorTestResult< T > & toRet ) {
   constexpr T alpha = 1.3;
   constexpr T beta0 = 0.42;
   constexpr T gamma = 0.0005;
   constexpr T delta = 0.3;
   constexpr T xi = 0.06;
   constexpr T lambda = 1.5;
   constexpr T mu = 0.0;
   constexpr T zeta = 0.18;
   constexpr T Vto = -2.4;

   T V_gs = V_gs_in;
   T V_ds = V_ds_in;

   auto beta = beta0;
   auto Vgst = V_gs - ( 1 + beta*beta ) * Vto + gamma * V_ds;
   auto Veff = 0.5 * ( Vgst + std::pow( std::pow( Vgst, 2 ) + delta*delta, 0.5 ) );
   auto power = lambda / ( 1 + mu * std::pow( V_ds, 2 ) + xi * Veff );
   auto area = alpha * V_ds * ( 1 + zeta * Veff );
   auto f1 = std::tanh( area );
   auto Ids_lim = beta * std::pow( Veff, power );
   auto Idrain = Ids_lim * f1;
   toRet.var = Idrain;

   auto dVeff_dVgs = 0.5*( 1 + Vgst * std::pow( Vgst * Vgst + delta * delta, -0.5 ) );
   auto dpower_dVgs = -lambda * xi * dVeff_dVgs * std::pow( power / lambda, 2 );
   auto df_dVgs = std::pow(1 / cosh(area), 2) * alpha * V_ds * zeta * dVeff_dVgs;
   toRet.diff1 = Idrain * ( power * ( dVeff_dVgs / Veff ) + std::log( Veff ) * dpower_dVgs ) + Ids_lim * df_dVgs ;

   auto dVeff_dVds = 0.5 * ( gamma + std::pow( Vgst * Vgst + delta * delta, -0.5 ) * Vgst *gamma );
   auto dpower_dVds = -lambda * ( 2 * mu * V_ds + xi * dVeff_dVds ) * std::pow( power / lambda, 2 );
   auto df_dVds = std::pow(1 / cosh(area), 2) * alpha * ( 1 + zeta * ( V_ds * dVeff_dVds + Veff ) );
   toRet.diff2 = Idrain *( power * ( dVeff_dVds / Veff ) + std::log( Veff ) * dpower_dVds ) + Ids_lim * df_dVds;
}

template< typename T >
void transistorTestAudoDiff( T V_gs_in, T V_ds_in, AD::DiffVar< T, 2 > & toRet ) {
   constexpr T alpha = 1.3;
   constexpr T beta0 = 0.42;
   constexpr T gamma = 0.0005;
   constexpr T delta = 0.3;
   constexpr T xi = 0.06;
   constexpr T lambda = 1.5;
   constexpr T mu = 0.0;
   constexpr T zeta = 0.18;
   constexpr T Vto = -2.4;

   using ADT = AD::DiffVar< T, 2 >;
   ADT V_gs( V_gs_in, 1, 0 );
   ADT V_ds( V_ds_in, 0, 1 );

   auto beta = beta0;
   auto Vgst = V_gs - ( 1 + beta*beta ) * Vto + gamma * V_ds;
   auto Veff = 0.5 * ( Vgst + AD::sqrt( AD::pow( Vgst, 2 ) + delta*delta ) );
   auto power = lambda / ( 1 + mu * AD::pow( V_ds, 2 ) + xi * Veff );
   auto area = alpha * V_ds * ( 1 + zeta * Veff );
   auto f1 = AD::tanh( area );
   auto Ids_lim = beta * AD::pow( Veff, power );
   toRet = Ids_lim * f1;
}

template< typename T >
[[clang::optnone]]
void transistorTest() {
   using ADT = AD::DiffVar< T, 2 >;
   std::vector< T > V_gs_in_vec = { 0, 1, 2, 3 };
   std::vector< T > V_ds_in_vec = { 0, 1, 2, 3 };
   size_t controlAccumulate = 0;
   size_t autoDiffAccumulate = 0;
   for ( auto V_gs_in : V_gs_in_vec ) {
      for ( auto V_ds_in : V_gs_in_vec ) {
         ADT autoDiffResult( 0 );
         TransistorTestResult< T > res;
         auto startTime = std::chrono::high_resolution_clock::now();
         for ( size_t i = 0; i < 10000; i++ ) {
            transistorTestAudoDiff( V_gs_in, V_ds_in, autoDiffResult );
         }
         auto endTime = std::chrono::high_resolution_clock::now();
         auto timeTaken = std::chrono::duration_cast< std::chrono::microseconds >(
                        endTime - startTime )
                        .count();
         autoDiffAccumulate += timeTaken;
         std::cout << std::setw( 15 ) << timeTaken << " us | ";

         startTime = std::chrono::high_resolution_clock::now();
         for ( size_t i = 0; i < 10000; i++ ) {
            transistorTestControl( V_gs_in, V_ds_in, res );
         }
         endTime = std::chrono::high_resolution_clock::now();
         timeTaken = std::chrono::duration_cast< std::chrono::microseconds >(
                        endTime - startTime )
                        .count();
         controlAccumulate += timeTaken;
         std::cout << std::setw( 15 ) <<  timeTaken << " us" << std::endl;

         std::cout << "V_gs: " << V_gs_in << " V_ds: " << V_ds_in << " " << std::endl;
         std::cout << std::setw( 15 ) << res.var << " ";
         std::cout << std::setw( 15 ) << res.diff1 << " ";
         std::cout << std::setw( 15 ) << res.diff2 << std::endl;
         std::cout << std::setw( 15 ) << autoDiffResult[ 0 ] << " ";
         std::cout << std::setw( 15 ) << autoDiffResult[ 1 ] << " ";
         std::cout << std::setw( 15 ) << autoDiffResult[ 2 ] << std::endl << std::endl;
         assert( std::abs( res.var - autoDiffResult[ 0 ] ) < 1e-12 );
         assert( std::abs( res.diff1 - autoDiffResult[ 1 ] ) < 1e-12 );
         assert( std::abs( res.diff2 - autoDiffResult[ 2 ] ) < 1e-12 );
      }
   }
   std::cout << std::setw( 15 ) << autoDiffAccumulate << " us | ";
   std::cout << std::setw( 15 ) << controlAccumulate << " us" << std::endl;

}

int main(int argc, char *argv[])
{
   testBasicOutput_NoCheck();

   testBJTModel< double >();

   transistorTest< double >();

   return 0;
}
