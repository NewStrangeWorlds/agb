
#ifndef _PHYSICAL_CONSTANTS_H
#define	_PHYSICAL_CONSTANTS_H


namespace agb{  namespace constants{

//physical constants in cgs units from the GSL library
constexpr double planck_h = 6.62606896e-27;                     //Planck constant in g cm^2 / s
constexpr double light_c = 2.99792458e10;                       //speed of light in cm / s
constexpr double boltzmann_k = 1.3806504e-16;                   //Boltzmann constant in g cm^2 / K s^2
constexpr double stefan_boltzmann = 5.67040047374e-5;           //Stefan-Boltzmann constant in g / K^4 s^3
constexpr double mass_proton = 1.67262158e-24;                  //proton mass in g
constexpr double mass_electron = 9.10938188e-28;                //electron mass in g
constexpr double gravitation_const = 6.673e-8;                  //Gravitational constant in cm^3 / g s^2
constexpr double amu = 1.6605e-24;                              //atomic mass unit
constexpr double atomic_weight_c = 12.01;                       //atomic weight of carbon
constexpr double gas_constant = 8.314472e7;                     //universal gas constant erg K^-1 mol^-1

//additional constants
constexpr double radius_earth = 6.37123e8;                      //mean Earth radius in cm
constexpr double radius_sun = 696342.*100000.;                  //Solar radius in cm
constexpr double radius_jupiter = 7.1492e9;                     //mean Jupiter radius in cm
constexpr double mass_earth = 5.97219e24 * 1000.;               //Earth mass in g
constexpr double mass_jupiter = 1.89813e27 * 1000.;             //Jupiter mass in g
constexpr double mass_sun = 1.98847e30*1000;                    //Solar mass in g
constexpr double luminosity_sun = 3.828e26*1e7;                 //Solar luminosity in erg s-1


constexpr double parsec = 3.08567758135e18;                     //1 parsec in cm

constexpr double pi = 3.14159265358979323846;
constexpr double deg_to_rad = pi / 180.;
constexpr double rad_to_deg = 180./pi;


}}

#endif
