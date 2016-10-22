#ifndef RGEN_LIB_HEADER
#define RGEN_LIB_HEADER

#include <iostream>

#include "ext/phys_units/phys/units/quantity.hpp"
#include "ext/phys_units/phys/units/io.hpp"

using namespace phys::units;
using namespace phys::units::io;

#include "ext/cpplatex/latex.hpp"

#define NUM(x)      latex::math::make_num(x)
#define ADD(x,y)    latex::math::make_add(x,y)
#define SUB(x,y)    latex::math::make_sub(x,y)
#define MULT(x,y)   latex::math::make_mult(x,y)
#define DIV(x,y)    latex::math::make_fraction(x,y)
#define POW(x,y)    latex::math::make_exp(x,y)
#define ROOT(x,y)   latex::math::make_root(x,y)

const static auto GRAVITY = 9.81 * (meter / (second * second));
const static auto UNIVERSAL_GAS_CONSTANT = (8.31446 * joule) / (mole * kelvin);

struct Fuel {
    std::string name;
    double molar_mass;
};

struct Oxidizer {
    std::string name;
    double molar_mass;
};

struct Propellants {
    Fuel     fuel;
    Oxidizer oxidizer;
    double   mixture_ratio;
};

class EngineBasis {
public:
    std::string                name;
    std::string                version;

    quantity<force_d>          thrust;
    quantity<time_interval_d>  isp;
    quantity<pressure_d>       pressure;
    Propellants                propellants;

    EngineBasis(
        const std::string& name, const std::string& version,
        double thrust, double isp, double pressure,
        Propellants prop
    )
        : name(name), version(version)
        , thrust(thrust * newton), isp(isp * second), pressure(pressure * pow(10, 6) /*to MPa*/ * pascal)
        , propellants(prop)
    {}
};

#endif
