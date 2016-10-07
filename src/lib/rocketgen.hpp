#ifndef RGEN_LIB_HEADER
#define RGEN_LIB_HEADER

#include <iostream>

#include "../ext/phys_units/phys/units/io.hpp"
#include "../ext/phys_units/phys/units/quantity.hpp"

using namespace phys::units;
using namespace phys::units::io;

const static auto GRAVITY = 9.81 * (meter / (second * second));
const static auto UNIVERSAL_GAS_CONSTANT = (8.31446 * joule) / (mole * kelvin);

#endif
