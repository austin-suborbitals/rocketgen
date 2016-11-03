#ifndef RGEN_LIB_HEADER
#define RGEN_LIB_HEADER

#include <iostream>

#include "ext/phys_units/phys/units/io.hpp"
#include "ext/phys_units/phys/units/quantity.hpp"

using namespace phys::units;
using namespace phys::units::io;

template <typename P, typename D, typename X, typename std::enable_if<std::is_floating_point<P>::value, char>::type = 0>
constexpr auto pow(const quantity<D,X>& val, P pwr) {
    if (pwr != 0.5) { // TODO: sigh...
        throw std::runtime_error("The dependency phys_units does not allow the runtime flexibility needed. This is a TODO.");
    }
    return sqrt(val);
}


auto mm = meter * milli; // TODO: explain why mult not div
auto mm_sq = mm*mm;
auto mm_cu = mm*mm*mm;

auto cm = meter * centi;
auto cm_cu = cm*cm*cm;

auto MPa = pascal * 1'000'000;




#include "ext/cpplatex/latex.hpp"

using namespace latex::doc;
using namespace latex::math;

constexpr static const char* vspace = "\\vspace{10 mm}";
constexpr static const char* svspace = "\\vspace{6 mm}";


// styling
using VarText = latex::Text<latex::math::style::Bold, latex::math::style::Italic>;
using ConstText = latex::Text<latex::math::style::Italic>;
using StyledVar = latex::math::Variable<VarText>;

template <typename T>
class StyledValVar: public latex::math::ValuedVariable<T, VarText> {
public:
    StyledValVar(const T& t, const std::string& str) : latex::math::ValuedVariable<T, VarText>(t, str) {}
    StyledValVar(const T& t, const char* str) : latex::math::ValuedVariable<T, VarText>(t, str) {}
};

template <typename T>
class StyledConst: public latex::math::ValuedVariable<T, ConstText> {
public:
    StyledConst(const T& t, const std::string& str) : latex::math::ValuedVariable<T, ConstText>(t, str) {}
    StyledConst(const T& t, const char* str) : latex::math::ValuedVariable<T, ConstText>(t, str) {}
};


#include "util.hpp"

#define NUM(x)      make_num(x)
#define ADD(x,y)    make_add(x,y)
#define SUB(x,y)    make_sub(x,y)
#define MULT(x,y)   make_mult(x,y)
#define DIV(x,y)    make_fraction(x,y)
#define POW(x,y)    make_pow(x,y)
#define EXP(x)      make_exp(x)
#define LOG(x,y)    make_log(x,y)
#define LOGN(x)     make_ln(x)
#define ROOT(x,y)   make_root(x,y)
#define PAREN(x)    make_paren(x)
#define SIN(x)      make_sin(x)
#define COS(x)      make_cos(x)
#define TAN(x)      make_tan(x)

#define TORAD(x)    (x*(M_PI/180))
#define TODEG(x)    (x*(180/M_PI))

const static auto GRAVITY = 9.81 * (meter / (second * second));
const static auto GAS_CONSTANT = 8.31446 * (((meter * meter * meter) * pascal) / (mole * kelvin));

template <typename T, typename C, typename N>
std::string scale_to_string(const T& val, const C& cutoff, const N& new_scale, const char* unit_str) {
    return magnitude(val) > cutoff
        ? to_string(val)
        : strjoin(val.to(new_scale) /* returns a num */, unit_str);
}


template <typename T, typename I>
auto calc_mass_flow(const T& thrust, const I& isp) {
    return thrust / (GRAVITY * isp);
}


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
    double   gamma;     // TODO: feasible to calculate this?
    quantity<mass_d>   avg_combustion_weight; // TODO: molar mass
    quantity<thermodynamic_temperature_d> flame_temperature; // TODO: feasible to calculate this?
};

struct EngineBasis {
    std::string                name;
    std::string                version;

    quantity<force_d>          thrust;
    quantity<time_interval_d>  isp;
    quantity<pressure_d>       pressure;
    quantity<pressure_d>       external_pressure;
    double                     expansion_ratio;
    quantity<length_d>         l_star;
    double                     converging_angle;
    double                     diverging_angle;
    Propellants                propellants;
};



//
// throat
//

struct Throat {
    quantity<pressure_d> pressure;
    quantity<thermodynamic_temperature_d> temperature;
    quantity<length_d> radius;

    auto area() const { return (nth_power<2>(radius) * M_PI); }
    auto diameter() const  { return radius * 2; }
};

template <typename Gamma, typename Temp>
auto throat_temperature_eqn(const Gamma& gamma, const Temp& flame_temp) {
    auto temperature_coefficient = 
        NUM(1)
        /
        (NUM(1) + ((StyledValVar<Gamma>(gamma, "\\gamma") - NUM(1))/NUM(2)));
    return MULT(flame_temp, temperature_coefficient);
}

template <typename Gamma, typename Pressure>
auto throat_pressure_eqn(const Gamma& gamma, const Pressure& press) {
    auto gamma_var = StyledValVar<Gamma>(gamma, "\\gamma");
    auto temperature_coefficient = 
        (1 + DIV(gamma_var - 1, 2)).pow(-1 * DIV(gamma_var, gamma_var - 1));
    return MULT(temperature_coefficient, press);
}

template <typename T, typename P, typename O>
auto build_throat_area_section(O& os, const EngineBasis& rocket, const T& temp, const P& pressure) {
    auto mass_flow = calc_mass_flow(rocket.thrust, rocket.isp);

    StyledValVar<decltype(mass_flow)> mflow_var(mass_flow, "M_{flow}");

    auto gas_var = StyledValVar<decltype(GAS_CONSTANT)>(GAS_CONSTANT, "R");
    auto r_eqn = gas_var / NUM(rocket.propellants.avg_combustion_weight / mole);
    StyledValVar<decltype(r_eqn)> r_var(r_eqn, "R'");

    os << vspace
       << strjoin("First, we find $", r_var, "$ using the following equation: ", vspace)
       << make_aligned_eqn(
            r_var,
            r_eqn,
            r_eqn.solve()
       )
       << svspace
       << "where: "
       << svspace
       << make_eqn(gas_var, gas_var.solve())
       << vspace
       << strjoin("This allows us to find the area of the throat cross-section via:", vspace);


    StyledValVar<double> gamma_var(rocket.propellants.gamma, "\\gamma");

    auto flow_half = (NUM(mass_flow)/pressure);
    auto temp_part = r_var * NUM(temp);
    auto sqrt_half = (temp_part / gamma_var).sqrt();
    auto area_eqn = flow_half * sqrt_half;
    auto area = area_eqn.solve();

    os << make_aligned_eqn(
        StyledVar("A_{throat}"),
        MULT(flow_half, sqrt_half),
        MULT(make_paren(flow_half.solve()), make_paren(sqrt_half.solve())),
        scale_to_string(area, 1, mm_sq, "mm+2")
    )
    << vspace;


    auto pi_var = StyledConst<double>(M_PI, "\\pi");
    auto radius_eqn = (NUM(area) / pi_var).sqrt();
    auto radius = radius_eqn.solve();
    auto diameter_eqn = StyledValVar<decltype(radius)>(radius, "r_{throat}") * NUM(2);
    auto diameter = diameter_eqn.solve();
    os << "From here, basic geometry gets us the diameter/radius:\n\n"
       << vspace
       << make_aligned_eqn(
        StyledVar("r_{throat}"),
        radius_eqn,
        scale_to_string(radius, 1, mm, "mm")
       )
       << vspace
       << make_aligned_eqn(
        StyledVar("d_{throat}"),
        diameter_eqn,
        scale_to_string(diameter, 1, mm, "mm")
       );

    return radius;
}


Throat build_throat(latex::doc::Report& doc, const EngineBasis& rocket) {
    using namespace latex;
    using PrefixText = Text<latex::style::Bold>;


    // open the section
    latex::doc::Section section("Throat", true);
    section << "Defining characteristics and calculations for the throat joining the nozzle and combustion chamber.\n";

    // list our critical constants
    Subsection critical_vars("Critical Constants");
    critical_vars << "Underlying constants that form the throat profile:";
    critical_vars << (UnorderedList()
                        << strjoin(PrefixText("Gamma (prop heat capacity ratio): "), rocket.propellants.gamma)
                        << strjoin(PrefixText("Propellant Combustion Temp (Kelvin): "), rocket.propellants.flame_temperature)
                        << strjoin(PrefixText("Chamber Pressure: "), eng::to_string(rocket.pressure))
                     );
    section << critical_vars;

    //
    // temperature
    //

    Subsection temp_sub("Temperature Profile");
    auto tmp_eqn = throat_temperature_eqn(rocket.propellants.gamma, rocket.propellants.flame_temperature);
    quantity<thermodynamic_temperature_d> temperature = tmp_eqn.solve();
    auto temp_throat_eqn = make_aligned_eqn(
        StyledVar("T_{throat}"),
        StyledVar("T_{coefficient}") * StyledVar("T_{chamber}"),
        tmp_eqn,
        temperature
    );
    temp_sub
        << "The following equation gives the temperature to be expected at the throat plane of the nozzle."
        << temp_throat_eqn.latex();

    //
    // pressure
    //

    Subsection pressure_sub("Pressure Profile");
    auto pressure_eqn = throat_pressure_eqn(rocket.propellants.gamma, rocket.pressure);
    quantity<pressure_d> pressure = pressure_eqn.solve();
    auto pressure_throat_eqn = make_aligned_eqn(
        StyledVar("P_{throat}"),
        StyledVar("P_{coefficient}") * StyledVar("P_{chamber}"),
        pressure_eqn,
        pressure
    );
    pressure_sub
        << "The following equation gives the pressure to be expected at the throat plane of the nozzle."
        << pressure_throat_eqn.latex();


    //
    // geometry
    //

    Subsection geometry("Throat Geometry");
    geometry << "The following equations define the geometrical/physical dimensions of the throat.\n\n\n";
    auto radius = build_throat_area_section(geometry, rocket, temperature, pressure);
    Throat throat_desc{pressure, temperature, radius};
    geometry << "\n\n" << vspace << "\n\n";


    section << temp_sub;
    section << pressure_sub;
    section << geometry
            << (Subsection("Summary")
                << (UnorderedList() 
                    << strjoin(PrefixText("Throat Pressure: "),     eng::to_string(throat_desc.pressure))
                    << strjoin(PrefixText("Throat Temperature: "),  throat_desc.temperature)
                    << strjoin(PrefixText("Throat Area: "),         scale_to_string(throat_desc.area(), 1, mm_sq, "mm+2"))
                    << strjoin(PrefixText("Throat Radius: "),       scale_to_string(throat_desc.radius, 1, mm, "mm"))
                )
            );

    doc << section;
    return throat_desc;
}




struct Nozzle {
    quantity<length_d> throat_radius;
    quantity<length_d> exit_radius;

    quantity<length_d> diverging_length;

    auto exit_area() const { return (nth_power<2>(exit_radius) * M_PI); }
    auto exit_diameter() const  { return exit_radius * 2; }

    auto throat_area() const { return (nth_power<2>(throat_radius) * M_PI); }
    auto throat_diameter() const  { return throat_radius * 2; }
};

Nozzle build_nozzle(latex::doc::Report& doc, const EngineBasis& rocket, const Throat& throat) {
    using PrefixText = latex::Text<latex::style::Bold>;

    // open the section
    latex::doc::Section section("Nozzle", true);
    section << "Defining characteristics and calculations for the nozzle.\n";


    // list our critical constants
    Subsection critical_vars("Critical Constants");
    critical_vars << "Underlying constants that form the nozzle profile:";
    critical_vars << (UnorderedList()
                        << strjoin(PrefixText("Throat Area: "), scale_to_string(throat.area(), 1, mm_sq, "mm+2"))
                        << strjoin(PrefixText("Expansion Ratio: "), rocket.expansion_ratio)
                        // TODO: external pressure?
                     );

    //
    // geometry
    //

    Subsection geometry("Nozzle Geometry");
    geometry << "The following equations define the geometrical/physical dimensions of the nozzle and exit area.\n\n\n";

    StyledValVar<decltype(throat.area())> area_throat(throat.area(), "A_{throat}");

    auto exit_area_eqn = area_throat * NUM(rocket.expansion_ratio);
    auto exit_area = exit_area_eqn.solve();

    auto pi_var = StyledConst<double>(M_PI, "\\pi");
    auto radius_eqn = (NUM(exit_area) / pi_var).sqrt();
    auto radius = radius_eqn.solve();
    auto diameter_eqn = StyledValVar<decltype(radius)>(radius, "r") * NUM(2);
    auto diameter = diameter_eqn.solve();

    // the properties of the diverging angle
    StyledConst<double> ninty_deg_as_rad(TORAD(90), "90\\degree");
    StyledConst<decltype(rocket.diverging_angle)> angle_var(rocket.diverging_angle, "\\theta");
    auto diverging_height_eqn = StyledValVar<decltype(radius)>(radius, "R_{nozzle}") - StyledValVar<decltype(throat.radius)>(throat.radius, "R_{throat}");
    auto diverging_height = diverging_height_eqn.solve();
    StyledValVar<decltype(diverging_height)> height_var(diverging_height, "R_{diff}");
    auto diverging_length_eqn = (height_var / SIN(angle_var)) * SIN(ninty_deg_as_rad-angle_var);
    auto diverging_length = diverging_length_eqn.solve();

    geometry
        << vspace
        << make_aligned_eqn(
            StyledVar("A_{exit}"),
            exit_area_eqn,
            scale_to_string(exit_area, 1, mm_sq, "mm+2")
        )
        << vspace
        << "Just as in the throat calculaions, basic geometry gives us radius and diameter.\n\n\n"
        << svspace
        << make_aligned_eqn(
            StyledVar("r_{exit}"),
            radius_eqn,
            scale_to_string(radius, 1, mm, "mm")
        )
        << vspace
        << make_aligned_eqn(
            StyledVar("d_{exit}"),
            diameter_eqn,
            scale_to_string(diameter, 1, mm, "mm")
        )
        << svspace
        << "Via trigonometry we find the nozzle length:\n\n"
        << svspace
        << make_aligned_eqn(
            height_var,
            diverging_height_eqn,
            scale_to_string(diverging_height, 1, mm, "mm")
        )
        << svspace
        << make_aligned_eqn(
            StyledVar("L_{nozzle}"),
            diverging_length_eqn,
            scale_to_string(diverging_length, 1, mm, "mm")
        )
        << svspace
        << strjoin("For the given diverging angle of:", TODEG(rocket.diverging_angle), "\\degree\n\n")
    ;


    section
        << critical_vars
        << geometry;

    doc << section;
    return Nozzle{throat.radius, radius, diverging_length};
}








struct Chamber {
    quantity<length_d> flatwall_radius;
    quantity<length_d> flatwall_length;

    quantity<length_d> converging_length;
    quantity<length_d> converging_height;

    const quantity<length_d> throat_radius;

    quantity<length_d> diameter() const { return flatwall_radius * 2; }
    quantity<length_d> total_length() const { return flatwall_length + converging_length; }

    // TODO: express m^2 in return type
    auto flatwall_area() const { return nth_power<2>(flatwall_radius) * M_PI; }

    // TODO: express m^3 in return type
    auto flatwall_volume() const { return flatwall_area() * flatwall_length; }

    // TODO: express m^3 in return type
    auto converging_volume() const {
        return converging_length /* h */ * (M_PI / 3) * (nth_power<2>(flatwall_radius) + (flatwall_radius * throat_radius) + nth_power<2>(throat_radius));
    }

    // TODO: express m^3 in return type
    auto volume() const { return converging_volume() + flatwall_volume(); }
};

Chamber build_chamber(latex::doc::Report& doc, const EngineBasis& rocket, const Throat& throat) {
    using PrefixText = latex::Text<latex::style::Bold>;

    // open the section
    latex::doc::Section section("Combustion Chamber", true);
    section << "Defining characteristics and calculations for the combustion chamber.\n";

    // list our critical constants
    Subsection critical_vars("Critical Constants");
    critical_vars << "Underlying constants that form the combustion chamber profile:";
    critical_vars << (UnorderedList()
                        << strjoin(PrefixText("Throat Area: "), scale_to_string(throat.area(), 1, mm_sq, "mm+2"))
                        << strjoin(PrefixText("Throat Dimeter: "), scale_to_string(throat.diameter(), 1, mm, "mm"))
                        << strjoin(PrefixText("Converging Angle ($\\theta$): "), TODEG(rocket.converging_angle), "\\degree")
                        << strjoin(PrefixText("L*: "), scale_to_string(rocket.l_star, 1, mm, "mm"))
                     );

    //
    // geometry
    //

    Subsection geometry("Combustion Chamber Geometry");
    geometry << "The following equations define the geometrical/physical dimensions of the chamber and converging portion of the nozzle.\n\n\n";

    // make valued variables for the throat and L* variables
    StyledValVar<decltype(throat.area())> area_throat(throat.area(), "A_{throat}");
    StyledValVar<decltype(throat.diameter())> diameter_throat(throat.diameter(), "D_{throat}");
    StyledValVar<quantity<length_d>> l_star_var(rocket.l_star, "L*");

    // make any other known-constants
    StyledConst<double> pi_var(M_PI, "\\pi");
    StyledConst<double> angle_var(rocket.converging_angle, "\\theta");

    // solve for volume
    auto volume_eqn = area_throat * l_star_var;
    auto volume = volume_eqn.solve();
    StyledValVar<decltype(volume)> volume_var(volume, "V_{chamber}");

    // make an estimate for the chamber length
    auto initial_len_est_eqn = EXP(NUM(0.029)*LOGN(diameter_throat).pow(2) + NUM(0.047)*LOGN(diameter_throat) + NUM(1.94));
    auto initial_len_est = (initial_len_est_eqn.solve() * (meter / 100)); // TODO: why centimeters?
    StyledValVar<decltype(initial_len_est)> initial_length_var(initial_len_est, "L_{estimate}");

    // get a diameter from this estimate
    auto initial_radius = sqrt((volume/initial_len_est)/M_PI);
    StyledValVar<decltype(initial_len_est)> est_diameter_var(2*initial_radius, "D_{estimate}");

    // the equation we iterate on to approach a better D_chamber
    auto iter_eqn = [&](double iter_val){
        return
        ((diameter_throat.pow(3) + ((NUM(24)/pi_var) * TAN(angle_var) * volume_var))
        /
        (iter_val + (NUM(6) * TAN(angle_var) * initial_length_var))) // TODO: always use initial? or Lc=Vc/calc(Dc)?
        .sqrt();
    };

    // do the iteration
    std::size_t num_iters = 100;
    auto initial_iter_eqn = iter_eqn(est_diameter_var.solve()); // TODO: explain the static value
    auto diam_iter = iter_eqn(est_diameter_var.solve()).solve();
    for (std::size_t i = 1; i < num_iters; ++i) {
        diam_iter = iter_eqn(diam_iter).solve();
    }

    // solve for the actual diameter and radius
    auto diameter = diam_iter*meter;
    auto radius = diameter / 2;

    // now calculate the area of the (non-converging section) chamber
    auto area_eqn = NUM(nth_power<2>(radius)) * pi_var;
    auto area = area_eqn.solve();
    auto area_var = StyledValVar<decltype(area)>(area, "A_{chamber}");

    // we can also now calc the contraction ratio
    auto contraction_eqn = area_var / area_throat;
    auto contraction_ratio = contraction_eqn.solve();

    // as well as the chamber length
    auto length = volume / area;

    // ... and the properties of the converging angle
    auto converging_height_eqn = StyledValVar<decltype(radius)>(radius, "R_{conv}") - StyledValVar<decltype(throat.radius)>(throat.radius, "R_{throat}");
    auto converging_height = converging_height_eqn.solve();
    StyledValVar<decltype(converging_height)> height_var(converging_height, "R_{diff}");

    StyledConst<double> ninty_deg_as_rad(TORAD(90), "90\\degree");
    auto converging_length_eqn = (height_var / SIN(angle_var)) * SIN(ninty_deg_as_rad-angle_var);
    auto converging_length = converging_length_eqn.solve();

    // and misc
    auto flat_eqn = StyledValVar<decltype(length)>(length, "L_{chamber}") - StyledValVar<decltype(converging_length)>(converging_length, "L_{conv}");


    geometry
        << vspace
        << make_aligned_eqn(
            StyledVar("V_{chamber}"),
            volume_eqn,
            MULT(PAREN(area_throat.solve()), PAREN(l_star_var.solve())),
            scale_to_string(volume, 1, mm_cu, "mm+3")
        )
        << vspace
        << "We achieve an initial approximation of the chamber length from the following formula:\n\n"
        << svspace
        << make_aligned_eqn(
            StyledVar("L_{estimate}"),
            initial_len_est_eqn,
            scale_to_string(initial_len_est, 1, mm, "mm")
        )
        << svspace
        << make_aligned_eqn(
            StyledVar("D_{estimate}"),
            scale_to_string(initial_radius*2, 1, mm, "mm")
        )
        << "Which we can further refine by solving the following via iteration:\n\n"
        << vspace
        << make_eqn(est_diameter_var, initial_iter_eqn)
        << vspace
        << strjoin("Which yields (after", num_iters, "iterations):\n\n")
        << svspace
        << make_eqn(StyledVar("D_{chamber}"), scale_to_string(diameter, 1, mm, "mm"))
        << make_eqn(StyledVar("R_{chamber}"), scale_to_string(radius, 1, mm, "mm"))
        << vspace
        << "Giving:\n\n"
        << svspace
        << make_aligned_eqn(
            area_var,
            area_eqn,
            scale_to_string(area, 1, mm_sq, "mm+2")
        )
        << svspace
        << "Giving a contraction ratio of:\n\n"
        << svspace
        << make_eqn(
            contraction_eqn,
            contraction_ratio
        )
        << vspace
        << "We can now find a length for the chamber:\n\n"
        << svspace
        << make_aligned_eqn(
            StyledVar("L_{chamber}"),
            DIV(volume_var, area_var),
            DIV(
                scale_to_string(volume_var.solve(), 1, mm_cu, "mm+3"),
                scale_to_string(area_var.solve(), 1, mm_sq, "mm+2")
            ),
            scale_to_string(length, 1, mm, "mm")
        )
        << vspace
        << "Via trigonometry we find the converging section length:\n\n"
        << svspace
        << make_aligned_eqn(
            height_var,
            converging_height_eqn,
            scale_to_string(converging_height, 1, mm, "mm")
        )
        << vspace
        << make_aligned_eqn(
            StyledVar("L_{conv}"),
            converging_length_eqn,
            scale_to_string(converging_length, 1, mm, "mm")
        )
        << svspace
        << "Which then yields:\n\n"
        << svspace
        << make_aligned_eqn(
            StyledVar("L_{flatwall}"),
            flat_eqn,
            scale_to_string(flat_eqn.solve(), 1, mm, "mm")
        )
    ;

    Chamber chamber{radius, length, converging_length, converging_height, throat.radius};

    section
        << critical_vars
        << geometry
        << (Subsection("Summary")
            << vspace
            << (UnorderedList()
                << strjoin(PrefixText("Chamber Length:"), scale_to_string(length, 1, mm, "mm"))
                << strjoin(PrefixText("Chamber Diameter:"), scale_to_string(diameter, 1, mm, "mm"))
                << strjoin(PrefixText("Chamber Area:"), scale_to_string(area, 1, mm_sq, "mm+2"))
                << strjoin(PrefixText("Chamber Flatwall Length:"), scale_to_string(length-converging_length, 1, mm, "mm"))
                << strjoin(PrefixText("Converging Section Height:"), scale_to_string(converging_height, 1, mm, "mm"))
                << strjoin(PrefixText("Converging Section Length:"), scale_to_string(converging_length, 1, mm, "mm"))
            )
            << svspace
            << (UnorderedList()
                << strjoin(PrefixText("Length/Width Ratio:"), length/diameter)
                << strjoin(PrefixText("Contraction Ratio:"), contraction_ratio)
            )
            << svspace
            << (UnorderedList()
                << strjoin(PrefixText("Flatwall Volume:"),              scale_to_string(chamber.flatwall_volume(), 1, cm_cu, "cm+3"))
                << strjoin(PrefixText("Converging Section Volume:"),    scale_to_string(chamber.converging_volume(), 1, cm_cu, "cm+3"))
                << strjoin(PrefixText("Combustible Volume:"),           scale_to_string(chamber.volume(), 1, cm_cu, "cm+3"))
            )
        )
    ;

    doc << section;
    return chamber;
}


#endif
