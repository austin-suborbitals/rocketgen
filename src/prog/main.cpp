#include "ext/cxxopts/src/cxxopts.hpp"
#include "ext/json/src/json.hpp"

#include "rocketgen.hpp"

using namespace nlohmann; // get access to bare json:: namespace

void sanitize_options(cxxopts::Options& opts) {
    if (opts.count("help")) {
        std::cout << opts.help({"General", "Input", "Output"}) << std::endl;
        exit(0);
    }

    if (not opts.count("config")) {
        std::cout << "expected a config to be given!\n";
        exit(1); // TODO: proper return codes
    }
}

cxxopts::Options handle_options(int argc, char** argv) {
    cxxopts::Options opts("rocketgen", "paramteric rocket engine designer");
    opts.add_options("General")
        ("h,help", "prints this help dialog")
        ;
    opts.add_options("Input")
        ("c,config", "config file (JSON) containing the needed parameters", cxxopts::value<std::string>(), "FILE")
        ;
    opts.add_options("Output")
        ;

    opts.parse(argc, argv);
    sanitize_options(opts);

    return opts;
}

EngineBasis marshal_config(const json& conf) {
    Fuel fuel{
        conf["propellant"]["fuel"]["name"].get<std::string>(),
        conf["propellant"]["fuel"]["molar_mass"].get<double>()
    };

    Oxidizer oxidizer{
        conf["propellant"]["oxidizer"]["name"].get<std::string>(),
        conf["propellant"]["oxidizer"]["molar_mass"].get<double>()
    };

    auto thrust = (conf["engine"]["thrust"].get<double>() * newton);
    auto isp = conf["engine"]["isp"].get<double>();

    return EngineBasis{
        conf["name"].get<std::string>(),
        conf["version"].get<std::string>(),
        thrust, isp*second,
        (conf["engine"]["chamber_pressure"].get<double>() * MPa),
        (conf["engine"]["external_pressure"].get<double>() * MPa),
        conf["engine"]["expansion_ratio"].get<double>(),
        (conf["engine"]["l_star"].get<double>() * meter),
        TORAD(conf["engine"]["converging_angle"].get<double>()),  // we take in as deg, but calcs need rad
        TORAD(conf["engine"]["diverging_angle"].get<double>()),  // we take in as deg, but calcs need rad
        Propellants{
            fuel, oxidizer,
            conf["propellant"]["mixture_ratio"].get<double>(),
            conf["propellant"]["gamma"].get<double>(),
            conf["propellant"]["avg_combustion_weight"].get<double>() * gram,
            (conf["propellant"]["flame_temp"].get<double>() + 273.15) * kelvin
        }
    };
}

void build_overview(latex::doc::Report& doc, const EngineBasis& rocket) {
    using PrefixText = latex::Text<latex::style::Large, latex::style::Bold>;

    auto mass_flow = calc_mass_flow(rocket.thrust, rocket.isp);
    auto fuel_flow_rate = mass_flow / (rocket.propellants.mixture_ratio + 1); // TODO: +1?

    latex::doc::Section overview("Overview", true);
    overview
        << "A brief overview of the target-metrics for the engine." << "\n\n"
        << vspace
        << (latex::doc::UnorderedList()
                << strjoin(PrefixText("Thrust: "), rocket.thrust)
                << strjoin(PrefixText("Isp: "), rocket.isp)
                << strjoin(PrefixText("Chamber Pressure: "), eng::to_string(rocket.pressure))
                << PrefixText("Propellants: ")
                << (latex::doc::UnorderedList()
                    << strjoin(PrefixText("Mixture Ratio: "), rocket.propellants.mixture_ratio)
                    << strjoin(PrefixText("Total Mass Flow Rate: "), mass_flow)
                    << PrefixText("Fuel: ")
                    << (latex::doc::UnorderedList()
                        << strjoin(PrefixText("Type: "), rocket.propellants.fuel.name)
                        << strjoin(PrefixText("Mass Flow Rate: "), fuel_flow_rate)
                    )
                    << PrefixText("Oxidizer: ")
                    << (latex::doc::UnorderedList()
                        << strjoin(PrefixText("Type: "), rocket.propellants.oxidizer.name)
                        << strjoin(PrefixText("Mass Flow Rate: "), (mass_flow - fuel_flow_rate))
                    )
                )
           );

    doc << overview;
}


int main(int argc, char** argv) {
    auto opts = handle_options(argc, argv);
    auto conf_str = read_file(opts["config"].as<std::string>());
    auto rocket = marshal_config(json::parse(conf_str));

    auto doc = Report(
        strjoin(rocket.name, rocket.version),
        strjoin(
            rocket.thrust,
            rocket.propellants.fuel.name, "/", rocket.propellants.oxidizer.name,
            "Liquid Rocket Engine"
        )
    )
    .use("amsmath")
    .use("gensymb");

    // build and append the overview
    build_overview(doc, rocket);
    auto throat = build_throat(doc, rocket);
    build_nozzle(doc, rocket, throat);
    build_chamber(doc, rocket, throat);

    std::cout << std::setprecision(4) << doc << std::endl;

    return 0;
}
