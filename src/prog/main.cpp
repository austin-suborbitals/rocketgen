#include "ext/cxxopts/src/cxxopts.hpp"
#include "ext/json/src/json.hpp"

#include "lib/rocketgen.hpp"

#include "util.hpp"

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

    return EngineBasis(
        conf["name"].get<std::string>(),
        conf["version"].get<std::string>(),
        conf["engine"]["thrust"].get<double>(),
        conf["engine"]["isp"].get<double>(),
        conf["engine"]["chamber_pressure"].get<double>(),
        Propellants{fuel, oxidizer, conf["propellant"]["mixture_ratio"].get<double>()}
    );
}

void build_overview(latex::doc::Report& doc, const EngineBasis& rocket) {
    using namespace latex;
    using PrefixText = Text<style::Larger, style::Bold>;

    auto total_flow_rate = rocket.thrust / (GRAVITY * rocket.isp);
    auto fuel_flow_rate = total_flow_rate / (rocket.propellants.mixture_ratio + 1); // TODO: +1?

    latex::doc::Section overview("Overview");
    overview
        << "A brief overview of the target-metrics for the engine." << "\n\n"
        << (latex::doc::UnorderedList()
                << strjoin(PrefixText("Thrust: "), rocket.thrust)
                << strjoin(PrefixText("Isp: "), rocket.isp)
                << strjoin(PrefixText("Chamber Pressure: "), eng::to_string(rocket.pressure))
                << PrefixText("Propellants: ")
                << (doc::UnorderedList()
                    << strjoin(PrefixText("Mixture Ratio: "), rocket.propellants.mixture_ratio)
                    << strjoin(PrefixText("Total Mass Flow Rate: "), total_flow_rate)
                    << PrefixText("Fuel: ")
                    << (doc::UnorderedList()
                        << strjoin(PrefixText("Type: "), rocket.propellants.fuel.name)
                        << strjoin(PrefixText("Mass Flow Rate: "), fuel_flow_rate)
                    )
                    << PrefixText("Oxidizer: ")
                    << (doc::UnorderedList()
                        << strjoin(PrefixText("Type: "), rocket.propellants.oxidizer.name)
                        << strjoin(PrefixText("Mass Flow Rate: "), (total_flow_rate - fuel_flow_rate))
                    )
                )
           );

    doc << overview;
}


int main(int argc, char** argv) {
    auto opts = handle_options(argc, argv);
    auto conf_str = read_file(opts["config"].as<std::string>());
    auto rocket = marshal_config(json::parse(conf_str));

    latex::doc::Report doc(
        strjoin(rocket.name, rocket.version),
        strjoin(
            rocket.thrust,
            rocket.propellants.fuel.name, "/", rocket.propellants.oxidizer.name,
            "Liquid Rocket Engine"
        )
    );

    // build and append the overview
    build_overview(doc, rocket);

    std::cout << doc << std::endl;

    return 0;
}
