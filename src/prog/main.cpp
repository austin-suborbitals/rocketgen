#include "../ext/cxxopts/src/cxxopts.hpp"
#include "../ext/json/src/json.hpp"

#include "util.hpp"
#include "../lib/rocketgen.hpp"

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

int main(int argc, char** argv) {
    auto opts = handle_options(argc, argv);
    auto conf_str = read_file(opts["config"].as<std::string>());
    auto conf = json::parse(conf_str);

    auto thrust = conf["engine"]["thrust"].get<double>() * newton;
    auto isp = conf["propellant"]["expected_isp"].get<double>() * second;

    auto total_flow_rate = thrust / (isp * GRAVITY);
    auto fuel_flow_rate = total_flow_rate / (conf["propellant"]["mixture_ratio"].get<double>() + 1); // TODO: +1?

    std::cout << "Total propellant flow rate: " << total_flow_rate << std::endl;
    std::cout << "Fuel flow rate: " << fuel_flow_rate << std::endl;
    std::cout << "Oxidizer flow rate: " << (total_flow_rate - fuel_flow_rate) << std::endl;

    return 0;
}
