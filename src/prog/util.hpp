#ifndef RGEN_BIN_UTILS
#define RGEN_BIN_UTILS


#include <string>
#include <fstream>


std::string read_file(const std::string& path){
    std::ifstream ifs(path);
    return std::string((std::istreambuf_iterator<char>(ifs)),
                  (std::istreambuf_iterator<char>()));
}


#endif
