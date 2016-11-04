#ifndef RGEN_BIN_UTILS
#define RGEN_BIN_UTILS


#include <string>
#include <fstream>


std::string read_file(const std::string& path){
    std::ifstream ifs(path);
    return std::string((std::istreambuf_iterator<char>(ifs)),
                  (std::istreambuf_iterator<char>()));
}


/** Sequentially join strings with the given separator (first template arg).
 *
 * All arguments must be ostream compatible.
 *
 * The null character '\0' as a separator acts as a sentinel for no separator.
 */
template <char sep = ' ', typename... Args>
std::string strjoin(const Args& ...args) {
    std::stringstream ss;
    using expander = int[];
    if (sep == '\0') { // TODO: gross
        (void) expander{ (ss << std::forward<const Args&>(args), void(), 0)... };
    } else {
        (void) expander{ (ss << std::forward<const Args&>(args) << sep, void(), 0)... };
    }
    return ss.str();
}


#endif
