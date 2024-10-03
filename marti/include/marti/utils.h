#ifndef MARTI_INCLUDE_MARTI_UTILS_H_
#define MARTI_INCLUDE_MARTI_UTILS_H_

#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>

#include "marti/definitions.h"

namespace marti::utils {

    // String conversion utilities for printing
    template<typename T>
    inline auto to_string(const T &t) -> typename std::enable_if<std::is_same<T, std::string>::value, std::string>::type {
        return static_cast<std::string>(t);
    }

    template<typename T>
    inline auto to_string(const T &t) -> typename std::enable_if<!std::is_same<T, std::string>::value, std::string>::type {
        return std::to_string(t);
    }

    inline std::string to_string(const std::string_view &sv) {
        return std::string(sv);
    }

    template<typename T>
    inline std::string implode(const T& list, const std::string& delimiter=",") {
        if (list.empty()) return "";
        return std::accumulate(std::next(list.begin()), list.end(), to_string(*list.begin()),
                               [&delimiter](const std::string &a, const auto &b) {
                                   return a + delimiter + to_string(b);
                               }
        );
    }

    // Convert a map to a newline-separated string of key:value pairs
    template<typename T>
    inline std::string implode_map(const T &map) {
        if (map.empty()) return "";
        return std::accumulate(std::next(map.begin()), map.end(),
                               to_string(map.begin()->first) + ": " + to_string(map.begin()->second),
                               [](const std::string &a, const auto &b) {
                                   return a + "\n" + to_string(b.first) + ": " + to_string(b.second);
                               }
        );
    }


    // String/sequence utilities

    // Compute the reverse complement of a DNA sequence
    inline std::string reverse_complement(const std::string &seq) {
        std::string rc(seq.rbegin(), seq.rend());
        std::transform(rc.begin(), rc.end(), rc.begin(),
                       [](char &c) {return kDNAComplement.at(static_cast<char>(std::toupper(c)));});
        return rc;
    }

    // Pad a string with spaces on the right until it is at least max_width long
    inline std::string ljust(const std::string& s, int width) {
        std::ostringstream strm;
        strm.width(width);
        strm << std::left << s;
        return strm.str();
    }

    // Convenience methods
    template<class C, typename T>
    inline bool contains(C&& c, T e) { return std::find(std::begin(c), std::end(c), e) != std::end(c); };
}


#endif //MARTI_INCLUDE_MARTI_UTILS_H_
