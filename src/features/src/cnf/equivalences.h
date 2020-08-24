
#pragma once

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>

namespace sltp::cnf {

struct transition_denotation {
    uint8_t value_s:6;
    uint8_t increases:1;
    uint8_t decreases:1;

    transition_denotation()
        : value_s(0), increases(0), decreases(0)
    {}

    transition_denotation(unsigned value_s, int sign)
        : value_s(value_s), increases((sign>0) ? 1 : 0), decreases((sign<0) ? 1 : 0)
    {
        if (value_s >= 64) {  // i.e. 2^6-1
            throw std::runtime_error("Value given for transition_denotation is too large - Change the struct definition and recompile");
        }
    }

    [[nodiscard]] bool nils() const { return increases == 0 && decreases == 0; }
};
static_assert(sizeof(transition_denotation) == 1);

inline bool operator==(const transition_denotation& x, const transition_denotation& y) {
    return x.value_s == y.value_s && x.increases == y.increases && x.decreases == y.decreases;
}


std::size_t hash_value(const transition_denotation& x);




struct transition_trace {
//    std::array<transition_denotation, 40> data; // This could be useful if we decide to compile the code on the fly

    //! A list of denotations, one denotation for every transition
    std::vector<transition_denotation> denotations;

    explicit transition_trace(unsigned nfeatures)
        : denotations(nfeatures)
    {}
};

inline bool operator==(const transition_trace& x, const transition_trace& y) { return x.denotations == y.denotations; }

} // namespaces

// Specializations of std::hash
namespace std {

template<> struct hash<sltp::cnf::transition_denotation> {
    std::size_t operator()(const sltp::cnf::transition_denotation& x) const {
        return hash_value(x);
    }
};

template<> struct hash<sltp::cnf::transition_trace> {
    std::size_t operator()(const sltp::cnf::transition_trace& x) const {
        return boost::hash_range(x.denotations.begin(), x.denotations.end());
    }
};
}