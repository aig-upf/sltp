
#include "equivalences.h"


namespace sltp::cnf {

std::size_t hash_value(const transition_denotation& x) {
    std::size_t seed = 0;
    boost::hash_combine(seed, x.value_s);
    boost::hash_combine(seed, x.increases);
    boost::hash_combine(seed, x.decreases);
    return seed;
}

} // namespaces


