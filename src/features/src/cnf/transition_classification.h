
#pragma once

#include "generator.h"
#include "equivalences.h"

namespace sltp::cnf {

class TransitionClassificationEncoding : public CNFEncoding {
public:
    enum class transition_type : bool {
        alive_to_solvable,
        alive_to_dead
    };

    using state_pair = std::pair<uint16_t, uint16_t>;

    TransitionClassificationEncoding(const Sample::Sample& sample, const Options& options) :
            CNFEncoding(sample, options),
            transition_ids_(),
            transition_ids_inv_(),
            from_trace_to_class_repr_(),
            from_transition_to_eq_class_(),
            types_()
    {
        compute_equivalence_relations();
    }


    std::pair<bool, CNFWriter> write(std::ostream &os) override;


protected:
    // A mapping from pairs of state to the assigned transition id
    std::unordered_map<state_pair, unsigned, boost::hash<state_pair>> transition_ids_;
    // The reverse mapping
    std::vector<state_pair> transition_ids_inv_;

    // A mapping from a full transition trace to the ID of the corresponding equivalence class
    std::unordered_map<transition_trace, unsigned> from_trace_to_class_repr_;

    // A mapping from the ID of the transition to the ID of its equivalence class
    std::vector<unsigned> from_transition_to_eq_class_;

    // A mapping from the ID of the transition to its transition type
    std::vector<transition_type> types_;

    //!
    void compute_equivalence_relations();

    //!
    std::vector<bool> check_feature_dominance();


    //!
    boost::container::flat_set<transition_pair> compute_dt(unsigned f);

};

} // namespaces

