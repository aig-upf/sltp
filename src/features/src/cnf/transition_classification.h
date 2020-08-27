
#pragma once

#include "generator.h"

//#include <boost/container/flat_set.hpp>

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
            class_representatives_(),
            from_transition_to_eq_class_(),
            types_(),
            necessarily_bad_transitions_()
    {
        compute_equivalence_relations();
    }


    bool write(CNFWriter& wr) override;

    inline unsigned get_transition_id(uint16_t s, uint16_t t) const { return transition_ids_.at(state_pair(s, t)); }

    inline unsigned get_representative_id(unsigned tx) const { return from_transition_to_eq_class_.at(tx); }

    inline unsigned get_class_representative(uint16_t s, uint16_t t) const {
        return get_representative_id(get_transition_id(s, t));
    }

    inline const state_pair& get_state_pair(unsigned tx) const { return transition_ids_inv_.at(tx); }

    inline bool is_necessarily_bad(unsigned tx) const {
        return necessarily_bad_transitions_.find(tx) != necessarily_bad_transitions_.end();
    }

    inline int get_max_v(unsigned s) const {
        int vstar = sample_.transitions().vstar(s);
        if (vstar < 0) return -1;
        return std::ceil(options.v_slack * vstar);
    }

    inline unsigned compute_D() const {
        // return 20;
        // D will be the maximum over the set of alive states of the upper bounds on V_pi
        unsigned D = 0;
        for (const auto s:all_alive()) {
            auto max_v_s = get_max_v(s);
            if (max_v_s > D) D = max_v_s;
        }
        return D;
    }

    bool check_validity_of_existing_solution();

protected:
    // A mapping from pairs of state to the assigned transition id
    std::unordered_map<state_pair, unsigned, boost::hash<state_pair>> transition_ids_;
    // The reverse mapping
    std::vector<state_pair> transition_ids_inv_;

    // A list of transition IDs that are the representative of their class
    std::vector<unsigned> class_representatives_;

    // A mapping from the ID of the transition to the ID of its equivalence class
    std::vector<unsigned> from_transition_to_eq_class_;

    // A mapping from the ID of the transition to its transition type
    std::vector<transition_type> types_;

    std::unordered_set<unsigned> necessarily_bad_transitions_;

    //!
    void compute_equivalence_relations();

    //!
    std::vector<bool> check_feature_dominance();


    //! Return DT(f), the set of pairs of transitions that are distinguished by the given feature f
    std::vector<transition_pair> compute_dt(unsigned f);

};

} // namespaces

