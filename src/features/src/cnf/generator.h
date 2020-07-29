
#pragma once

#include <queue>
#include <blai/sample.h>
#include "cnfwriter.h"

//! A feature index
using feature_t = uint32_t;

namespace sltp::cnf {

struct Options {
    enum class Encoding {
        Basic,
        D2Tree,
        TransitionSeparation
    };

    std::string workspace;
    bool prune_redundant_states;
    bool verbose;
    Encoding encoding;
    bool use_only_alive_unmarked_states;
    bool distinguish_transitions_locally;
    std::vector<unsigned> enforced_features;


    bool use_d2tree() const { return encoding == Encoding::D2Tree; }
    bool use_separation_encoding() const { return encoding == Encoding::TransitionSeparation; }
};

} // namespaces

inline void undist_goal_warning(unsigned s, unsigned t) {
    std::cout << Utils::warning()
        <<  "No feature can distinguish state " << s << " from state " << t << ", but only one of them is a goal"
        <<  ". The MAXSAT encoding will be UNSAT" << std::endl;
}

inline void undist_deadend_warning(unsigned s, unsigned t) {
    std::cout << Utils::warning()
        <<  "No feature can distinguish state " << s << " from state " << t << ", but (only) one of them is a"
        <<  " dead-end. The MAXSAT encoding will be UNSAT" << std::endl;
}

//! Return a sorted vector with those features that d1-distinguish s from t
std::vector<feature_t> compute_d1_distinguishing_features(const Sample::Sample& sample, unsigned s, unsigned t);

//! Return a sorted vector with those features that d2-distinguish transition (s, s') from (t, t')
std::vector<feature_t> compute_d2_distinguishing_features(const Sample::Sample& sample,
        unsigned s, unsigned sprime, unsigned t, unsigned tprime);


class CNFGenerator {
public:
    //!
    using d1_key = std::tuple<unsigned, unsigned>;

    //! A map from 4 state IDs to the ID of the corresponding D2(s, s', t, t') variable
    using d2_key = std::tuple<unsigned, unsigned, unsigned, unsigned>;
    using d2map_t = std::unordered_map<d2_key, cnfvar_t, boost::hash<d2_key>>;

    using transition_t = Sample::TransitionSample::transition_t;
    using transition_set_t = Sample::TransitionSample::transition_set_t;
    using transition_list_t = Sample::TransitionSample::transition_list_t;

protected:
    //! The transition sample data
    const Sample::Sample& sample_;

    //! The CNF encoding options
    const sltp::cnf::Options& options;

    //! The number of states in the encoding
    const std::size_t ns_;

    //! The number of features in the encoding
    const std::size_t nf_;

    //! For convenient and performant access, a list of goal and non-goal states
    std::vector<unsigned> goals_, nongoals_;

    d2map_t d2ids_;
    std::vector<cnfvar_t> d2vars_;

    //!
    std::unordered_map<d1_key, std::vector<feature_t>, boost::hash<d1_key>> d1_features_cache_;
    std::vector<std::vector<feature_t>> d2_features_cache_;

public:
    //! Some statistics
    unsigned n_selected_clauses;
    unsigned n_d2_clauses;
    unsigned n_separation_clauses;
    unsigned n_bridge_clauses;
    unsigned n_goal_clauses;
    unsigned n_deadend_clauses;
    unsigned n_good_tx_clauses;
    unsigned n_leq_clauses;

    explicit CNFGenerator(const Sample::Sample& sample, const sltp::cnf::Options& options) :
        sample_(sample),
        options(options),
        ns_(sample.matrix().num_states()),
        nf_(sample.matrix().num_features()),
        n_selected_clauses(0),
        n_d2_clauses(0),
        n_separation_clauses(0),
        n_bridge_clauses(0),
        n_goal_clauses(0),
        n_deadend_clauses(0),
        n_good_tx_clauses(0),
        n_leq_clauses(0)
    {
        for (unsigned s = 0; s < ns_; ++s) {
            if (is_goal(s)) goals_.push_back(s);
            else nongoals_.push_back(s);
        }
    }

    const std::vector<unsigned>& all_alive() const { return sample_.transitions().all_alive(); }

    bool is_goal(unsigned s) const { return sample_.matrix().goal(s); }

    bool is_alive(unsigned s) const { return sample_.transitions().is_alive(s); }

    bool is_solvable(unsigned s) const { return is_alive(s) || is_goal(s); }

    bool is_sound_transition(unsigned s, unsigned sprime) const {
        return is_marked_transition(s, sprime);
    }

    bool is_marked_transition(unsigned s, unsigned sprime) const {
        return sample_.transitions().marked(s, sprime);
    }

    const transition_set_t& sound_transitions() const {
        return marked_transitions();
    }

    const transition_set_t& marked_transitions() const {
        return sample_.transitions().marked_transitions();
    }

    const transition_list_t& unmarked_transitions() const {
        return sample_.transitions().unmarked_transitions();
    }

    const transition_list_t& unmarked_and_alive_transitions() const {
        return sample_.transitions().unmarked_and_alive_transitions();
    }

    const transition_list_t& get_relevant_unmarked_transitions(unsigned s) const {
        if (options.distinguish_transitions_locally) {
            assert(is_alive(s));
            return sample_.transitions().unmarked_transitions_starting_at(s);
        } else {
            if (options.use_only_alive_unmarked_states) {
                return sample_.transitions().unmarked_and_alive_transitions();
            } else {
                return sample_.transitions().unmarked_transitions();
            }
        }
    }

    unsigned feature_weight(unsigned f) const {
        return sample_.matrix().feature_cost(f);
    }

    const std::vector<unsigned>& successors(unsigned s) const {
        return sample_.transitions().successors(s);
    }

    //! Return possibly-cached set of features that d1-distinguish s from t
    const std::vector<feature_t>& d1_distinguishing_features(unsigned s, unsigned t) {
        const auto idx = d1idx(s, t);
        const auto it = d1_features_cache_.find(idx);
        if (it != d1_features_cache_.end()) return it->second;

        std::vector<unsigned>& features = d1_features_cache_[idx];
        features = compute_d1_distinguishing_features(sample_, s, t);
        return features;
    }

    //! Return possibly-cached set of features that d2-distinguish (s, s') from (t, t')
    const std::vector<feature_t>& d2_distinguishing_features(unsigned d2id) {
        // key is guaranteed to be on the cache, as we have inserted upfront all necessary sets
        return d2_features_cache_.at(d2id);
    }


    static d1_key d1idx(unsigned s, unsigned t) {
        assert(s != t);
        return (s < t) ? std::make_tuple(s, t) : std::make_tuple(t, s);
    }

    static d2_key d2idx(unsigned s, unsigned sprime, unsigned t, unsigned tprime) {
        assert(s != t);
        return (s < t) ? std::make_tuple(s, sprime, t, tprime) :
               std::make_tuple(t, tprime, s, sprime);
    }

    cnfvar_t get_d2var(unsigned s, unsigned sprime, unsigned t, unsigned tprime) const {
        auto d2id = d2ids_.at(d2idx(s, sprime, t, tprime));
        return d2vars_.at(d2id);
    }

    //! Write the desired CNF encoding to the specified file
    std::pair<bool, CNFWriter> write_encoding(std::ofstream& os);

    //! Generate and write the CNF instance for the standard encoding as we go
    std::pair<bool, CNFWriter> write_basic_maxsat(std::ostream &os);

    //! Generate and write the CNF instance for the transition-separation encoding as we go
    std::pair<bool, CNFWriter> write_transition_classification_maxsat(std::ostream &os);
};

using isomorphism_t = std::unordered_map<unsigned, unsigned>;

//! Check whether all transitions starting in s have some transition starting in t with same qualitative nature
//! on the set of all features in the given feature matrix
bool all_tx_have_analogs(const Sample::Sample& sample, unsigned s, unsigned t);

//! Check whether t appears isomorphic to s, and in that case, add it to the given list of isomorphisms
void check_isomorphic(const Sample::Sample& sample, unsigned s, unsigned t, isomorphism_t& isomorphisms);

//!
isomorphism_t compute_redundant_states(const Sample::Sample& sample);