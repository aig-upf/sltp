
#pragma once

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

    //! The path of the workspace where output files will be left
    std::string workspace;

    //! The type of encoding we want to use
    Encoding encoding;

    //! In the transition-separation encoding, whether we want to exploit the equivalence relation
    //! among transitions given by the feature pool
    bool use_equivalence_classes;

    //! In the transition-separation encoding, whether we want to exploit the dominance among features to ignore
    //! dominated features and reduce the size of the encoding.
    bool use_feature_dominance;

    //! In the AAAI19 encoding, whether to prune states that appear redundant for the given feature pool
    bool prune_redundant_states;

    //! Whether to be more verbose in the generation of the encoding
    bool verbose;

    //! The slack value for the maximum allowed value for V_\pi(s) = slack * V^*(s)
    double v_slack;

    //! In the transition-separation encoding, whether to use the incremental refinement approach
    bool use_incremental_refinement;

    //! In the transition-separation CNF encoding, whether to distinguish good transitions *only from*
    //! unmarked transitions that start in the same state as the good transition
    bool distinguish_transitions_locally;

    //! A list of user-provided feature IDs for which we want to enforce selection
    std::vector<unsigned> enforced_features;


    [[nodiscard]] bool use_d2tree() const { return encoding == Encoding::D2Tree; }
    [[nodiscard]] bool use_separation_encoding() const { return encoding == Encoding::TransitionSeparation; }
};


enum class CNFGenerationOutput : unsigned {
    Success = 0,
    UnsatTheory = 1,
    ValidationCorrectNoRefinementNecessary = 2
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

//! Return a sorted vector with those features that d2-distinguish transition (s, s') from (t, t')
std::vector<feature_t> compute_d1d2_distinguishing_features(const Sample::Sample& sample,
                                                          unsigned s, unsigned sprime, unsigned t, unsigned tprime);

bool are_transitions_d1d2_distinguished(int s_f, int sprime_f, int t_f, int tprime_f);

class CNFEncoding {
public:

    using transition_t = Sample::TransitionSample::transition_t;
    using transition_set_t = Sample::TransitionSample::transition_set_t;
    using transition_list_t = Sample::TransitionSample::transition_list_t;

    CNFEncoding(const Sample::Sample& sample, const sltp::cnf::Options& options) :
        sample_(sample),
        options(options),
        ns_(sample.matrix().num_states()),
        nf_(sample.matrix().num_features())
    {
        for (unsigned s = 0; s < ns_; ++s) {
            if (is_goal(s)) goals_.push_back(s);
            else nongoals_.push_back(s);
        }
    }

    [[nodiscard]] const std::vector<unsigned>& all_alive() const { return sample_.transitions().all_alive(); }

    [[nodiscard]] bool is_goal(unsigned s) const { return sample_.matrix().goal(s); }

    [[nodiscard]] bool is_alive(unsigned s) const { return sample_.transitions().is_alive(s); }

    [[nodiscard]] bool is_solvable(unsigned s) const { return is_alive(s) || is_goal(s); }

    [[nodiscard]] unsigned feature_weight(unsigned f) const {
        return sample_.matrix().feature_cost(f);
    }

    [[nodiscard]] const std::vector<unsigned>& successors(unsigned s) const {
        return sample_.transitions().successors(s);
    }

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

};


using isomorphism_t = std::unordered_map<unsigned, unsigned>;

//! Check whether all transitions starting in s have some transition starting in t with same qualitative nature
//! on the set of all features in the given feature matrix
bool all_tx_have_analogs(const Sample::Sample& sample, unsigned s, unsigned t);

//! Check whether t appears isomorphic to s, and in that case, add it to the given list of isomorphisms
void check_isomorphic(const Sample::Sample& sample, unsigned s, unsigned t, isomorphism_t& isomorphisms);

//!
isomorphism_t compute_redundant_states(const Sample::Sample& sample);