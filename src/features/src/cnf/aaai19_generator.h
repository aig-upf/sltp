
#pragma once

#include "generator.h"

namespace sltp::cnf {


class AAAI19Generator : public CNFEncoding {
public:
    //!
    using d1_key = std::tuple<unsigned, unsigned>;

    //! A map from 4 state IDs to the ID of the corresponding D2(s, s', t, t') variable
    using d2_key = std::tuple<unsigned, unsigned, unsigned, unsigned>;
    using d2map_t = std::unordered_map<d2_key, cnfvar_t, boost::hash<d2_key>>;


protected:
    d2map_t d2ids_;
    std::vector<cnfvar_t> d2vars_;

    //!
    std::unordered_map<d1_key, std::vector<feature_t>, boost::hash<d1_key>> d1_features_cache_;
    std::vector<std::vector<feature_t>> d2_features_cache_;

public:
    //! Some statistics
    unsigned n_selected_clauses;
    unsigned n_d2_clauses;
    unsigned n_bridge_clauses;
    unsigned n_goal_clauses;
    unsigned n_deadend_clauses;

    AAAI19Generator(const TrainingSet& sample, const sltp::cnf::Options& options) :
            CNFEncoding(sample, options),
            n_selected_clauses(0),
            n_d2_clauses(0),
            n_bridge_clauses(0),
            n_goal_clauses(0),
            n_deadend_clauses(0)
    {}

    static TrainingSet preprocess_sample(const TrainingSet& sample, const sltp::cnf::Options& options);


    bool is_sound_transition(unsigned s, unsigned sprime) const {
        return is_marked_transition(s, sprime);
    }

    bool is_marked_transition(unsigned s, unsigned sprime) const {
        return tr_set_.transitions().marked(s, sprime);
    }

    const transition_set_t& sound_transitions() const {
        return marked_transitions();
    }

    const transition_set_t& marked_transitions() const {
        return tr_set_.transitions().marked_transitions();
    }

    const transition_list_t& unmarked_transitions() const {
        return tr_set_.transitions().unmarked_transitions();
    }

    //! Return possibly-cached set of features that d1-distinguish s from t
    const std::vector<feature_t>& d1_distinguishing_features(unsigned s, unsigned t) {
        const auto idx = d1idx(s, t);
        const auto it = d1_features_cache_.find(idx);
        if (it != d1_features_cache_.end()) return it->second;

        std::vector<unsigned>& features = d1_features_cache_[idx];
        features = compute_d1_distinguishing_features(tr_set_, s, t);
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

    sltp::cnf::CNFGenerationOutput write(CNFWriter& wr);
};

} // namespaces

