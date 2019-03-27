
#pragma once

#include <blai/sample.h>
#include "cnfwriter.h"

void undist_goal_warning(unsigned s, unsigned t) {
    std::cout <<  "No feature can distinguish state " << s << " from state " << t << ", but only one of tthem is a goal"
              <<  ". The MAXSAT encoding will be UNSAT" << std::endl;
}


class CNFGenerator {
public:
    //!
    using d1_key = std::tuple<unsigned, unsigned>;

    //! A map from 4 state IDs to the ID of the corresponding D2(s, s', t, t') variable
    using d2_key = std::tuple<unsigned, unsigned, unsigned, unsigned>;
    using d2map_t = std::unordered_map<d2_key, cnfvar_t, boost::hash<d2_key>>;

    using transition_t = Sample::Transitions::transition_t;
    using transition_set_t = Sample::Transitions::transition_set_t;
    using transition_list_t = Sample::Transitions::transition_list_t;

protected:
    //! The transition sample data
    const Sample::Sample& sample_;

    //! The number of states in the encoding
    const unsigned ns_;

    //! The number of features in the encoding
    const unsigned nf_;

    //! For convenient and performant access, a list of goal and non-goal states + a list of expanded states
    std::vector<unsigned> goals_, nongoals_;
    std::vector<unsigned> expanded_states_;

    //!
    std::unordered_map<d1_key, std::vector<unsigned>, boost::hash<d1_key>> d1_features_cache_;
    std::unordered_map<d2_key, std::vector<unsigned>, boost::hash<d2_key>> d2_features_cache_;

public:
    //! Some statistics
    unsigned n_selected_clauses;
    unsigned n_d2_clauses;
    unsigned n_bridge_clauses;
    unsigned n_goal_clauses;

public:
    explicit CNFGenerator(const Sample::Sample& sample) :
        sample_(sample),
        ns_(sample.matrix().num_states()),
        nf_(sample.matrix().num_features()),
        n_selected_clauses(0),
        n_d2_clauses(0),
        n_bridge_clauses(0),
        n_goal_clauses(0)
    {
        for (unsigned s = 0; s < ns_; ++s) {
            if (is_goal(s)) goals_.push_back(s);
            else nongoals_.push_back(s);

            if (sample_.transitions_.num_transitions(s) > 0) {
                // TODO This check is not fully correct, we should annotate the planner output with this information
                // TODO so that we do not take a dead-end for a not expanded state
                expanded_states_.push_back(s);
            }
        }
    }

    bool is_goal(unsigned s) const { return sample_.matrix().goal(s); }

    bool is_sound_transition(unsigned s, unsigned sprime) const {
        return sample_.transitions().marked(s, sprime);
    }

    const transition_set_t& sound_transitions() const {
        return sample_.transitions().marked_transitions();
    }

    const transition_list_t& all_transitions() const {
        return sample_.transitions().all_transitions();
    }

    unsigned feature_weight(unsigned f) const {
        return sample_.matrix().feature_cost(f);
    }

    const std::vector<unsigned>& successors(unsigned s) const {
        return sample_.transitions().successors(s);
    }

    const std::vector<unsigned>& d1_distinguishing_features(unsigned s, unsigned t) {
        const auto idx = d1idx(s, t);
        const auto it = d1_features_cache_.find(idx);
        if (it != d1_features_cache_.end()) return it->second;

        std::vector<unsigned>& features = d1_features_cache_[idx];
        for(unsigned f = 0; f < nf_; ++f) {
            // Store those features that d1-distinguish s from t
            auto sf = sample_.matrix().entry(s, f);
            auto tf = sample_.matrix().entry(t, f);
            if ((sf == 0) != (tf == 0)) {
                features.push_back(f);
            }
        }
        return features;
    }

    const std::vector<unsigned>& d2_distinguishing_features(const d2_key& key) {
        unsigned s = std::get<0>(key), sprime = std::get<1>(key), t = std::get<2>(key), tprime = std::get<3>(key);
        const auto it = d2_features_cache_.find(key);
        if (it != d2_features_cache_.end()) return it->second;

        std::vector<unsigned>& features = d2_features_cache_[key];
        for(unsigned f = 0; f < nf_; ++f) {
            // Store those features that d2-distinguish (s, s') from (t, t'), but do _not_ d1-distinguish s from t
            int sf = sample_.matrix().entry(s, f);
            int tf = sample_.matrix().entry(t, f);
            if ((sf == 0) != (tf == 0)) continue; // f d1-distinguishes s from t

            int sprime_f = sample_.matrix().entry(sprime, f);
            int tprime_f = sample_.matrix().entry(tprime, f);

            int type_s = sprime_f - sf; // <0 if DEC, =0 if unaffected, >0 if INC
            int type_t = tprime_f - tf; // <0 if DEC, =0 if unaffected, >0 if INC

            // Get the sign
            type_s = (type_s > 0) ? 1 : ((type_s < 0) ? -1 : 0);
            type_t = (type_t > 0) ? 1 : ((type_t < 0) ? -1 : 0);

            if(type_s != type_t) {
                features.push_back(f);
            }
        }
        return features;
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

    std::pair<bool, CNFWriter> write_maxsat(std::ostream &os, bool verbose) {
        unsigned num_states = ns_;
        unsigned num_features = nf_;

        CNFWriter writer(os);

        /////// Create the CNF variables ///////

        // Selected(f) for each feature f
        std::vector<cnfvar_t> var_selected;
        for (unsigned f = 0; f < nf_; ++f) {
//            std::cout << "#" << f << ": " << sample_.matrix().feature_name(f) << std::endl;
            var_selected.push_back(writer.variable());
        }

        // D2(s, s', t, t')
        d2map_t d2vars;
        using bridge_clause_t = std::tuple<unsigned, unsigned, unsigned>;
        std::unordered_set<bridge_clause_t, boost::hash<bridge_clause_t>> bridge_clauses;
        for (const transition_t& tx1:sound_transitions()) {
            unsigned s = tx1.first, sprime = tx1.second;

            for (unsigned t = 0; t < ns_; ++t) {
                // If s and t have to be distinguished (because only one of them is a goal),
                // then no need to consider bridge clauses for them, as they will be trivially true
                if (is_goal(s) != is_goal(t)) continue; // NOTE: This wasn't being checked in the Python version
                // Bridge constraints not reflexive
                if (s == t) continue;


                for (unsigned tprime:successors(t)) {
                    if (!(is_sound_transition(t, tprime) && t < s)) {
                        // Symmetry-breaking: no need to define two D2 variables for permutations of s, t
                        auto res = d2vars.emplace(d2idx(s, sprime, t, tprime), writer.variable());
                        assert(res.second); // Make sure D2 was not already declared
                    }

                    bridge_clauses.emplace(s, sprime, t);
                }
            }
        }
        // No more variables will be created. Print total count.
        std::cout << "A total of " << writer.nvars() << " variables were created" << std::endl;


        /////// Create the CNF constraints ///////

        std::cout << "Generating weighted selected constraints for " << var_selected.size() << " features" << std::endl;
        for (unsigned f = 0; f < nf_; ++f) {
            unsigned w = feature_weight(f);
            writer.print_clause({CNFWriter::literal(var_selected[f], false)}, w);
            n_selected_clauses += 1;
        }

        std::cout << "Generating bridge constraints for " << bridge_clauses.size() << " triplets" << std::endl;
        for (const auto& bc:bridge_clauses) {
            unsigned s = std::get<0>(bc), sprime = std::get<1>(bc), t = std::get<2>(bc);

            const auto& d1feats = d1_distinguishing_features(s, t);

            // Start with OR_i Selected(f_i)
            cnfclause_t clause;
            for (unsigned f:d1feats) {
                clause.push_back(CNFWriter::literal(var_selected.at(f), true));
            }

            // And add a literal -D2(s,s',t,t') for each child t' of t
            for (unsigned tprime:successors(t)) {
                auto d2_var = d2vars.at(d2idx(s, sprime, t, tprime));
                clause.push_back(CNFWriter::literal(d2_var, false));
            }

            writer.print_clause(clause);
            n_bridge_clauses += 1;
        }

        // Force D1(s1, s2) to be true if exactly one of the two states is a goal state
        std::cout << "Generating goal constraints for " << goals_.size() * nongoals_.size() << " state pairs" << std::endl;
        for (unsigned s:goals_) {
            for (unsigned t:nongoals_) {
                const auto& d1feats = d1_distinguishing_features(s, t);
                if (d1feats.empty()) {
                    undist_goal_warning(s, t);
                    return {true, writer};
                }

                cnfclause_t clause;
                for (unsigned f:d1feats) {
                    clause.push_back(CNFWriter::literal(var_selected.at(f), true));
                }

                writer.print_clause(clause);
                n_goal_clauses += 1;
            }
        }

        std::cout << "Generating D2 constraints for " << d2vars.size() << " variables" << std::endl;
        for (const auto& d2elem:d2vars) {
            cnflit_t d2lit = CNFWriter::literal(d2elem.second, true);

            // D2(s0,s1,t0,t2) <-- OR_f selected(f), where f ranges over features that d2-distinguish the transition
            // but do _not_ d1-distinguish the two states at the origin of each transition.
            for (unsigned f:d2_distinguishing_features(d2elem.first)) {
                writer.print_clause({d2lit, CNFWriter::literal(var_selected.at(f), false)});
                n_d2_clauses += 1;
            }
        }

        return {false, writer};
    }
};