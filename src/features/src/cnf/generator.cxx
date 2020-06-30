
#include "generator.h"
#include "d2tree.h"


//! Generate and write the actual CNF instance as we go
std::pair<bool, CNFWriter> CNFGenerator::write_basic_maxsat(std::ostream &os, const std::vector<unsigned>& enforce_features, bool use_d2tree) {
    CNFWriter writer(os);

    /////// Create the CNF variables ///////

    // Selected(f) for each feature f
    std::vector<cnfvar_t> var_selected;
    for (unsigned f = 0; f < nf_; ++f) {
//            std::cout << "#" << f << ": " << sample_.matrix().feature_name(f) << std::endl;
        var_selected.push_back(writer.variable());
    }

    // D2(s, s', t, t')
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

            // If t is a dead-end (and s, because it is a source of a transition, is not), then we'll have to
            // necessarily D1-distinguish s from t. For this reason no bridge clause is necessary.
            // OTOH, if the sample is incomplete, we don't want to place bridge clauses for states t
            // that have NOT been expanded. In either case, the following loop will have zero iterations,
            // but the difference is that we'll enforce later on the D1-distinguishability of dead-ends

            for (unsigned tprime:successors(t)) {
                // Symmetry-breaking: no need to define two D2 variables for permutations of s, t:
                if (!(is_sound_transition(t, tprime) && t < s)) {
                    auto d2_id = d2ids_.size();
                    assert(d2_id == d2_features_cache_.size());

                    auto res = d2ids_.emplace(d2idx(s, sprime, t, tprime), d2_id);
                    assert(res.second); // Make sure this D2 variable was not already declared

                    // Register the set of d2-distinguishing features for this pair of transitions
                    d2_features_cache_.push_back(compute_d2_distinguishing_features(sample_, s, sprime, t, tprime));
                }

                bridge_clauses.emplace(s, sprime, t);
            }
        }
    }

    auto nd2vars = (unsigned) d2ids_.size(); // The number of D2 variables is final at this point

    /////// Build the D2 tree ///////
    d2tree::D2TreeBuilder treebuilder(writer, nd2vars, d2_features_cache_, 5);
    if (use_d2tree) {
        treebuilder.build();
        d2vars_ = treebuilder.generate_variables();
    } else { // Simply initialize all D2 variables with a fresh variable
        for (unsigned i = 0; i < nd2vars; ++i) {
            d2vars_.push_back(writer.variable());
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
            auto d2_var = get_d2var(s, sprime, t, tprime);
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

    // Force D1(s1, s2) to be true if exactly one of the two states is a dead-end
    for (const transition_t& tx1:sound_transitions()) {
        unsigned s = tx1.first;
        for (unsigned t:sample_.matrix().deadends()) {
            const auto& d1feats = d1_distinguishing_features(s, t);
            if (d1feats.empty()) {
                undist_deadend_warning(s, t);
                return {true, writer};
            }

            cnfclause_t clause;
            for (unsigned f:d1feats) {
                clause.push_back(CNFWriter::literal(var_selected.at(f), true));
            }

            writer.print_clause(clause);
            n_deadend_clauses += 1;
        }
    }

    // For goal-identifying features that we want to enforce in the solution, we add a unary clause "selected(f)"
    for (auto f:enforce_features) {
        writer.print_clause({CNFWriter::literal(var_selected.at(f), true)});
        n_goal_clauses += 1;
    }

    if (use_d2tree) {
        n_d2_clauses += treebuilder.generate_constraints(var_selected);
    } else {
        std::cout << "Generating D2 constraints for " << nd2vars << " variables" << std::endl;
        for (unsigned i = 0; i < nd2vars; ++i) {
            cnflit_t d2lit = CNFWriter::literal(d2vars_[i], true);

//            std::cout << "D2 variable " << i << " has " << d2_distinguishing_features(i).size() << " distinguishing features" << std::endl;

            // D2(s0,s1,t0,t2) <-- OR_f selected(f), where f ranges over features that d2-distinguish the transition
            // but do _not_ d1-distinguish the two states at the origin of each transition.
            for (feature_t f:d2_distinguishing_features(i)) {
                writer.print_clause({d2lit, CNFWriter::literal(var_selected.at(f), false)});
                n_d2_clauses += 1;
            }
        }
    }

    return {false, writer};
}


//! Return a sorted vector with those features that d1-distinguish s from t
std::vector<feature_t> compute_d1_distinguishing_features(const Sample::Sample& sample, unsigned s, unsigned t) {
    std::vector<unsigned> features;
    const auto& mat = sample.matrix();
    for(unsigned f = 0; f < mat.num_features(); ++f) {
        auto sf = mat.entry(s, f);
        auto tf = mat.entry(t, f);
        if ((sf == 0) != (tf == 0)) {
            features.push_back(f);
        }
    }
    return features;
}

//! Return a sorted vector with those features that d2-distinguish transition (s, s') from (t, t')
std::vector<feature_t> compute_d2_distinguishing_features(const Sample::Sample& sample,
                                                          unsigned s, unsigned sprime, unsigned t, unsigned tprime) {

    std::vector<unsigned> features;
    const auto& mat = sample.matrix();

    for(unsigned f = 0; f < mat.num_features(); ++f) {
        // Store those features that d2-distinguish (s, s') from (t, t'), but do _not_ d1-distinguish s from t
        int sf = mat.entry(s, f);
        int tf = mat.entry(t, f);
        if ((sf == 0) != (tf == 0)) continue; // f d1-distinguishes s from t

        int sprime_f = mat.entry(sprime, f);
        int tprime_f = mat.entry(tprime, f);

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


//! Check whether t appears isomorphic to s, and in that case, add it to the given list of isomorphisms
void check_isomorphic(const Sample::Sample& sample, unsigned s, unsigned t, isomorphism_t& isomorphisms) {
    // if either s or t are isomorphic of some other state no need to recheck, will be detected in due time
    if (isomorphisms.find(s) != isomorphisms.end() || isomorphisms.find(t) != isomorphisms.end()) return;

    auto distinguishing = compute_d1_distinguishing_features(sample, s, t);
    if (!distinguishing.empty()) return; // s and t are distinguishable, ergo not isomorphic

    if (sample.matrix().goal(s) != sample.matrix().goal(t)) {
        // Only one of the two is a goal: the SAT theory will be unsat
        undist_goal_warning(s, t);
        return;
    }

    if (sample.is_deadend(s) != sample.is_deadend(t)) {
        // Only one of the two is a deadend: the SAT theory will be unsat
        undist_deadend_warning(s, t);
        return;
    }

    if (all_tx_have_analogs(sample, s, t) && all_tx_have_analogs(sample, t, s)) {
        isomorphisms.emplace(t, s);  // t can be pruned in favor of s
    }
}


isomorphism_t compute_redundant_states(const Sample::Sample& sample) {
    isomorphism_t isomorphisms;

    // Collect optimal and non-optimal states
    std::unordered_set<unsigned> optimal_states;
    for (const auto& tx:sample.transitions().marked_transitions()) {
        optimal_states.insert(tx.first);
    }

    std::vector<unsigned> nonoptimal_states;
    for (unsigned s=0; s < sample.transitions().num_states(); ++s) {
        if (optimal_states.find(s) == optimal_states.end()) nonoptimal_states.push_back(s);
    }


    // Check isomorphism between an optimal and a non-optimal state
    for (unsigned s:optimal_states) {
        for (unsigned t:nonoptimal_states) {
            if (s != t) {
                check_isomorphic(sample, s, t, isomorphisms);
            }
        }
    }

    // Check isomorphism between two non-optimal states
    for (unsigned i=0; i < nonoptimal_states.size(); ++i) {
        for (unsigned j=i+1; j < nonoptimal_states.size(); ++j) {
            unsigned s = nonoptimal_states[i], t = nonoptimal_states[j];
            check_isomorphic(sample, s, t, isomorphisms);
        }
    }
    return isomorphisms;
}

bool all_tx_have_analogs(const Sample::Sample& sample, unsigned s, unsigned t) {

    for (unsigned sprime:sample.transitions().successors(s)) {
        bool tx_has_analog = false;
        for (unsigned tprime:sample.transitions().successors(t)) {
            auto d2distinguishing = compute_d2_distinguishing_features(sample, s, sprime, t, tprime);
            if (d2distinguishing.empty()) tx_has_analog = true;
        }

        if (!tx_has_analog) return false;
    }
    return true;
}


//! Generate and write the actual CNF instance as we go
std::pair<bool, CNFWriter> CNFGenerator::write_separation_maxsat(std::ostream &os, const std::vector<unsigned>& enforce_features) {
    CNFWriter writer(os);

    /////// Create the CNF variables ///////
    // Selected(f) for each feature f
    std::vector<cnfvar_t> var_selected;
    for (unsigned f = 0; f < nf_; ++f) {
//            std::cout << "#" << f << ": " << sample_.matrix().feature_name(f) << std::endl;
        var_selected.push_back(writer.variable());
    }

    // Good(s, s') for each optimal transition (s, s')
    // map from transitions to SAT variables:
    std::unordered_map<transition_t, cnfvar_t, boost::hash<transition_t>> var_good_txs;
    for (const transition_t& tx:marked_transitions()) {
        var_good_txs.insert(std::make_pair(tx, writer.variable()));
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

    /*
    TODO: Looks this might not be necessary anymore?
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
    */

    // For goal-identifying features that we want to enforce in the solution, we add a unary clause "selected(f)"
    assert (enforce_features.empty());
    for (auto f:enforce_features) {
        writer.print_clause({CNFWriter::literal(var_selected.at(f), true)});
        n_goal_clauses += 1;
    }


    // Compute which pairs (t1, t2) of transitions need to be distinguished, and for that pair, which
    // features do actually distinguish them
    std::cout << "Generating transition-separation constraints for " << marked_transitions().size()
              << " positive transitions" << std::endl;

    for (const transition_t& tx1:marked_transitions()) {
        unsigned s = tx1.first, sprime = tx1.second;
        const auto& tx1_is_good = var_good_txs.at(tx1);

        // TODO Precompute the non-marked transitions
        for (unsigned t = 0; t < ns_; ++t) {
            for (unsigned tprime:successors(t)) {
                if (!is_marked_transition(t, tprime)) {
                    // We want to distinguish (s, s') from (t, t')
                    cnfclause_t clause;
                    clause.push_back(CNFWriter::literal(tx1_is_good, false));
                    for (feature_t f:compute_d2_distinguishing_features(sample_, s, sprime, t, tprime)) {
                        clause.push_back(CNFWriter::literal(var_selected.at(f), true));
                    }

                    writer.print_clause(clause);
                    n_d2_clauses += 1;
                }
            }
        }
    }

    // For each alive state, make sure that at least one of the outgoing marked transitions is labeled as good.
    for (unsigned s = 0; s < ns_; ++s) {
        if (!is_alive(s)) continue;

        cnfclause_t clause;
        bool at_least_one_successor = false;
        for (unsigned sprime:successors(s)) {
            if (!is_marked_transition(s, sprime)) continue;

            at_least_one_successor = true;
            const auto& tx_is_good = var_good_txs.at(std::make_pair(s, sprime));
            clause.push_back(CNFWriter::literal(tx_is_good, true));
        }
        assert(at_least_one_successor);
        writer.print_clause(clause);
        ++n_good_tx_clauses;
    }


    return {false, writer};
}