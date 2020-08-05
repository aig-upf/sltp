
#include "generator.h"
#include "d2tree.h"
#include <algorithm>


//! Generate and write the actual CNF instance as we go
std::pair<bool, CNFWriter> CNFGenerator::write_basic_maxsat(std::ostream &os) {
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
    if (options.use_d2tree()) {
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
    for (auto f:options.enforced_features) {
        writer.print_clause({CNFWriter::literal(var_selected.at(f), true)});
        n_goal_clauses += 1;
    }

    if (options.use_d2tree()) {
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


bool are_transitions_distinguished(int s_f, int sprime_f, int t_f, int tprime_f) {
    if ((s_f == 0) != (t_f == 0)) return true;

    int type_s = sprime_f - s_f; // <0 if DEC, =0 if unaffected, >0 if INC
    int type_t = tprime_f - s_f; // <0 if DEC, =0 if unaffected, >0 if INC

    // Get the sign
    type_s = (type_s > 0) ? 1 : ((type_s < 0) ? -1 : 0);
    type_t = (type_t > 0) ? 1 : ((type_t < 0) ? -1 : 0);

    return type_s != type_t;
}

bool operator<(const transition_pair& x, const transition_pair& y) {
    return std::tie(x.s, x.sprime, x.t, x.tprime) < std::tie(y.s, y.sprime, y.t, y.tprime);
}

boost::container::flat_set<transition_pair> CNFGenerator::
compute_d2_prime(unsigned f) {
    const auto& mat = sample_.matrix();

    boost::container::flat_set<transition_pair> d2prime;

    for (const auto s:all_alive()) {
        for (const auto t:all_alive()) {
            for (unsigned sprime:successors(s)) {
                if (!is_solvable(sprime)) continue;

                for (unsigned tprime:successors(t)) {
                    if (are_transitions_distinguished(
                            mat.entry(s, f), mat.entry(sprime, f), mat.entry(t, f), mat.entry(tprime, f))) {
                        d2prime.emplace(s, sprime, t, tprime);
                    }
                }
            }
        }
    }
//    std::cout << "DT(" << f << ") has " << d2prime.size() << " elements." << std::endl;
    return d2prime;
}

std::vector<bool> CNFGenerator::
check_feature_dominance() {
    const auto& mat = sample_.matrix();
    auto nfeatures = mat.num_features();

    std::vector<bool> dominated(nfeatures, false);
//    return dominated;

    std::cout << "Compute dominance relations... " << std::endl;
    unsigned ndominated = 0;
    for (unsigned f1 = 0; f1 < nfeatures; ++f1) {
        if (dominated[f1]) continue;

        std::cout << "f=" << f1 << std::endl;
        const auto d2_f1 = compute_d2_prime(f1);
        for (unsigned f2 = f1+1; f2 < nfeatures; ++f2) {
            if (dominated[f2]) continue;
            if (feature_weight(f1) > feature_weight(f2)) throw std::runtime_error("Features not ordered by complexity");

            const auto d2_f2 = compute_d2_prime(f2);
            if (d2_f2.size() <= d2_f1.size() && std::includes(d2_f1.begin(), d2_f1.end(), d2_f2.begin(), d2_f2.end())) {
                std::cout << "Feat. " << mat.feature_name(f1) << " dominates " << mat.feature_name(f2) << std::endl;
                ++ndominated;
                dominated[f2] = true;
            } else if (feature_weight(f1) == feature_weight(f2) &&
                    d2_f1.size() <= d2_f2.size() && std::includes(d2_f2.begin(), d2_f2.end(), d2_f1.begin(), d2_f1.end())
            ) {
                std::cout << "Feat. " << mat.feature_name(f1) << " dominates " << mat.feature_name(f2) << std::endl;
                ++ndominated;
                dominated[f1] = true;
            }
        }
    }

    std::cout << "A total of " << ndominated << " are dominated and can be ignored" << std::endl;
    return dominated;
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
std::pair<bool, CNFWriter>
CNFGenerator::write_transition_classification_maxsat(std::ostream &os)
{
    using Wr = CNFWriter;

    auto ignore_features = check_feature_dominance();

    auto varmapstream = get_ofstream(options.workspace + "/varmap.wsat");
    auto selected_map_stream = get_ofstream(options.workspace + "/selecteds.wsat");
    auto allvarsstream = get_ofstream(options.workspace + "/allvars.wsat");

    CNFWriter wr(os, &allvarsstream);

    const unsigned max_d = 10; // TODO Adjust this

    // Keep a map `good_tx_vars` from transitions to SAT variable IDs:
    std::unordered_map<transition_t, cnfvar_t, boost::hash<transition_t>> good_vars;

    // Keep a map `good_and_vleq_vars` from transitions to SAT variable IDs:
    using GV_idx = std::tuple<unsigned, unsigned, unsigned>;
    std::unordered_map<GV_idx, cnfvar_t, boost::hash<GV_idx>> good_and_vleq_vars;

    // Keep a map from pairs (s, d) to SAT variable ID of the variable Vleq(s, d)
    std::vector<std::vector<cnfvar_t>> vleqs(ns_);

    // Keep a map from each feature index to the SAT variable ID of Selected(f)
    std::vector<cnfvar_t> var_selected;

    unsigned n_vleq_vars = 0;
    unsigned n_upper_bound_clauses = 0;
    unsigned n_justification_clauses = 0;
    unsigned n_gv_aux_clauses = 0;


    /////// CNF variables ///////
    for (unsigned f = 0; f < nf_; ++f) {
        if (!ignore_features[f]) {
            auto v = wr.var("Select(" + sample_.matrix().feature_name(f) + ")");
            var_selected.push_back(v);
            selected_map_stream << f << "\t" << v << "\t" << sample_.matrix().feature_name(f) << std::endl;

        } else {
            var_selected.push_back(std::numeric_limits<uint32_t>::max());
        }
    }

    // Create variables Vleq(s, k) that denote that V(s) <= d, for all d in 0..D and all alive state s
    for (const auto s:all_alive()) {

        const auto v_s_0 = wr.var("Vleq(" + std::to_string(s) + ", 0)");
        cnfclause_t at_least_one_v_clause{Wr::lit(v_s_0, true)};

        vleqs[s].reserve(max_d + 1);
        vleqs[s].push_back(v_s_0);
        for (unsigned d = 1; d <= max_d; ++d) {
            const auto v_s_d = wr.var("Vleq(" + std::to_string(s) + ", " + std::to_string(d) + ")");
            vleqs[s].push_back(v_s_d);
            at_least_one_v_clause.push_back(Wr::lit(v_s_d, true));

            // Add clauses (7):  V(s) <= d-1 --> V(s) <= d
            wr.cl({Wr::lit(vleqs[s][d - 1], false), Wr::lit(v_s_d, true)});
        }
        n_leq_clauses += max_d;
        n_vleq_vars += max_d + 1;

        // Add clauses (6)
        wr.cl(at_least_one_v_clause);
        n_leq_clauses += 1;
    }


    for (const auto s:all_alive()) {
        // Good(s, s') for each transition (s, s') such that both s and s' are alive
        // We will in the same loop create all constraints of the form
        //     OR Good(s, s'),
        // where s is alive and s' iterates over all children of s that are solvable
        cnfclause_t clause;
        for (unsigned sprime:successors(s)) {
            if (!is_solvable(sprime)) continue; // alive-to-unsolvable transitions cannot be good

            // Create the Good(s, s') variable
            const auto good_s_sprime = wr.var("Good(" + std::to_string(s) + ", " + std::to_string(sprime) + ")");
            good_vars.insert(std::make_pair(std::make_pair(s, sprime), good_s_sprime));
            //        std::cout << "GOOD(" << tx.first << ", " << tx.second << "): " << vid << std::endl;
            varmapstream << good_s_sprime << " " << s << " " << sprime << std::endl;

            // Push it into the clause
            clause.push_back(Wr::lit(good_s_sprime, true));

            if (is_alive(sprime)) {
                for (unsigned d=0; d < max_d; ++d) {
                    const auto var = wr.var("GV(" + std::to_string(s) + ", " + std::to_string(sprime) + ", " + std::to_string(d) + ")");
                    good_and_vleq_vars.insert({{s, sprime, d}, var});

                    // GV(s, s', d) -> Good(s, s') and Vleq(s', d)
                    wr.cl({Wr::lit(var, false), Wr::lit(good_s_sprime, true)});
                    wr.cl({Wr::lit(var, false), Wr::lit(vleqs[sprime][d], true)});
                    n_gv_aux_clauses += 2;
                }
            }
        }

        // Add clauses (1) for this state
        wr.cl(clause);
        ++n_good_tx_clauses;
    }


    // From this point on, no more variables will be created. Print total count.
    std::cout << "A total of " << wr.nvars() << " variables were created" << std::endl;
    std::cout << "\tSelect(f): " << var_selected.size() << std::endl;
    std::cout << "\tGood(s, s'): " << good_vars.size() << std::endl;
    std::cout << "\tV(s) <= d: " << n_vleq_vars << std::endl;
    std::cout << "\tGV(s, s', d): " << good_and_vleq_vars.size() << std::endl;

    assert(wr.nvars() == var_selected.size() + good_vars.size() + n_vleq_vars + good_and_vleq_vars.size());

    /////// Rest of CNF constraints ///////
    std::cout << "Generating CNF constraints for " << all_alive().size() << " alive states" << std::endl;
    for (const auto s:all_alive()) {

        const auto& succs = successors(s);
        for (unsigned i=0; i < succs.size(); ++i) {
            unsigned sprime = succs[i];
            if (!is_solvable(sprime)) continue;

            if (is_goal(sprime)) {
                for (unsigned d=0; d < max_d; ++d) {
                    // (3) Good(s, s') -> V(s) <= d+1
                    wr.cl({
                        Wr::lit(good_vars.at({s, sprime}), false),
                        Wr::lit(vleqs[s][d+1], true)});
                    ++n_upper_bound_clauses;
                }
            }

            if (is_alive(sprime)) {
                for (unsigned d=0; d < max_d; ++d) {
                    // (2) Good(s, s') and V(s') <= d -> V(s) <= d+1
                    wr.cl({
                        Wr::lit(good_vars.at({s, sprime}), false),
                        Wr::lit(vleqs[sprime][d], false),
                        Wr::lit(vleqs[s][d+1], true)});
                    ++n_upper_bound_clauses;

                    // (3') if Vleq(s,d+1) and -Vleq(s',d) then -Good(s,s')
                    // i.e.: -Good(s, s') OR -Vleq(s, d+1) OR Vleq(s', d)
                    wr.cl({
                        Wr::lit(good_vars.at({s, sprime}), false),
                        Wr::lit(vleqs[s][d+1], false),
                        Wr::lit(vleqs[sprime][d], true)});
                    ++n_upper_bound_clauses;
                }
            }

            // -Good(s, s') or -Good(s, s'') - deterministic policy
            // Uncomment to enforce a deterministic policy
            /*
            for (unsigned j=i+1; j < succs.size(); ++j) {
                unsigned sprimeprime = succs[j];
                wr.cl({Wr::lit(good_vars.at({s, sprime}), false), Wr::lit(good_vars.at({s, sprimeprime}), false)});
                ++n_good_tx_clauses;
            }
            */
        }

        // Clauses (4), (5):
        for (unsigned d=0; d < max_d; ++d) {
            cnfclause_t clause{Wr::lit(vleqs[s][d+1], false)};
            for (unsigned sprime:successors(s)) {
                if (is_goal(sprime)) clause.push_back(Wr::lit(good_vars.at({s, sprime}), true));
                else if (is_alive(sprime)) clause.push_back(Wr::lit(good_and_vleq_vars.at({s, sprime, d}), true));
            }
            wr.cl(clause);
            ++n_justification_clauses;
        }
        wr.cl({Wr::lit(vleqs[s][0], false)});  // s is not a goal
        ++n_justification_clauses;


        // Clauses (8), (9):
        for (const auto t:all_alive()) {
            for (unsigned sprime:successors(s)) {
                if (!is_solvable(sprime)) continue;

                for (unsigned tprime:successors(t)) {
                    cnfclause_t clause{Wr::lit(good_vars.at({s, sprime}), false)};

                    // Either some feature that D1-distinguishes s and t is true
                    for (feature_t f:compute_d1_distinguishing_features(sample_, s, t)) {
                        if (!ignore_features[f]) {
                            clause.push_back(Wr::lit(var_selected.at(f), true));
                        }
                    }

                    // ... or some feature that d2-distinguishes the transitions is true
                    for (feature_t f:compute_d2_distinguishing_features(sample_, s, sprime, t, tprime)) {
                        if (!ignore_features[f]) {
                            clause.push_back(Wr::lit(var_selected.at(f), true));
                        }
                    }

                    if (is_solvable(tprime)) {
                        clause.push_back(Wr::lit(good_vars.at({t, tprime}), true));
                    }

                    wr.cl(clause);
                    n_separation_clauses += 1;
                }
            }
        }
    }


    std::cout << "Generating weighted selected constraints for " << var_selected.size() << " features" << std::endl;
    for (unsigned f = 0; f < nf_; ++f) {
        if (!ignore_features[f]) {
            wr.cl({Wr::lit(var_selected[f], false)}, feature_weight(f));
        }
    }
    n_selected_clauses += nf_;


    // Print a breakdown of the clauses
    std::cout << "A total of " << wr.nclauses() << " clauses were created" << std::endl;
    std::cout << "\t(Weighted) Select(f): " << n_selected_clauses << std::endl;
    std::cout << "\tPolicy completeness [1]: " << n_good_tx_clauses << std::endl;
    std::cout << "\tUpper-bounding V(s) clauses [2,3]: " << n_upper_bound_clauses << std::endl;
    std::cout << "\tClauses justifying V(s) bounds [4]: " << n_justification_clauses << std::endl;
    std::cout << "\tV(s)<=d consistency [5]: " << n_leq_clauses << std::endl;
    std::cout << "\tTransition-separation clauses [6,7]: " << n_separation_clauses << std::endl;
    std::cout << "\tGV(s, s', d) auxiliary clauses: " << n_gv_aux_clauses << std::endl;
    assert(wr.nclauses() == n_selected_clauses + n_good_tx_clauses + n_upper_bound_clauses + n_justification_clauses
    + n_leq_clauses + n_separation_clauses + n_gv_aux_clauses);

    // For goal-identifying features that we want to enforce in the solution, we add a unary clause "selected(f)"
    assert (options.enforced_features.empty()); // ATM haven't really thought whether this feature makes sense for this encoding
//    for (auto f:options.enforced_features) {
//        writer.print_clause({Wr::lit(var_selected.at(f), true)});
//        n_goal_clauses += 1;
//    }

    return {false, wr};
}

std::pair<bool, CNFWriter> CNFGenerator::write_encoding(std::ofstream& os) {
    if (options.use_separation_encoding()) {
        return write_transition_classification_maxsat(os);
    } else {
        return write_basic_maxsat(os);
    }
}
