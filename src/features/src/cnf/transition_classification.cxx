
#include "transition_classification.h"

#include <iostream>
#include <vector>
#include <unordered_map>

#include <boost/functional/hash.hpp>



namespace sltp::cnf {


bool are_transitions_distinguished(int s_f, int sprime_f, int t_f, int tprime_f) {
    if ((s_f == 0) != (t_f == 0)) return true;

    int type_s = sprime_f - s_f; // <0 if DEC, =0 if unaffected, >0 if INC
    int type_t = tprime_f - t_f; // <0 if DEC, =0 if unaffected, >0 if INC

    // Get the sign
    type_s = (type_s > 0) ? 1 : ((type_s < 0) ? -1 : 0);
    type_t = (type_t > 0) ? 1 : ((type_t < 0) ? -1 : 0);

    return type_s != type_t;
}


using feature_value_t = Sample::FeatureMatrix::feature_value_t;
sltp::cnf::transition_denotation compute_transition_denotation(feature_value_t s_f, feature_value_t sprime_f) {
    int type_s = (int) sprime_f - (int) s_f; // <0 if DEC, =0 if unaffected, >0 if INC
    int sign = (type_s > 0) ? 1 : ((type_s < 0) ? -1 : 0); // Get the sign
    return sltp::cnf::transition_denotation(s_f, sign);
}


void TransitionClassificationEncoding::compute_equivalence_relations() {
    const auto& mat = sample_.matrix();

    unsigned num_eq_classes = 0;

    for (const auto s:all_alive()) {
        for (unsigned sprime:successors(s)) {
            auto tx = std::make_pair((uint16_t)s, (uint16_t) sprime);
            auto id = (unsigned) transition_ids_inv_.size(); // Assign a sequential ID to the transition

            transition_ids_inv_.push_back(tx);
            transition_ids_.emplace(tx, id);

            // Store the type of the transition
            types_.push_back(is_solvable(sprime) ? transition_type::alive_to_solvable : transition_type::alive_to_dead);

            // Compute the trace of the transition for all features
            transition_trace trace(nf_);
            for (unsigned f = 0; f < nf_; ++f) {
                trace.denotations[f] = compute_transition_denotation(mat.entry(s, f), mat.entry(sprime, f));
            }

            // Check whether some previous transition has the same transition trace
            auto it = from_trace_to_class_repr_.find(trace);
            if (it == from_trace_to_class_repr_.end()) {
                // We have a new equivalence class, to which we assign the ID of the representative transition
                from_transition_to_eq_class_.push_back(id);
                from_trace_to_class_repr_.emplace(trace, id);
                num_eq_classes++;
            } else {
                // We already have other transitions undistinguishable from this one
                assert(it->second < id);
                from_transition_to_eq_class_.push_back(it->second);

//                if (types_[it->second] != types_[id]) {
//                    // We have to non-distinguishable transitions, one from alive to solvable, and one from alive
//                    // to dead; hence, both will need to be labeled as not Good
//                    throw std::runtime_error("Found two non-distinguishable transitions with different types");
//                }
            }
        }
    }

    std::cout << "Number of transitions: " << transition_ids_.size() << std::endl;
    std::cout << "Number of equivalence classes: " << num_eq_classes << std::endl;
}


std::vector<bool> TransitionClassificationEncoding::
check_feature_dominance() {
    const auto& mat = sample_.matrix();

    std::vector<bool> dominated(nf_, false);
    return dominated;

//    std::cout << "Computing sets DT(f)... " << std::endl;
//    std::vector<boost::container::flat_set<transition_pair>> dt(nfeatures);
//    for (unsigned f1 = 0; f1 < nfeatures; ++f1) {
//        dt[f1] = compute_dt(f1);
//    }

    std::cout << "Computing feature-dominance relations... " << std::endl;
    unsigned ndominated = 0;
    for (unsigned f1 = 0; f1 < nf_; ++f1) {
        if (dominated[f1]) continue;

        std::cout << "f=" << f1 << std::endl;
        const auto d2_f1 = compute_dt(f1);
//        const auto& d2_f1 = dt[f1];
        for (unsigned f2 = f1+1; f2 < nf_; ++f2) {
            if (dominated[f2]) continue;
            if (feature_weight(f1) > feature_weight(f2)) throw std::runtime_error("Features not ordered by complexity");

            const auto d2_f2 = compute_dt(f2);
//            const auto& d2_f2 = dt[f2];
            if (d2_f2.size() <= d2_f1.size() && std::includes(d2_f1.begin(), d2_f1.end(), d2_f2.begin(), d2_f2.end())) {
//                std::cout << "Feat. " << mat.feature_name(f1) << " dominates " << mat.feature_name(f2) << std::endl;
                ++ndominated;
                dominated[f2] = true;
            } else if (feature_weight(f1) == feature_weight(f2) &&
                       d2_f1.size() <= d2_f2.size() && std::includes(d2_f2.begin(), d2_f2.end(), d2_f1.begin(), d2_f1.end())
                    ) {
//                std::cout << "Feat. " << mat.feature_name(f1) << " dominates " << mat.feature_name(f2) << std::endl;
                ++ndominated;
                dominated[f1] = true;
            }
        }
    }

    std::cout << "A total of " << ndominated << " features are dominated by some less complex feature and can be ignored" << std::endl;
    return dominated;
}

boost::container::flat_set<transition_pair> TransitionClassificationEncoding::
compute_dt(unsigned f) {
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



std::pair<bool, CNFWriter> TransitionClassificationEncoding::write(std::ostream &os)
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
    unsigned n_leq_clauses = 0;
    unsigned n_good_tx_clauses = 0;
    unsigned n_selected_clauses = 0;
    unsigned n_separation_clauses = 0;


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





} // namespaces