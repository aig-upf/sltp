
#include "transition_classification.h"
#include "equivalences.h"

#include <iostream>
#include <vector>
#include <unordered_map>

#include <boost/functional/hash.hpp>



namespace sltp::cnf {


using feature_value_t = Sample::FeatureMatrix::feature_value_t;
sltp::cnf::transition_denotation compute_transition_denotation(feature_value_t s_f, feature_value_t sprime_f) {
    int type_s = (int) sprime_f - (int) s_f; // <0 if DEC, =0 if unaffected, >0 if INC
    int sign = (type_s > 0) ? 1 : ((type_s < 0) ? -1 : 0); // Get the sign
    return sltp::cnf::transition_denotation(bool(s_f > 0), sign);
}


void TransitionClassificationEncoding::compute_equivalence_relations() {
    const auto& mat = sample_.matrix();

    // A mapping from a full transition trace to the ID of the corresponding equivalence class
    std::unordered_map<transition_trace, unsigned> from_trace_to_class_repr;

    for (const auto s:all_alive()) {
        for (unsigned sprime:successors(s)) {
            auto tx = std::make_pair((uint16_t)s, (uint16_t) sprime);
            auto id = (unsigned) transition_ids_inv_.size(); // Assign a sequential ID to the transition

            transition_ids_inv_.push_back(tx);
            transition_ids_.emplace(tx, id);

            // Store the type of the transition
            types_.push_back(is_solvable(sprime) ? transition_type::alive_to_solvable : transition_type::alive_to_dead);

            if (!is_solvable(sprime)) { // An alive-to-dead transition cannot be Good
                necessarily_bad_transitions_.emplace(id);
            }

            if (!options.use_equivalence_classes) {
                // If we don't want to use equivalence classes, we simply create one fictitious equivalence class
                // for each transition, and proceed as usual
                from_transition_to_eq_class_.push_back(id);
                class_representatives_.push_back(id);
                continue;
            }

            // Compute the trace of the transition for all features
            transition_trace trace(nf_);
            for (unsigned f = 0; f < nf_; ++f) {
                trace.denotations[f] = compute_transition_denotation(mat.entry(s, f), mat.entry(sprime, f));
            }

            // Check whether some previous transition has the same transition trace
            auto it = from_trace_to_class_repr.find(trace);
            if (it == from_trace_to_class_repr.end()) {
                // We have a new equivalence class, to which we assign the ID of the representative transition
                from_transition_to_eq_class_.push_back(id);
                from_trace_to_class_repr.emplace(trace, id);
                class_representatives_.push_back(id);
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

    // All transitions that belong to some class where at least one transition must be bad, must be bad
    std::unordered_set<unsigned> necessarily_bad_classes;
    for (const auto id:necessarily_bad_transitions_) {
        necessarily_bad_classes.insert(get_representative_id(id));
    }

    for (unsigned id=0; id < transition_ids_.size(); ++id) {
        auto repr = get_representative_id(id);
        if (necessarily_bad_classes.find(repr) != necessarily_bad_classes.end()) {
            necessarily_bad_transitions_.insert(id);
        }
    }

    // Print some stats
    std::cout << "Number of transitions/equivalence classes: " << transition_ids_.size()
              << "/" << class_representatives_.size() << std::endl;
    std::cout << "Number of necessarily bad transitions/classes: " << necessarily_bad_transitions_.size()
              << "/" << necessarily_bad_classes.size() << std::endl;
}


std::vector<bool> TransitionClassificationEncoding::
check_feature_dominance() {
    const auto& mat = sample_.matrix();

    std::vector<bool> dominated(nf_, false);
    if (!options.use_feature_dominance) return dominated;  // No feature will be considered as dominated

    std::cout << "Computing sets DT(f)... " << std::endl;
    std::vector<std::vector<transition_pair>> dt(nf_);
//    std::vector<boost::container::flat_set<transition_pair>> dt(nf_);
    for (unsigned f1 = 0; f1 < nf_; ++f1) {
        dt[f1] = compute_dt(f1);
    }

    std::cout << "Computing feature-dominance relations... " << std::endl;

    //! The computation and manipulation of the sets DT(f) below, notably the use of std::includes, requires that
    //! the vector of transition IDs is sorted
    assert(std::is_sorted(class_representatives_.begin(), class_representatives_.end()));

    unsigned ndominated = 0;
    for (unsigned f1 = 0; f1 < nf_; ++f1) {
        if (dominated[f1]) continue;

//        std::cout << "f=" << f1 << std::endl;
//        const auto d2_f1 = compute_dt(f1);
        const auto& d2_f1 = dt[f1];
        for (unsigned f2 = f1+1; f2 < nf_; ++f2) {
            if (dominated[f2]) continue;
            if (feature_weight(f1) > feature_weight(f2)) throw std::runtime_error("Features not ordered by complexity");

//            const auto d2_f2 = compute_dt(f2);
            const auto& d2_f2 = dt[f2];
            if (d2_f2.size() <= d2_f1.size() && std::includes(d2_f1.begin(), d2_f1.end(), d2_f2.begin(), d2_f2.end())) {
//                std::cout << "Feat. " << mat.feature_name(f1) << " dominates " << mat.feature_name(f2) << std::endl;
                ++ndominated;
                dominated[f2] = true;

            } else if (feature_weight(f1) == feature_weight(f2) && d2_f1.size() <= d2_f2.size()
                      && std::includes(d2_f2.begin(), d2_f2.end(), d2_f1.begin(), d2_f1.end())) {
//                std::cout << "Feat. " << mat.feature_name(f1) << " dominates " << mat.feature_name(f2) << std::endl;
                ++ndominated;
                dominated[f1] = true;
            }
        }
    }

    std::cout << "A total of " << ndominated << " features are dominated by some less complex feature and can be ignored" << std::endl;
    return dominated;
}

std::vector<transition_pair> TransitionClassificationEncoding::
compute_dt(unsigned f) {
    const auto& mat = sample_.matrix();

//    boost::container::flat_set<transition_pair> dt;
    std::vector<transition_pair> dt;

    for (const auto tx1:class_representatives_) {
        if (is_necessarily_bad(tx1)) continue;
        const auto& tx1pair = get_state_pair(tx1);
        const auto s = tx1pair.first;
        const auto sprime = tx1pair.second;


        for (const auto tx2:class_representatives_) {
            const auto& tx2pair = get_state_pair(tx2);
            const auto t = tx2pair.first;
            const auto tprime = tx2pair.second;

            if (are_transitions_d1d2_distinguished(
                    mat.entry(s, f), mat.entry(sprime, f), mat.entry(t, f), mat.entry(tprime, f))) {
                dt.emplace_back(tx1, tx2);
            }
        }
    }
//    std::cout << "DT(" << f << ") has " << dt.size() << " elements." << std::endl;
    return dt;
}



sltp::cnf::CNFGenerationOutput TransitionClassificationEncoding::write(
        CNFWriter& wr, const std::vector<transition_pair>& transitions_to_distinguish)
{
    using Wr = CNFWriter;

    auto ignore_features = check_feature_dominance();

    auto varmapstream = get_ofstream(options.workspace + "/varmap.wsat");
    auto selected_map_stream = get_ofstream(options.workspace + "/selecteds.wsat");

    const unsigned max_d = compute_D();
    std::cout << "Using an upper bound for V_pi(s) values of " << max_d << std::endl;

    // Keep a map `good_tx_vars` from transition IDs to SAT variable IDs:
    std::unordered_map<unsigned, cnfvar_t> good_vars;

    // Keep a map `good_and_vleq_vars` from transitions to SAT variable IDs:
    using GV_idx = std::tuple<unsigned, unsigned, unsigned>;
    std::unordered_map<GV_idx, cnfvar_t, boost::hash<GV_idx>> good_and_vleq_vars;

    // Keep a map from pairs (s, d) to SAT variable ID of the variable Vleq(s, d)
    std::vector<std::vector<cnfvar_t>> vleqs(ns_);

    // Keep a map from each feature index to the SAT variable ID of Selected(f)
    std::vector<cnfvar_t> var_selected;

    unsigned n_select_vars = 0;
    unsigned n_vleq_vars = 0;
    unsigned n_upper_bound_clauses = 0;
    unsigned n_justification_clauses = 0;
    unsigned n_gv_aux_clauses = 0;
    unsigned n_leq_clauses = 0;
    unsigned n_good_tx_clauses = 0;
    unsigned n_selected_clauses = 0;
    unsigned n_separation_clauses = 0;
    unsigned n_max_v_s_clauses = 0;


    /////// CNF variables ///////
    for (unsigned f = 0; f < nf_; ++f) {
        if (!ignore_features[f]) {
            auto v = wr.var("Select(" + sample_.matrix().feature_name(f) + ")");
            var_selected.push_back(v);
            selected_map_stream << f << "\t" << v << "\t" << sample_.matrix().feature_name(f) << std::endl;
            ++n_select_vars;

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
            auto tx = get_transition_id(s, sprime);
            auto repr = get_representative_id(tx);

            if (is_necessarily_bad(tx)) continue; // This includes  alive-to-dead transitions

            cnfvar_t good_s_sprime = 0;
            if (tx == repr) {
                // Create the Good(s, s') variable
                good_s_sprime = wr.var("Good(" + std::to_string(s) + ", " + std::to_string(sprime) + ")");
                good_vars.emplace(tx, good_s_sprime);
                //        std::cout << "GOOD(" << tx.first << ", " << tx.second << "): " << vid << std::endl;
                varmapstream << good_s_sprime << " " << s << " " << sprime << std::endl;
            } else {
                good_s_sprime = good_vars.at(repr);
            }

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
        if (clause.empty()) {
            throw std::runtime_error("State #" + std::to_string(s) + " has no successor that can be good");
        }
        wr.cl(clause);
        ++n_good_tx_clauses;
    }


    // From this point on, no more variables will be created. Print total count.
    std::cout << "A total of " << wr.nvars() << " variables were created" << std::endl;
    std::cout << "\tSelect(f): " << n_select_vars << std::endl;
    std::cout << "\tGood(s, s'): " << good_vars.size() << std::endl;
    std::cout << "\tV(s) <= d: " << n_vleq_vars << std::endl;
    std::cout << "\tGV(s, s', d): " << good_and_vleq_vars.size() << std::endl;

    assert(wr.nvars() == n_select_vars + good_vars.size() + n_vleq_vars + good_and_vleq_vars.size());

    /////// Rest of CNF constraints ///////
    std::cout << "Generating CNF encoding for " << all_alive().size() << " alive states, "
              <<  transition_ids_.size() << " alive-to-solvable and alive-to-dead transitions and "
              << class_representatives_.size() << " transition equivalence classes" << std::endl;

    for (const auto s:all_alive()) {

        const auto& succs = successors(s);
        for (unsigned i=0; i < succs.size(); ++i) {
            unsigned sprime = succs[i];
            auto tx = get_transition_id(s, sprime);
            if (is_necessarily_bad(tx)) continue; // This includes  alive-to-dead transitions

            auto good_s_sprime = good_vars.at(get_class_representative(s, sprime));

            if (is_goal(sprime)) {
                for (unsigned d=0; d < max_d; ++d) {
                    // (3) Good(s, s') -> V(s) <= d+1
                    wr.cl({
                                  Wr::lit(good_s_sprime, false),
                                  Wr::lit(vleqs[s][d+1], true)});
                    ++n_upper_bound_clauses;
                }
            }

            if (is_alive(sprime)) {
                for (unsigned d=0; d < max_d; ++d) {
                    // (2) Good(s, s') and V(s') <= d -> V(s) <= d+1
                    wr.cl({
                                  Wr::lit(good_s_sprime, false),
                                  Wr::lit(vleqs[sprime][d], false),
                                  Wr::lit(vleqs[s][d+1], true)});
                    ++n_upper_bound_clauses;

                    // (3') if Vleq(s,d+1) and -Vleq(s',d) then -Good(s,s')
                    // i.e.: -Good(s, s') OR -Vleq(s, d+1) OR Vleq(s', d)
                    wr.cl({
                                  Wr::lit(good_s_sprime, false),
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
                wr.cl({Wr::lit(good_s_sprime, false), Wr::lit(good_s_sprime, false)});
                ++n_good_tx_clauses;
            }
            */
        }

        // Clauses (4), (5):
        for (unsigned d=0; d < max_d; ++d) {
            cnfclause_t clause{Wr::lit(vleqs[s][d+1], false)};
            for (unsigned sprime:successors(s)) {
                auto tx = get_transition_id(s, sprime);
                if (is_necessarily_bad(tx)) continue;

                auto good_s_sprime = good_vars.at(get_class_representative(s, sprime));


                if (is_goal(sprime)) clause.push_back(Wr::lit(good_s_sprime, true));
                else if (is_alive(sprime)) clause.push_back(Wr::lit(good_and_vleq_vars.at({s, sprime, d}), true));
            }
            wr.cl(clause);
            ++n_justification_clauses;
        }
        wr.cl({Wr::lit(vleqs[s][0], false)});  // (5) if s is not a goal, then -V(s) <= 0
        ++n_justification_clauses;
    }

    // Clauses (8), (9):
    std::cout << "Posting distinguishability constraints for " << transitions_to_distinguish.size()
              << " pairs of transitions" << std::endl;
    for (const auto& tpair:transitions_to_distinguish) {
        assert (!is_necessarily_bad(tpair.tx1));
        const auto& tx1pair = get_state_pair(tpair.tx1);
        const auto s = tx1pair.first;
        const auto sprime = tx1pair.second;
        const auto& tx2pair = get_state_pair(tpair.tx2);
        const auto t = tx2pair.first;
        const auto tprime = tx2pair.second;

        cnfclause_t clause{Wr::lit(good_vars.at(tpair.tx1), false)};

        // Compute first the Selected(f) terms
        for (feature_t f:compute_d1d2_distinguishing_features(sample_, s, sprime, t, tprime)) {
            if (!ignore_features[f]) {
                clause.push_back(Wr::lit(var_selected.at(f), true));
            }
        }

        if (!is_necessarily_bad(tpair.tx2)) {
            auto good_t_tprime = good_vars.at(tpair.tx2);
            clause.push_back(Wr::lit(good_t_tprime, true));
        }
        wr.cl(clause);
        n_separation_clauses += 1;
    }

    // Clauses (10)
    for (const auto s:all_alive()) {
        auto max_v_s = get_max_v(s);
        assert(max_v_s <= max_d);
        if (max_v_s < 0) throw std::runtime_error("State #" + std::to_string(s) + " has infinite V^* value");
        wr.cl({Wr::lit(vleqs[s][max_v_s], true)});
        n_max_v_s_clauses += 1;
    }


    std::cout << "Posting (weighted) soft constraints for " << var_selected.size() << " features" << std::endl;
    for (unsigned f = 0; f < nf_; ++f) {
        if (!ignore_features[f]) {
            wr.cl({Wr::lit(var_selected[f], false)}, feature_weight(f));
            n_selected_clauses += 1;
        }
    }


    // Print a breakdown of the clauses
    std::cout << "A total of " << wr.nclauses() << " clauses were created" << std::endl;
    std::cout << "\t(Weighted) Select(f): " << n_selected_clauses << std::endl;
    std::cout << "\tPolicy completeness [1]: " << n_good_tx_clauses << std::endl;
    std::cout << "\tUpper-bounding V(s) clauses [2,3,3']: " << n_upper_bound_clauses << std::endl;
    std::cout << "\tClauses justifying V(s) bounds [4,5]: " << n_justification_clauses << std::endl;
    std::cout << "\tV(s)<=d consistency [6,7]: " << n_leq_clauses << std::endl;
    std::cout << "\tTransition-separation clauses [8,9]: " << n_separation_clauses << std::endl;
    std::cout << "\tMax V(s) clauses [10]: " << n_max_v_s_clauses << std::endl;
    std::cout << "\tGV(s, s', d) auxiliary clauses: " << n_gv_aux_clauses << std::endl;
    assert(wr.nclauses() == n_selected_clauses + n_good_tx_clauses + n_upper_bound_clauses + n_justification_clauses
                            + n_leq_clauses + n_separation_clauses + n_gv_aux_clauses + n_max_v_s_clauses);

    // For goal-identifying features that we want to enforce in the solution, we add a unary clause "selected(f)"
    assert (options.enforced_features.empty()); // ATM haven't really thought whether this feature makes sense for this encoding
//    for (auto f:options.enforced_features) {
//        writer.print_clause({Wr::lit(var_selected.at(f), true)});
//        n_goal_clauses += 1;
//    }

    return sltp::cnf::CNFGenerationOutput::Success;
}

CNFGenerationOutput TransitionClassificationEncoding::refine_theory(CNFWriter& wr) {
    std::vector<transition_pair> flaws;
    bool previous_solution = check_existing_solution_for_flaws(flaws);
    if (previous_solution && flaws.empty()) {
        return CNFGenerationOutput::ValidationCorrectNoRefinementNecessary;
    }

    if (options.use_incremental_refinement) {
        return write(wr, compute_transitions_to_distinguish(flaws));
    }

    return write(wr, distinguish_all_transitions());
}

std::vector<transition_pair> TransitionClassificationEncoding::distinguish_all_transitions() const {
    std::vector<transition_pair> transitions_to_distinguish;
    transitions_to_distinguish.reserve(class_representatives_.size() * class_representatives_.size());

    for (const auto tx1:class_representatives_) {
        if (is_necessarily_bad(tx1)) continue;
        for (const auto tx2:class_representatives_) {
            transitions_to_distinguish.emplace_back(tx1, tx2);
        }
    }
    return transitions_to_distinguish;
}

std::vector<transition_pair>
TransitionClassificationEncoding::compute_transitions_to_distinguish(
        const std::vector<transition_pair> &flaws) const {

    if (flaws.size()>1) throw std::runtime_error("Code below needs to be adapted");

    transition_pair onlyflaw = {std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max()};
    if (flaws.size()==1) onlyflaw = flaws[0];

    std::vector<transition_pair> transitions_to_distinguish;
    for (const auto tx1:class_representatives_) {
        if (is_necessarily_bad(tx1)) continue;

        const auto& tx1pair = get_state_pair(tx1);
        const auto s = tx1pair.first;

        std::unordered_set<uint32_t> transitions2;
        for (unsigned sprime:successors(s)) {
            if (sprime == tx1pair.second) continue;
            transitions2.insert(get_class_representative(s, sprime));
        }

        // This will need to be adapted if we want to refine with more than one flaw at a time
        if (onlyflaw.tx1 == tx1) transitions2.insert(get_representative_id(onlyflaw.tx2));

        for (auto tx2:transitions2) {
            transitions_to_distinguish.emplace_back(tx1, tx2);
        }
    }
    return transitions_to_distinguish;
}

    bool TransitionClassificationEncoding::check_existing_solution_for_flaws(
        std::vector<transition_pair>& flaws)  const {
    auto ifs_good_transitions = get_ifstream(options.workspace + "/good_transitions.io");
    auto ifs_good_features = get_ifstream(options.workspace + "/good_features.io");

    std::vector<unsigned> good_features;
    int featureid = -1;
    while (ifs_good_features >> featureid) {
        good_features.push_back(featureid);
    }

    std::vector<unsigned> good_transitions_repr;
    int s = -1, sprime = -1;
    while (ifs_good_transitions >> s >> sprime) {
        good_transitions_repr.emplace_back(get_transition_id(s, sprime));
    }

    ifs_good_transitions.close();
    ifs_good_features.close();

    if (good_features.empty()) return false;

    // Let's exploit the equivalence classes between transitions. The transitions that have been read off as Good from
    // the SAT solution are already class representatives by definition of the SAT theory
    std::vector<unsigned> bad_transitions_repr;
    std::unordered_set<unsigned> good_set(good_transitions_repr.begin(), good_transitions_repr.end());

    for (auto repr:class_representatives_) {
        if (good_set.find(repr) == good_set.end()) {
            bad_transitions_repr.push_back(repr);
        }
    }

    // Let's check whether the policy is indeed able to distinguish between all pairs of good and bad transitions
    for (auto gtx:good_transitions_repr) {
        const auto& tx1pair = get_state_pair(gtx);

        for (auto btx:bad_transitions_repr) {
            const auto& tx2pair = get_state_pair(btx);

            if (!are_transitions_d1d2_distinguishable(tx1pair.first, tx1pair.second, tx2pair.first, tx2pair.second, good_features)) {
                // We found a flaw in the computed policy: Transitions gtx, which is labeled as Good, cannot be
                // distinguished from transition btx, labeled as bad, based only on the selected ("good") features.
                flaws.emplace_back(gtx, btx);
                break;
            }
        }
    }

    std::cout << "Refinement of computed policy found " << flaws.size() << " flaws" << std::endl;
    return true;
}

bool TransitionClassificationEncoding::are_transitions_d1d2_distinguishable(
        uint16_t s, uint16_t sprime, uint16_t t, uint16_t tprime, const std::vector<unsigned>& features) const {
    const auto& mat = sample_.matrix();
    for (unsigned f:features) {
        if (are_transitions_d1d2_distinguished(mat.entry(s, f), mat.entry(sprime, f),
                                               mat.entry(t, f), mat.entry(tprime, f))) {
            return true;
        }
    }
    return false;
}



} // namespaces