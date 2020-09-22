
#include "aaai19_generator.h"
#include "d2tree.h"


namespace sltp::cnf {


TrainingSet AAAI19Generator::preprocess_sample(const TrainingSet& sample, const sltp::cnf::Options& options) {
    auto isomorphisms = compute_redundant_states(sample);

    std::cout << isomorphisms.size() << " / " << sample.matrix().num_states()
              << " were found to be isomorphic and were discarded" << std::endl;

    std::unordered_set<unsigned> nonisomorphic;
    for (unsigned s = 0; s < sample.transitions().num_states(); ++s) {
        if (isomorphisms.find(s) == isomorphisms.end()) nonisomorphic.insert(s);
    }

    // TODO
    throw std::runtime_error("Latest refactorings require reimplementing the resample method");
//    auto resampled = TrainingSet(sample.resample(nonisomorphic));
//    std::cout << "Pruned training sample: " << resampled << std::endl;
//    return resampled;
    return sample;
}

//! Generate and write the actual CNF instance as we go
sltp::cnf::CNFGenerationOutput AAAI19Generator::write(CNFWriter& writer) {

    /////// Create the CNF variables ///////

    // Selected(f) for each feature f
    std::vector<cnfvar_t> var_selected;
    var_selected.reserve(nf_);
    for (unsigned f = 0; f < nf_; ++f) {
//            std::cout << "#" << f << ": " << tr_set_.matrix().feature_name(f) << std::endl;
        var_selected.push_back(writer.variable());
    }

    // D2(s, s', t, t')
    using bridge_clause_t = std::tuple<unsigned, unsigned, unsigned>;
    std::unordered_set<bridge_clause_t, boost::hash<bridge_clause_t>> bridge_clauses;
    for (const auto& tx1:sound_transitions()) {
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
                    d2_features_cache_.push_back(compute_d2_distinguishing_features(tr_set_, s, sprime, t, tprime));
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
                return sltp::cnf::CNFGenerationOutput::UnsatTheory;
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
    for (const auto& tx1:sound_transitions()) {
        unsigned s = tx1.first;
        for (unsigned t:tr_set_.matrix().deadends()) {
            const auto& d1feats = d1_distinguishing_features(s, t);
            if (d1feats.empty()) {
                undist_deadend_warning(s, t);
                return sltp::cnf::CNFGenerationOutput::UnsatTheory;
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

    return sltp::cnf::CNFGenerationOutput::Success;
}

} // namespaces