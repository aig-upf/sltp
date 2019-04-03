
#pragma once

#include <algorithm>

namespace d2tree {

using sorted_nodes_t = std::vector<std::tuple<unsigned, unsigned, unsigned>>;


class D2TreeBuilder {
protected:
    CNFWriter& writer_;
    unsigned nleaves_;
    std::vector<cnfvar_t> d2treevars_;

    const std::vector<std::vector<feature_t>>& leaf_features_;
    std::vector<std::vector<feature_t>> internal_features_;

    unsigned min_intersection_size_;

    //! tree structure
    std::unordered_map<unsigned, unsigned> parent_;
    std::unordered_map<unsigned, std::vector<unsigned>> children_;

    // Each d2-pair can have a different canonical pair, if the set of features that d2-distinguish them are the same
    std::vector<unsigned> canonicals_;


public:

    D2TreeBuilder(
            CNFWriter& writer,
            unsigned nleaves,
            const std::vector<std::vector<feature_t>>& leaf_features,
            unsigned min_intersection_size) :

            writer_(writer),

            nleaves_(nleaves),
            d2treevars_(nleaves, (cnfvar_t) -1), // initialize the vector with fake variable IDs

            leaf_features_(leaf_features),

            // Initialize with n empty elements that will act as offset for the deindexing of the vector
            internal_features_(nleaves),

            min_intersection_size_(min_intersection_size),
            parent_(),
            canonicals_()
    {
        for (unsigned i = 0; i < nleaves; ++i)  canonicals_.push_back(i);
    }

    const std::vector<feature_t>& get_distinguishing_features(unsigned i) const {
        return (i < leaf_features_.size()) ? leaf_features_.at(i) : internal_features_.at(i);
    }

    sorted_nodes_t sort_pairs(unsigned start, unsigned end) {
        std::vector<std::tuple<unsigned, unsigned, unsigned>> pairs;
        assert(canonicals_.size() == end);
        std::cout << "\tSorting pairs from " << start << " to " << end << "...";

        for (unsigned i = start; i < end; ++i) {
            if (canonicals_[i] != i) continue; // i is already subsumed by a previous variable

            for (unsigned  j = i+1; j < end; ++j) {
                if (canonicals_[j] != j) continue;

                auto& d_f1 = get_distinguishing_features(i);
                auto& d_f2 = get_distinguishing_features(j);
                std::vector<unsigned> res;
                std::set_intersection(d_f1.begin(), d_f1.end(), d_f2.begin(), d_f2.end(), std::back_inserter(res));
                auto intersection_size = res.size();
                if (intersection_size > min_intersection_size_) {
                    pairs.emplace_back(intersection_size, i, j);
                }
            }
        }
        std::sort(pairs.begin(), pairs.end(), std::greater<>());
        std::cout << " " << pairs.size() << " suitable pairs found" << std::endl;
        return pairs;
    }

    void mark_duplicities(unsigned start, unsigned end) {
        assert(canonicals_.size() == end);
        for (unsigned i = start; i < end; ++i) assert(canonicals_[i] == i);
        std::cout << "\tSearching duplicate D(n) sets from " << start << " to " << end << "...";
        unsigned num_found = 0;

        for (unsigned i = start; i < end; ++i) {
            if (canonicals_[i] != i) continue; // i is already subsumed by a previous variable

            for (unsigned j = i+1; j < end; ++j) {
                if (canonicals_[j] != j) continue;
                if (get_distinguishing_features(i) == get_distinguishing_features(j)) {
                    canonicals_[j] = i;
                    ++num_found;
                }
            }
        }
        std::cout << " " << num_found << " duplicities found" << std::endl;
    }

    //! Flatten the tree of canonical relations
    void flatten_tree() {
        for (bool fixpoint = false; !fixpoint;) {
            fixpoint = true;
            for (unsigned i = 0; i < canonicals_.size(); ++i) {
                auto can = canonicals_[i];
                if (can != i and can != canonicals_[can]) {
                    canonicals_[i] = canonicals_[can];
                    fixpoint = false;
                }
            }
        }
    }

    void register_parent(unsigned parent, unsigned child) {
        assert (parent_.find(child) == parent_.end());
        parent_[child] = parent;
        children_[parent].push_back(child);
    }

    void build() {
        std::cout << "D2-tree - Building tree..." << std::endl;
        unsigned n_subsuming_merges = 0, n_duplicates = 0;

        unsigned start = 0, end = (unsigned) d2treevars_.size();


        // Mark duplicities within the set of original D2 variables
        mark_duplicities(start, end);

        for (unsigned layer = 0; ; ++layer) {
            std::cout << "D2-tree - Layer #" << layer << " [" << start << ", " << end << "]:" << std::endl;


            // Compute pairs of indexes <i, j> sorted by (largest) size of the intersection S(i) \cap S(j) first
            sorted_nodes_t pairs = sort_pairs(start, end);

            for (const auto& elem:pairs) {
                unsigned n1 = std::get<1>(elem), n2 = std::get<2>(elem);
                if (parent_.find(n1) != parent_.end() || parent_.find(n2) != parent_.end())
                    continue; // Some node in the pair already has one parent, we ignore the pair

                // Compute the intersection of S(n1) and S(n2)
                auto& d_f1 = get_distinguishing_features(n1);
                auto& d_f2 = get_distinguishing_features(n2);
                std::vector<feature_t> inters;
                std::set_intersection(d_f1.begin(), d_f1.end(), d_f2.begin(), d_f2.end(), std::back_inserter(inters));
                auto isize = inters.size();

                unsigned merge_id = (unsigned)-1;

                bool is_duplicate = false;
                for (unsigned i = end; i < d2treevars_.size(); ++i) {
                    if (get_distinguishing_features(i) == inters) {
                        // The intersection of n1 and n2 is equal to the intersection of two other nodes in the same
                        // layer, hence we create one single parent with the four nodes as children
                        is_duplicate = true;
                        register_parent(i, n1);
                        register_parent(i, n2);
                        merge_id = i; // i is considered the result of the merge
                        ++n_duplicates;
                        break;
                    }
                }


                // If the merge of n1 and n2 is not a duplicate, we create a new node to represent it in the tree
                if (!is_duplicate) {
                    merge_id = (unsigned) d2treevars_.size();
                    d2treevars_.push_back((cnfvar_t) -1);
                    register_parent(merge_id, n1);
                    register_parent(merge_id, n2);

                    assert (canonicals_.size() == merge_id);
                    canonicals_.push_back(merge_id);
                    internal_features_.push_back(std::move(inters));
                }


                // If one node was a strict subset of the other, e.g. n1 \subset n2, then n1 is still the canonical
                // node for the merge, as they will both be identical. We detect this by checking size of intersection
                if (isize == d_f1.size() || isize == d_f2.size()) {
                    assert(isize != d_f1.size() || isize != d_f2.size()); // cannot be both equal, otherwise n1=n2

                    if (isize == d_f1.size()) {
                        canonicals_[merge_id] = canonicals_[n1];
                    } else {
                        canonicals_[merge_id] = canonicals_[n2];
                    }
                    ++n_subsuming_merges;
                }
            }

            // Get some stats for reporting
            auto nmerges = report_progress(start, end);
            if (!nmerges) break;

            // Update the layer bounds
            start = end;
            end = (unsigned) d2treevars_.size();
        }

        assert(canonicals_.size() == d2treevars_.size());

//        flatten_tree();

        report_final_stats(start, end, n_subsuming_merges, n_duplicates);
//        debug();

//        std::cout << "Nodes: " << std::endl;
//        for (unsigned i = 0; i < nd2vars; ++i) std::cout << i << ": " << d2_distinguishing_features(i).size() << std::endl;

//        std::cout << "Pairs: " << std::endl;
//        while (!pairs.empty()) {
//            const auto& x = pairs.top();
//            std::cout << "(" << std::get<1>(x) << ", " << std::get<2>(x) << ")[" << std::get<0>(x) << "], ";
//            if (pairs.size() % 10 == 0) std::cout << std::endl;
//            pairs.pop();
//        }
    }

    std::vector<cnfvar_t> generate_variables() {
        std::vector<cnfvar_t> leaf_variables;
        unsigned created = 0;
        // Create one single CNF variable for each canonical D2 variable
        for (unsigned i = 0; i < d2treevars_.size(); ++i) {
            auto can = canonicals_[i];
            assert(can <= i);
            assert (can == i || canonicals_[can] == can);  // the tree must have one single level
            assert(d2treevars_[i] == (cnfvar_t) -1);

            if (can==i) {
                d2treevars_[i] = writer_.variable();
                created++;
            } else {
                d2treevars_[i] = d2treevars_[can];
            }

            // Store the leaf variables, which are the ones visible to the external script
            if (i < nleaves_) leaf_variables.push_back(d2treevars_[i]);
        }
        std::cout << created << " CNF D2-variables were created" << std::endl;
        return leaf_variables;
    }

    //! Generate the D2 tree constraints. Return the number of generated clauses
    std::size_t generate_constraints(const std::vector<cnfvar_t>& var_selected) {
        std::cout << "Generating D2-Tree constraints for " << d2treevars_.size() << " nodes" << std::endl;
        ulong num_clauses_0 = writer_.nclauses(), nclauses_on_leaf_variables = 0;

        for (unsigned i = 0; i < d2treevars_.size(); ++i) {
            if (canonicals_[i] != i) continue; // The d2-pair is a duplicate of some other pair
            cnflit_t pnlit = CNFWriter::literal(d2treevars_[i], true);
            const auto& feats_i = get_distinguishing_features(i);

            auto it = parent_.find(i);
            if (it == parent_.end()) {
                // An orphan tree node
                // If OR_{f in S(n)} selected(f) --> p(n);  for nodes n without parents

                for (feature_t f:feats_i) {
                    writer_.print_clause({pnlit, CNFWriter::literal(var_selected.at(f), false)});

                    if (i < nleaves_) nclauses_on_leaf_variables++;
                }


            } else {
                // A parented node with parent n'
                // If OR_{f in S(n)\S(n')} selected(f) or p(n') --> p(n);   for parent n' of n
                unsigned parent = it->second;

                cnflit_t parentlit = CNFWriter::literal(d2treevars_.at(parent), false);
                writer_.print_clause({pnlit, parentlit});
                if (i < nleaves_) nclauses_on_leaf_variables++;

                const auto& feats_parent = get_distinguishing_features(parent);

                std::vector<unsigned> diff;
                std::set_difference(feats_i.begin(), feats_i.end(), feats_parent.begin(), feats_parent.end(),
                        std::back_inserter(diff));

                for (feature_t f:diff) {
                    writer_.print_clause({pnlit, CNFWriter::literal(var_selected.at(f), false)});
                    if (i < nleaves_) nclauses_on_leaf_variables++;
                }
            }
        }

        auto ngenerated = writer_.nclauses() - num_clauses_0;
        std::cout << "\t" << ngenerated << " D2-Tree constraints generated, of which "
                  << nclauses_on_leaf_variables << " correspond to leaf variables" << std::endl;
        return ngenerated;
    }

    unsigned report_progress(unsigned start, unsigned end) const {
        unsigned unparented = 0, duplicates = 0;
        for (unsigned i = start; i < end; ++i) {
            if (canonicals_[i] == i) {
                if (parent_.find(i) == parent_.end()) unparented++;
            } else {
                duplicates++;
            }
        }
        auto nmerges = d2treevars_.size() - end;

        std::cout << "\tOut of " << end-start << " nodes: " << 2*nmerges << " merged, "<< unparented << " orphaned, "
                  << duplicates << " duplicates." << std::endl;
        return nmerges;
    }

    void report_final_stats(unsigned start, unsigned end, unsigned n_subsuming_merges, unsigned n_duplicates) const {
        unsigned unparented = 0, duplicates = 0;
        for (unsigned i = 0; i < d2treevars_.size(); ++i) {
            if (canonicals_[i] == i) {
                if (parent_.find(i) == parent_.end()) unparented++;
            } else {
                duplicates++;
            }
        }

        auto total_merges = d2treevars_.size() - nleaves_;
        std::cout << "D2-tree - Total node vars: " << d2treevars_.size()
                  << ", of which " << total_merges
                  << " internal and " << nleaves_ << " leaves. "
                  << unparented << " nodes are orphan." << std::endl;
        std::cout << "\t" << n_subsuming_merges << "/" << total_merges << " ("
                  << n_subsuming_merges * 100.0 / (float) total_merges << "%) subsuming merges" << std::endl;
        std::cout << "\t" << n_duplicates << " duplicate merges" << std::endl;
    }

    void debug() const {

        std::cout << "Canonical relation: " << std::endl;
        for (unsigned i = 0; i < canonicals_.size(); ++i) {
            auto can = canonicals_[i];
            if (i != can) {
                std::cout << "\t" << i << " maps to " << can << std::endl;
            }
        }


        std::cout << "CNF variables: " << std::endl;
        for (unsigned i = 0; i < d2treevars_.size(); ++i) {
            std::cout << "\t" << i << " has CNF variable " << d2treevars_[i] << std::endl;
        }
    }



};


} // namespaces