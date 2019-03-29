
#pragma once

#include <queue>
#include <algorithm>

namespace d2tree {

using sorted_nodes_t = std::priority_queue<std::tuple<unsigned, unsigned, unsigned>>;
template<typename Function>
sorted_nodes_t sort_pairs(std::size_t start, std::size_t end, unsigned min_intersection_size, Function dist_features) {
    std::priority_queue<std::tuple<unsigned, unsigned, unsigned>> pairs;
    for (std::size_t i = start; i < end; ++i) {
//        if (i % 20 == 0) std::cout << "sort_pairs: i: " << i << std::endl;
        for (std::size_t  j = i+1; j < end; ++j) {
            auto& d_f1 = dist_features(i);
            auto& d_f2 = dist_features(j);
            std::vector<unsigned> res;
            std::set_intersection(d_f1.begin(), d_f1.end(), d_f2.begin(), d_f2.end(), std::back_inserter(res));
            auto intersection_size = res.size();
            if (intersection_size > min_intersection_size) {
                pairs.emplace(intersection_size, i, j);
            }
        }
    }
    return pairs;
}


class D2TreeBuilder {
protected:
    CNFWriter& writer_;
    const std::vector<cnfvar_t>& leaf_variables_;
    std::vector<cnfvar_t> d2treevars_;

    const std::vector<std::vector<feature_t>>& leaf_features_;
    std::vector<std::vector<feature_t>> internal_features_;

    unsigned min_intersection_size_;
    std::unordered_map<unsigned, unsigned> parent_; // The parenthood relation

public:

    using feature_getter_t = std::function<const std::vector<unsigned int>& (unsigned)>;

    D2TreeBuilder(
            CNFWriter& writer,
            const std::vector<cnfvar_t>& leaf_variables,
            const std::vector<std::vector<feature_t>>& leaf_features,
            unsigned min_intersection_size) :
            writer_(writer),

            leaf_variables_(leaf_variables),
            d2treevars_(leaf_variables_), // Start with all d2vars.

            leaf_features_(leaf_features),
            internal_features_(leaf_variables_.size()), // Initialize with n empty vectors

            min_intersection_size_(min_intersection_size),
            parent_()
    {}

    const std::vector<feature_t>& get_distinguishing_features(unsigned i) {
        return (i < leaf_features_.size()) ? leaf_features_.at(i) : internal_features_.at(i);
    }

    void build() {
        std::cout << "D2-tree - Building tree..." << std::endl;
        // Initialize with n empty vectors
        feature_getter_t fgetter = [&](unsigned i) -> const std::vector<unsigned int>& {
            return get_distinguishing_features(i);
        };

        d2tree::sorted_nodes_t pairs;
        std::size_t start = 0, end = d2treevars_.size();

        for (unsigned iteration = 0; ; ++iteration) {
            std::cout << "D2-tree - Layer #" << iteration << " [" << start << ", " << end << "]:" << std::flush;

            // Compute pairs of indexes <i, j> sorted by (largest) size of the intersection S(i) \cap S(j) first
            pairs = d2tree::sort_pairs(start, end, min_intersection_size_, fgetter);

            for (; !pairs.empty(); pairs.pop()) {
                const auto& elem = pairs.top();
                unsigned n1 = std::get<1>(elem), n2 = std::get<2>(elem);
                if (parent_.find(n1) != parent_.end() || parent_.find(n2) != parent_.end())
                    continue;

                auto treevar_id = d2treevars_.size();
                d2treevars_.push_back(writer_.variable());
                parent_[n1] = parent_[n2] = (unsigned) treevar_id;


                // Compute the intersection of S(n1) and S(n2)
                auto& d_f1 = fgetter(n1);
                auto& d_f2 = fgetter(n2);
                internal_features_.emplace_back();
                std::set_intersection(d_f1.begin(), d_f1.end(), d_f2.begin(), d_f2.end(), std::back_inserter(internal_features_.back()));
            }

            std::vector<unsigned> unparented;
            for (unsigned i = start; i < end; ++i) if (parent_.find(i) == parent_.end()) unparented.push_back(i);
            auto nmerges = d2treevars_.size() - end;
            std::cout << " " << nmerges << " merges, "<< unparented.size() << " orphaned" << std::endl;
            if (!nmerges) break;

            start = end;
            end = d2treevars_.size();
        }

        std::vector<unsigned> unparented;
        for (unsigned i = 0; i < d2treevars_.size(); ++i) if (parent_.find(i) == parent_.end()) unparented.push_back(i);

        std::cout << "D2-tree - Total node vars: " << d2treevars_.size()
                  << ", of which " << d2treevars_.size() - leaf_variables_.size()
                  << " internal and " << leaf_variables_.size() << " leafs. "
                  << unparented.size() << " nodes are orphan." << std::endl;


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

    //! Generate the D2 tree constraints. Return the number of generated clauses
    std::size_t generate_constraints(const std::vector<cnfvar_t>& var_selected) {
        std::cout << "Generating D2-Tree constraints for " << d2treevars_.size() << " nodes" << std::endl;
        auto num_clauses_0 = writer_.nclauses();
        unsigned nclauses_on_leaf_variables = 0;

        for (unsigned i = 0; i < d2treevars_.size(); ++i) {
            cnflit_t pnlit = CNFWriter::literal(d2treevars_[i], true);

            auto it = parent_.find(i);
            if (it == parent_.end()) {
                // If OR_{f in S(n)} selected(f) --> p(n);  for nodes n without parents

                for (feature_t f:get_distinguishing_features(i)) {
                    writer_.print_clause({pnlit, CNFWriter::literal(var_selected.at(f), false)});

                    if (i < leaf_variables_.size()) nclauses_on_leaf_variables++;
                }


            } else {
                // If OR_{f in S(n)\S(n')} selected(f) or p(n') --> p(n);   for parent n' of n
                unsigned parent = it->second;

                cnflit_t parentlit = CNFWriter::literal(d2treevars_.at(parent), false);
                writer_.print_clause({pnlit, parentlit});
                if (i < leaf_variables_.size()) nclauses_on_leaf_variables++;

                const auto& feats_i = get_distinguishing_features(i);
                const auto& feats_parent = get_distinguishing_features(parent);

                std::vector<unsigned> diff;
                std::set_difference(feats_i.begin(), feats_i.end(), feats_parent.begin(), feats_parent.end(),
                        std::back_inserter(diff));

//                if (i < leaf_variables_.size()) {
//                    std::cout << "D2 variable " << i << " has " << diff.size() << " distinguishing features ("
//                              << feats_i.size() << " - " << feats_parent.size() << ")" << std::endl;
//                }
                 
                for (feature_t f:diff) {
                    writer_.print_clause({pnlit, CNFWriter::literal(var_selected.at(f), false)});
                    if (i < leaf_variables_.size()) nclauses_on_leaf_variables++;
                }
            }
        }

        auto ngenerated = writer_.nclauses() - num_clauses_0;
        std::cout << "A total of " << ngenerated << " D2-Tree constraints were generated" << std::endl;
        std::cout << "Of which " << nclauses_on_leaf_variables << " correspond to leaf variables" << std::endl;
        return ngenerated;
    }



};


} // namespaces