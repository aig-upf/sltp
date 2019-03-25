
#pragma once

#include <cassert>
#include <iostream>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <string>
#include <unordered_set>
#include <vector>

namespace Sample {

    class Matrix {
    public:
        using feature_value_t = uint16_t;

    protected:
        const unsigned num_states_;
        const unsigned num_features_;
        const unsigned num_goals_;

        //! Contains pairs of feature name, feature cost
        std::vector<std::pair<std::string, unsigned>> feature_data_;
        std::unordered_map<std::string, unsigned> feature_name_to_id_;
        std::unordered_set<unsigned> goals_;
        std::vector<std::vector<feature_value_t>> rowdata_;
        std::vector<bool> binary_features_;
        std::vector<bool> numeric_features_;


    public:
        Matrix(unsigned num_states, unsigned num_features, unsigned num_goals)
                : num_states_(num_states),
                  num_features_(num_features),
                  num_goals_(num_goals),
                  binary_features_(num_features_, false),
                  numeric_features_(num_features_, false)
        {}

        virtual ~Matrix() = default;

        unsigned num_states() const { return num_states_; }

        unsigned num_features() const { return num_features_; }

        unsigned num_goals() const { return num_goals_; }

        unsigned feature_cost(unsigned i) const {
            assert(i < num_features_);
            return feature_data_[i].second;
        }

        bool goal(unsigned s) const {
            assert (s < rowdata_.size());
            return goals_.find(s) != goals_.end();
        }

        feature_value_t entry(unsigned s, unsigned f) const {
            return rowdata_[s][f];
        }

        feature_value_t operator()(unsigned s, unsigned f) const {
            return entry(s, f);
        }

        void print(std::ostream &os) const {
            os << "Matrix stats: #states=" << num_states_
               << ", #features=" << num_features_
               << ", #binary-features=" << std::count(binary_features_.begin(), binary_features_.end(), true)
               << ", #numeric-features=" << std::count(numeric_features_.begin(), binary_features_.end(), true)
               << std::endl;
            for (unsigned s = 0; s < num_states_; ++s) {
                os << "state " << s << ":";
                for (unsigned f = 0; f < num_features_; ++f) {
                    feature_value_t value = entry(s, f);
                    if (value > 0)
                        os << " " << f << ":" << value;
                }
                os << std::endl;
            }
        }

        // readers
        void read(std::istream &is) {
            // read features
            for (unsigned i = 0; i < num_features_; ++i) {
                std::string feature;
                is >> feature;
                feature_name_to_id_.emplace(feature, feature_data_.size());
                feature_data_.emplace_back(feature, 0);
            }

            // read feature costs
            for (unsigned i = 0; i < num_features_; ++i) {
                unsigned cost;
                is >> cost;
                assert(cost > 0);
                assert(feature_data_[i].second == 0);
                feature_data_[i].second = cost;
            }

            // read goals
            for (unsigned i = 0; i < num_goals_; ++i) {
                unsigned s;
                is >> s;
                goals_.insert(s);
            }
            assert(goals_.size() == num_goals_);

            // Read the actual feature matrix data
            rowdata_.reserve(num_states_);
            feature_value_t value;
            for (int i = 0; i < num_states_; ++i) {
                unsigned s, nentries;
                is >> s >> nentries;
                assert(i == s);  // Make sure states are listed in increasing order

                std::vector<feature_value_t> data(num_features_, 0);
                for(unsigned j = 0; j < nentries; ++j) {
                    char filler;
                    unsigned f;
                    is >> f >> filler >> value;
                    assert(filler == ':');
                    assert(f < num_features_);
                    assert(value > 0);
                    data[f] = value;
                }
//                for (unsigned f = 0; f < num_features_; ++f) {
//                    is >> value;
//                    data.push_back(value);
//                }
                rowdata_.push_back(std::move(data));
            }

            // Figure out which features are binary, which are numeric
            assert(numeric_features_.size() == num_features_ && binary_features_.size() == num_features_);
            for (unsigned f = 0; f < num_features_; ++f) {
                bool has_value_other_than_0_1 = false;
                for (unsigned s = 0; s < num_states_; ++s) {
                    if (entry(s, f) > 1) {
                        has_value_other_than_0_1 = true;
                        break;
                    }
                }
                if (has_value_other_than_0_1) {
                    numeric_features_[f] = true;
                } else {
                    binary_features_[f] = true;
                }
            }
        }

        static Matrix read_dump(std::istream &is, bool verbose) {
            unsigned num_states, num_features, num_goals;
            is >> num_states >> num_features >> num_goals;
            Matrix matrix(num_states, num_features, num_goals);
            matrix.read(is);
            if (verbose) {
                std::cout << "Matrix::read_dump: "
                          << "#states=" << matrix.num_states()
                          << ", #features=" << matrix.num_features()
                          << ", #goals=" << matrix.num_goals()
                          << std::endl;
            }
            return matrix;
        }
    };

} // namespaces
