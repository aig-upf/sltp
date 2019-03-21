
#pragma once

#include <cassert>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>

#include "utils.h"

namespace Sample {

class Matrix {
  protected:
    const int num_states_;
    const int num_features_;
    const int num_goals_;
    const int last_numerical_feature_;
    const int first_boolean_feature_;
    const int num_dummy_features_;
    int num_non_zero_entries_;

    std::vector<std::pair<std::string, int> > features_;
    std::map<std::string, int> feature_map_;
    std::set<int> goals_;
    std::vector<std::pair<int, std::pair<int, int>*> > rows_;

  public:
    Matrix(int num_states, int num_features, int num_goals, int last_numerical_feature, int first_boolean_feature)
      : num_states_(num_states),
        num_features_(num_features),
        num_goals_(num_goals),
        last_numerical_feature_(last_numerical_feature),
        first_boolean_feature_(first_boolean_feature),
        num_dummy_features_(first_boolean_feature_ - last_numerical_feature_),
        num_non_zero_entries_(0) {
    }
    virtual ~Matrix() {
        for( size_t i = 0; i < rows_.size(); ++i )
            delete[] rows_[i].second;
    }

    int num_states() const {
        return num_states_;
    }
    int num_features(bool remove_dummy_features = false) const {
        return remove_dummy_features ? num_features_ - num_dummy_features_ : num_features_;
    }
    int num_goals() const {
        return num_goals_;
    }
    int last_numerical_feature() const {
        return last_numerical_feature_;
    }
    int first_boolean_feature() const {
        return first_boolean_feature_;
    }
    int num_dummy_features() const {
        return num_dummy_features_;
    }
    int num_non_zero_entries() const {
        return num_non_zero_entries_;
    }

    const std::string& feature(int i) const {
        assert((0 <= i) && (i < num_features_));
        return features_[i].first;
    }
    int feature_cost(int i) const {
        assert((0 <= i) && (i < num_features_));
        return features_[i].second;
    }
    int feature(const std::string &f) const {
        std::map<std::string, int>::const_iterator it = feature_map_.find(f);
        return it != feature_map_.end() ? it->second : -1;
    }
    bool is_dummy(int i) const {
        return (last_numerical_feature_ <= i) && (i < first_boolean_feature_);
    }
    bool is_numerical(int i) const {
        return i < last_numerical_feature_;
    }

    // remap features when removing dummies
    int remap_feature(int f) const {
        assert((f < last_numerical_feature_) || (first_boolean_feature_ <= f));
        return f < last_numerical_feature_ ? f : f - num_dummy_features_;
    }

    // accessors
    bool goal(int s) const {
        assert((0 <= s) && (s < int(rows_.size())));
        return goals_.find(s) != goals_.end();
    }

    int binary_search(int s, int f) const {
        assert((0 <= s) && (s < int(rows_.size())));
        int start = 0;
        int end = rows_[s].first - 1;
        while( start < end ) {
            int mid = (start + end) / 2;
            if( rows_[s].second[mid].first == f ) {
                assert(rows_[s].second[mid].second > 0);
                return mid;
            } else if( rows_[s].second[mid].first < f ) {
                start = mid + 1;
            } else {
                end = mid - 1;
            }
        }
        return (start == end) && rows_[s].second[start].first == f ? start : -1;
    }
    int entry(int s, int f) const {
        int i = binary_search(s, f);
        return i == -1 ? 0 : rows_[s].second[i].second;
    }
    void get_row(int s, std::vector<int> &row) const {
        assert((0 <= s) && (s < int(rows_.size())));
        row = std::vector<int>(num_features_, 0);
        for( int i = 0; i < rows_[s].first; ++i )
            row[rows_[s].second[i].first] = rows_[s].second[i].second;
    }
    void get_column(int f, std::vector<int> &col) const {
        col = std::vector<int>(num_states_, 0);
        for( int s = 0; s < num_states_; ++s )
            col[f] = entry(s, f);;
    }
    int operator()(int s, int f) const {
        return entry(s, f);
    }

    // output
    void dump(std::ostream &os) const {
        os << num_states_ << " " << num_features_ << " " << last_numerical_feature_ << " " << first_boolean_feature_ << std::endl;
        for( size_t s = 0; s < rows_.size(); ++s ) {
            int count = rows_[s].first;;
            os << count;
            for( int i = 0; i < count; ++i )
                os << " " << rows_[s].second[i].first << ":" << rows_[s].second[i].second;
            os << std::endl;
        }
    }
    void print(std::ostream &os) const {
        os << "Matrix stats: #states=" << num_states_
           << ", #features=" << num_features_
           << ", last-numerical-feature=" << last_numerical_feature_
           << ", first-boolean-feature=" << first_boolean_feature_
           << ", #non-zero-entries=" << num_non_zero_entries_
           << std::endl;
        for( int s = 0; s < num_states_; ++s ) {
            os << "state " << s << ":";
            for( int f = 0; f < num_features_; ++f ) {
                int value = entry(s, f);
                if( value > 0 )
                    os << " " << f << ":" << value;
            }
            os << std::endl;
        }
    }

    // readers
    void read(std::istream &is) {
        // read features
        int n;
        is >> n;
        assert(n == num_features_);
        for( int i = 0; i < num_features_; ++i ) {
            std::string feature;
            is >> feature;
            feature_map_.insert(std::make_pair(feature, features_.size())); // there may be duplicate entries (e.g. 'dummy()')
            features_.emplace_back(feature, -1);
        }

        // read feature costs
        is >> n;
        assert(n == num_features_);
        for( int i = 0; i < num_features_; ++i ) {
            int cost;
            is >> cost;
            assert(is_dummy(i) || (cost > 0));
            assert(features_[i].second == -1);
            features_[i].second = cost;
        }

        // read goals
        is >> n;
        assert(n == num_goals_);
        for( int i = 0; i < num_goals_; ++i ) {
            int s;
            is >> s;
            goals_.insert(s);
        }
        assert(int(goals_.size()) == num_goals_);

        // read rows
        rows_ = std::vector<std::pair<int, std::pair<int, int>*> >(num_states_, std::make_pair(0, nullptr));
        for( int i = 0; i < num_states_; ++i ) {
            int s, count;
            is >> s >> count;
            rows_[s] = std::make_pair(count, new std::pair<int, int>[count]);
            for( int j = 0; j < count; ++j ) {
                size_t value;
                char filler;
                int f;

                is >> f >> filler >> value;
                assert(filler == ':');
                assert((0 <= f) && (f < num_features_));
                assert((j == 0) || (rows_[s].second[j - 1].first < f)); // verify unique and ordered entries
                rows_[s].second[j] = std::make_pair(f, value >= size_t(std::numeric_limits<int>::max()) ? std::numeric_limits<int>::max() : value);
                assert(rows_[s].second[j].second > 0);
                ++num_non_zero_entries_;
            }
        }
    }
    static const Matrix* read_dump(std::istream &is, bool verbose) {
        int num_states, num_features, num_goals, last_numerical_feature, first_boolean_feature;
        is >> num_states >> num_features >> num_goals >> last_numerical_feature >> first_boolean_feature;
        Matrix *M = new Matrix(num_states, num_features, num_goals, last_numerical_feature, first_boolean_feature);
        M->read(is);
        if( verbose ) {
            std::cout << "Matrix::read_dump: "
                      << "#states=" << M->num_states()
                      << ", #features=" << M->num_features()
                      << ", #goals=" << M->num_goals()
                      << ", #last-numerical-feature=" << M->last_numerical_feature()
                      << ", #first-boolean-feature=" << M->first_boolean_feature()
                      << ", #non-zero-entries=" << M->num_non_zero_entries()
                      << std::endl;
        }
        return M;
    }
};

class Transitions {
  protected:
    const int num_states_;
    const int num_transitions_;
    int num_marked_transitions_;
    int max_num_successor_states_;
    std::vector<std::pair<int, int*> > tr_;
    std::vector<int> offset_;

    using hash_for_marked_transitions_t = std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int> > >;
    using hash_for_marked_states_t = std::unordered_set<int>;
    hash_for_marked_transitions_t marked_transitions_;
    hash_for_marked_states_t marked_states_;

  public:
    Transitions(int num_states, int num_transitions, int num_marked_transitions)
      : num_states_(num_states),
        num_transitions_(num_transitions),
        num_marked_transitions_(num_marked_transitions),
        max_num_successor_states_(0) {
        tr_ = std::vector<std::pair<int, int*> >(num_states_, std::make_pair(0, nullptr));
        offset_ = std::vector<int>(num_states_, 0);
    }
    virtual ~Transitions() {
        for( int i = 0; i < int(tr_.size()); ++i )
            delete[] tr_[i].second;
    }

    int num_states() const {
        return num_states_;
    }
    int num_transitions() const {
        return num_transitions_;
    }
    int num_marked_transitions() const {
        return num_marked_transitions_;
    }
    int num_transitions(int s) const {
        assert((0 <= s) && (s < num_states_));
        return tr_[s].first;
    }
    int transition_state(int s, int j) const {
        assert((0 <= s) && (s < num_states_));
        assert((0 <= j) && (j < tr_[s].first));
        return tr_[s].second[j];
    }
    int offset(int s) const {
        assert((0 <= s) && (s < num_states_));
        return offset_[s];
    }

    // accessors
    int binary_search(int s, int t) const {
        assert((0 <= s) && (s < num_states_));
        int start = 0;
        int end = tr_[s].first - 1;
        while( start < end ) {
            int mid = (start + end) / 2;
            if( tr_[s].second[mid] == t ) {
                return mid;
            } else if( tr_[s].second[mid] < t ) {
                start = mid + 1;
            } else {
                end = mid - 1;
            }
        }
        return (start == end) && tr_[s].second[start] == t ? start : -1;
    }
    bool tr(int s, int t) const {
        return binary_search(s, t) != -1;
    }

    int successors(int s, std::vector<int> &succ) const {
        assert((0 <= s) && (s < num_states_));
        int count = tr_[s].first;
        if( count == 0 ) {
            succ.clear();
        } else {
            succ = std::vector<int>(count, 0);
            for( int i = 0; i < count; ++i )
                succ[i] = tr_[s].second[i];
        }
        return count;
    }
    int successors(int s, std::vector<std::pair<int, int> > &succ) const {
        assert((0 <= s) && (s < num_states_));
        int count = tr_[s].first;
        if( count == 0 ) {
            succ.clear();
        } else {
            succ = std::vector<std::pair<int, int> >(count, std::make_pair(0, 0));
            int offset = offset_[s];
            for( int i = 0; i < count; ++i )
                succ[i] = std::make_pair(offset + i, tr_[s].second[i]);
        }
        return count;
    }

    // marked states and transitions
    bool marked(int s) const {
        return marked_states_.find(s) != marked_states_.end();
    }
    bool marked(std::pair<int, int> p) const {
        return marked_transitions_.find(p) != marked_transitions_.end();
    }
    bool marked(int src, int dst) const {
        return marked(std::make_pair(src, dst));
    }
    std::pair<int, int> sample_transition(std::default_random_engine &random_engine, bool verbose) const {
        assert((num_states_ > 0) && (max_num_successor_states_ > 0));
        std::uniform_int_distribution<int> state_sampler(0, num_states_ - 1);
        std::uniform_int_distribution<int> successor_sampler(0, max_num_successor_states_ - 1);
        while( true ) {
            int s = state_sampler(random_engine);
            int j = successor_sampler(random_engine);
            assert((0 <= s) && (s < num_states_));
            assert((0 <= j) && (j < max_num_successor_states_));
            if( j < num_transitions(s) ) {
                return std::make_pair(s, transition_state(s, j));
            } else if( verbose ) {
                std::cout << "sample_transition:"
                          << " state s=" << s << " has no " << j << "th successors (#succ=" << num_transitions(s) << ")"
                          << std::endl;
            }
        }
    }
    void mark_transitions(int n, std::default_random_engine &random_engine, bool verbose) {
        int target_num_marked_transitions = num_marked_transitions_ + n;
        while( num_marked_transitions_ < target_num_marked_transitions ) {
            std::pair<int, int> p = sample_transition(random_engine, verbose);
            if( !marked(p) ) {
                marked_transitions_.emplace(p);
                marked_states_.emplace(p.first);
                marked_states_.emplace(p.second);
                ++num_marked_transitions_;
            }
        }
    }

    // output
    void dump(std::ostream &os) const {
        os << num_states_ << " " << num_transitions_ << std::endl;

        int count = 0;
        for( int s = 0; s < num_states_; ++s )
            count += tr_[s].first > 0 ? 1 : 0;
        os << count << std::endl;

        for( int s = 0; s < num_states_; ++s ) {
            int count = tr_[s].first;
            if( count > 0 ) {
                os << s << " " << count;
                for( int i = 0; i < count; ++i )
                    os << " " << tr_[s].second[i];
                os << std::endl;
            }
        }
    }
    void print(std::ostream &os) const {
        os << "Transitions stats: #states=" << num_states_
           << ", #transitions=" << num_transitions_
           << std::endl;
        for( int s = 0; s < num_states_; ++s ) {
            int count = tr_[s].first;;
            if( count > 0 ) {
                os << "state " << s << ":";
                for( int i = 0; i < count; ++i )
                    os << " " << tr_[s].second[i];
                os << std::endl;
            }
        }
    }

    // readers
    void read(std::istream &is) {
        // read marked transitions
        // (this is somewhat redundant, but it was added afterwards)
        int num_marked_transitions;
        is >> num_marked_transitions;
        assert(num_marked_transitions == num_marked_transitions_);

        std::vector<std::pair<int, int> > marked_transitions;
        marked_transitions.reserve(num_marked_transitions);
        for( int i = 0; i < num_marked_transitions_; ++i ) {
            int src, dst;
            is >> src >> dst;
            assert((0 <= src) && (src < num_states_));
            assert((0 <= dst) && (dst < num_states_));
            marked_transitions.emplace_back(src, dst);
        }
        assert(num_marked_transitions == marked_transitions.size());

        // read number of records in the rest of the file
        int num_records;
        is >> num_records;

        // read transitions
        max_num_successor_states_ = 0;
        for( int i = 0; i < num_records; ++i ) {
            int src, count, dst;
            is >> src >> count;
            assert((0 <= src) && (src < num_states_) && (0 <= count));
            max_num_successor_states_ = std::max(max_num_successor_states_, count);
            if( count > 0 ) {
                tr_[src] = std::make_pair(count, new int[count]);
                bool need_sort = false;
                for( int j = 0; j < count; ++j ) {
                    is >> dst;
                    tr_[src].second[j] = dst;
                    assert((0 <= dst) && (dst < num_states_));
                    need_sort = need_sort || ((j > 0) && (tr_[src].second[j - 1] > dst));
                }
                if( need_sort ) std::stable_sort(&tr_[src].second[0], &tr_[src].second[count]);
            }
        }

        // fill transition offsets
        for( int offset = 0, src = 0; 1 + src < num_states_; ++src ) {
            offset_[src] = offset;
            offset += tr_[src].first;
        }

        // store valid marked transitions
        for( int i = 0; i < int(marked_transitions.size()); ++i ) {
            int src = marked_transitions[i].first;
            int dst = marked_transitions[i].second;
            bool valid = false;
            for( int j = 0; !valid && (j < tr_[src].first); ++j )
                valid = tr_[src].second[j] == dst;
            if( valid ) marked_transitions_.emplace(src, dst);
        }

        // marked states are those appearing in marked transitions
        for( hash_for_marked_transitions_t::const_iterator it = marked_transitions_.begin(); it != marked_transitions_.end(); ++it ) {
            marked_states_.emplace(it->first);
            marked_states_.emplace(it->second);
        }
    }
    static const Transitions* read_dump(std::istream &is, bool verbose) {
        int num_states, num_transitions, num_marked_transitions;
        is >> num_states >> num_transitions >> num_marked_transitions;
        Transitions *T = new Transitions(num_states, num_transitions, num_marked_transitions);
        T->read(is);
        if( verbose ) {
            std::cout << "Transitions::read_dump: #states=" << T->num_states()
                      << ", #transitions=" << T->num_transitions()
                      << ", #marked-states=" << T->marked_states_.size()
                      << ", #marked-transitions=" << T->marked_transitions_.size()
                      << std::endl;
        }
        return T;
    }
};

class Sample {
public:
    const Matrix *matrix_;
    const Transitions *transitions_;

    Sample(const std::string &matrix_filename,
           const std::string &transitions_filename,
           int additional_marked_transitions = 0,
           float additional_marked_transitions_multiple = 0,
           std::default_random_engine *random_engine = nullptr,
           bool verbose = true)
      : matrix_(nullptr), transitions_(nullptr) {
        read(matrix_filename, transitions_filename, verbose);

        // add transitions as marked
        int num_marked_transitions_to_add = additional_marked_transitions;
        num_marked_transitions_to_add += int(additional_marked_transitions_multiple * T().num_marked_transitions());
        if( num_marked_transitions_to_add > 0 ) {
            assert(random_engine != nullptr);
            const_cast<Transitions*>(transitions_)->mark_transitions(num_marked_transitions_to_add, *random_engine, false);
            if( verbose ) std::cout << "sample: #marked-transitions=" << T().num_marked_transitions() << std::endl;
        }
    }
    virtual ~Sample() {
        delete matrix_;
        delete transitions_;
    }

    const Matrix& M() const {
        assert(matrix_ != nullptr);
        return *matrix_;
    }
    const Transitions& T() const {
        assert(transitions_ != nullptr);
        return *transitions_;
    }
    bool expanded(int s) const { // CHECK: HACK: this should be in matrix_, like goal()
        return transitions_->num_transitions(s) > 0;
    }

    // marked states and transitions
    bool marked(int s, bool complete_only_for_marked_transitions) const {
        return !complete_only_for_marked_transitions || transitions_->marked(s);
    }
    bool marked(int src, int dst, bool complete_only_for_marked_transitions) const {
        return !complete_only_for_marked_transitions || transitions_->marked(src, dst);
    }

    void read_matrix(std::istream &is, bool verbose) {
        delete matrix_;
        matrix_ = Matrix::read_dump(is, verbose);
    }
    void read_transitions(std::istream &is, bool verbose) {
        delete transitions_;
        transitions_ = Transitions::read_dump(is, verbose);
    }
    void read(const std::string &matrix_filename, const std::string &transitions_filename, bool verbose) {
        // read matrix
        if( verbose ) {
            std::cout << Utils::blue() << "reading" << Utils::normal()
                      << " '" << matrix_filename << "' ... "
                      << std::endl;
        }
        std::ifstream ifs_matrix(matrix_filename.c_str());
        if( !ifs_matrix.fail() ) {
            read_matrix(ifs_matrix, verbose);
            ifs_matrix.close();
        } else {
            throw std::runtime_error(Utils::error() + "opening file '" + matrix_filename + "'");
        }

        // read transitions
        if( verbose ) {
            std::cout << Utils::blue() << "reading" << Utils::normal()
                      << " '" << transitions_filename << "' ... "
                      << std::endl;
        }
        std::ifstream ifs_transitions(transitions_filename.c_str());
        if( !ifs_transitions.fail() ) {
            read_transitions(ifs_transitions, verbose);
            ifs_transitions.close();
        } else {
            throw std::runtime_error(Utils::error() + "opening file '" + transitions_filename + "'");
        }
    }
};

}; // end of namespace Sample
