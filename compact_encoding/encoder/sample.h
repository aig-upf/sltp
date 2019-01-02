#ifndef SAMPLE_H
#define SAMPLE_H

#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Sample {

class Matrix {
  protected:
    const int num_states_;
    const int num_features_;
    const int num_goals_;
    const int last_numerical_feature_;
    const int first_boolean_feature_;
    int num_non_zero_entries_;

    std::vector<std::string> features_;
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
        num_non_zero_entries_(0) {
    }
    virtual ~Matrix() {
        for( size_t i = 0; i < rows_.size(); ++i )
            delete[] rows_[i].second;
    }

    int num_states() const {
        return num_states_;
    }
    int num_features() const {
        return num_features_;
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
    int num_non_zero_entries() const {
        return num_non_zero_entries_;
    }
    const std::string& feature(int i) const {
        assert((0 <= i) && (i < num_features_));
        return features_[i];
    }
    int feature(const std::string &f) const {
        std::map<std::string, int>::const_iterator it = feature_map_.find(f);
        return it != feature_map_.end() ? it->second : -1;
    }
    bool numerical(int i) const {
        return i < last_numerical_feature_;
    }

    // accessors
    bool goal(int s) const {
        assert((0 <= s) && (s < rows_.size()));
        return goals_.find(s) != goals_.end();
    }
    int binary_search(int s, int f) const {
        assert((0 <= s) && (s < rows_.size()));
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
        assert((0 <= s) && (s < rows_.size()));
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
            features_.push_back(feature);
        }

        // read goals
        is >> n;
        assert(n == num_goals_);
        for( int i = 0; i < num_goals_; ++i ) {
            int s;
            is >> s;
            goals_.insert(s);
        }
        assert(goals_.size() == num_goals_);

        // read rows
        rows_ = std::vector<std::pair<int, std::pair<int, int>*> >(num_states_, std::make_pair(0, nullptr));
        for( int i = 0; i < num_states_; ++i ) {
            int s, count, f, value;
            is >> s >> count;
            rows_[s] = std::make_pair(count, new std::pair<int, int>[count]);
            for( int j = 0; j < count; ++j ) {
                char filler;
                is >> f >> filler >> value;
                assert(filler == ':');
                assert((0 <= f) && (f < num_features_));
                assert(value > 0);
                assert((j == 0) || (rows_[s].second[j - 1].first < f)); // verify unique and ordered entries
                rows_[s].second[j] = std::make_pair(f, value);
                ++num_non_zero_entries_;
            }
        }
    }
    static const Matrix* read_dump(std::istream &is) {
        int num_states, num_features, num_goals, last_numerical_feature, first_boolean_feature;
        is >> num_states >> num_features >> num_goals >> last_numerical_feature >> first_boolean_feature;
        Matrix *M = new Matrix(num_states, num_features, num_goals, last_numerical_feature, first_boolean_feature);
        M->read(is);
        std::cout << "Matrix::read_dump: "
                  << "#states=" << M->num_states()
                  << ", #features=" << M->num_features()
                  << ", #goals=" << M->num_goals()
                  << ", #last-numerical-feature=" << M->last_numerical_feature()
                  << ", #first-boolean-feature=" << M->first_boolean_feature()
                  << ", #non-zero-entries=" << M->num_non_zero_entries()
                  << std::endl;
        return M;
    }
};

class Transitions {
  protected:
    const int num_states_;
    const int num_transitions_;
    std::vector<std::pair<int, int*> > tr_;
    std::vector<int> offset_;

  public:
    Transitions(int num_states, int num_transitions)
      : num_states_(num_states), num_transitions_(num_transitions) {
        tr_ = std::vector<std::pair<int, int*> >(num_states_, std::make_pair(0, nullptr));
        offset_ = std::vector<int>(num_states_, 0);
    }
    virtual ~Transitions() {
        for( int i = 0; i < tr_.size(); ++i )
            delete[] tr_[i].second;
    }

    int num_states() const {
        return num_states_;
    }
    int num_transitions() const {
        return num_transitions_;
    }
    int num_transitions(int s) const {
        assert((0 <= s) && (s < num_states_));
        return tr_[s].first;
    }
    int transition_state(int s, int l) const {
        assert((0 <= s) && (s < num_states_));
        assert((0 <= l) && (l < tr_[s].first));
        return tr_[s].second[l];
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

    // output
    void dump(std::ostream &os) const {
        os << num_states_ << " " << num_transitions_ << std::endl;

        int count = 0;
        for( size_t s = 0; s < num_states_; ++s )
            count += tr_[s].first > 0 ? 1 : 0;
        os << count << std::endl;

        for( size_t s = 0; s < num_states_; ++s ) {
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
        for( size_t s = 0; s < num_states_; ++s ) {
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
        int num_records = 0, last_src = -1;
        is >> num_records;
        for( int i = 0; i < num_records; ++i ) {
            int src, count, dst;
            is >> src >> count;

            while( last_src < src ) {
                offset_[last_src + 1] = last_src == -1 ? 0 : offset_[last_src];
                ++last_src;
            }

            tr_[src] = std::make_pair(count, new int[count]);
            for( int j = 0; j < count; ++j ) {
                is >> dst;
                assert((0 <= dst) && (dst < num_states_));
                assert((j == 0) || (tr_[src].second[j - 1] < dst)); // verify unique and ordered transitions
                tr_[src].second[j] = dst;
            }

            if( src + 1 < num_states_ ) offset_[src + 1] = offset_[src] + count;
            assert(last_src == src);
            ++last_src;
        }

        // fill remaining offsets
        while( last_src < num_states_ ) {
            offset_[last_src + 1] = last_src == -1 ? 0 : offset_[last_src];
            ++last_src;
        }

        // check all transitions are accounted for
        assert(offset_[num_states_ - 1] + tr_[num_states_ - 1].first == num_transitions_);
    }
    static const Transitions* read_dump(std::istream &is) {
        int num_states, num_transitions;
        is >> num_states >> num_transitions;
        Transitions *T = new Transitions(num_states, num_transitions);
        T->read(is);
        std::cout << "Transitions::read_dump: #states=" << T->num_states() << ", #transitions=" << T->num_transitions() << std::endl;
        return T;
    }
};

class Sample {
  protected:
    const Matrix *matrix_;
    const Transitions *transitions_;

  public:
    Sample() : matrix_(nullptr), transitions_(nullptr) { }
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

    void read_matrix(std::istream &is) {
        delete matrix_;
        matrix_ = Matrix::read_dump(is);
    }
    void read_transitions(std::istream &is) {
        delete transitions_;
        transitions_ = Transitions::read_dump(is);
    }
    void read(const std::string &matrix_filename, const std::string &transitions_filename) {
        std::cout << "reading '" << matrix_filename << "' ... " << std::flush;
        std::ifstream ifs_matrix(matrix_filename.c_str());
        if( !ifs_matrix.fail() ) {
            read_matrix(ifs_matrix);
            ifs_matrix.close();
        } else {
            throw std::runtime_error("error: opening file '" + matrix_filename + "'");
        }

        std::cout << "reading '" << transitions_filename << "' ... " << std::flush;
        std::ifstream ifs_transitions(transitions_filename.c_str());
        if( !ifs_transitions.fail() ) {
            read_transitions(ifs_transitions);
            ifs_transitions.close();
        } else {
            throw std::runtime_error("error: opening file '" + transitions_filename + "'");
        }
    }
};

}; // end of namespace Sample

#endif

