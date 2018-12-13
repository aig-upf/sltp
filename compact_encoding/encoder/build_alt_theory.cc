#include <bitset>
#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

// TODO: * action must be applicable at some state
//       * further ordering of actions to reduce more symmetries
//       * separate goal from non-goal states (IMPORTANT)

using namespace std;

// returns k-th bit in integer l
inline bool bit(int k, int l) {
    return l & (1 << k);
}

// returns a string with binary representation of integer k (16 bits)
string bit_string(int k) {
    return bitset<8>(k).to_string();
}

string filename(const string &prefix, int K, int N, const string &suffix, bool opt = true) {
    string fn(prefix);
    if( opt ) {
        fn += string("_") + to_string(K);
        fn += string("_") + to_string(N);
    }
    fn += suffix;
    return fn;
}

class Var {
  protected:
    int index_;
    string str_;

  public:
    Var(int index, const string &str) : index_(index), str_(str) { }
    virtual ~Var() { }
    int index() const { return index_; }
    string str() const { return str_; }
    void print(ostream &os) const {
        os << str();
    }
};

class Literal {
  protected:
    const Var &var_;
    bool sign_;

  public:
    Literal(const Var &var, bool sign)
      : var_(var), sign_(sign) {
    }
    Literal(const Literal &L, bool negate = false)
      : var_(L.var()), sign_(negate ? !L.sign() : L.sign()) {
    }
    virtual ~Literal() { }
    const Var& var() const { return var_; }
    bool sign() const { return sign_; }
    int var_index() const { return var_.index(); }
    int as_int() const { return !sign_ ? var_index() : -var_index(); }
    string as_str() const {
        return !sign_ ? var_.str() : string("-") + var_.str();
    }
    void print(ostream &os) const {
        os << as_str();
    }
};

class Implication {
  protected:
    vector<int> antecedent_; // joined by AND
    vector<int> consequent_; // joined by OR

  public:
    Implication() { }
    ~Implication() { }

    void add_antecedent(int L) { antecedent_.push_back(L); }
    void add_consequent(int L) { consequent_.push_back(L); }

    void dump(ostream &os) const {
        for( size_t i = 0; i < antecedent_.size(); ++i ) {
            os << -antecedent_[i];
            if( i + 1 < antecedent_.size() ) os << " ";
        }
        if( !consequent_.empty() ) os << " ";
        for( size_t i = 0; i < consequent_.size(); ++i ) {
            os << consequent_[i];
            if( i + 1 < consequent_.size() ) os << " ";
        }
        os << " 0" << endl;
    }
    void print(ostream &os, const vector<const Literal*> &literals) const {
        for( size_t i = 0; i < antecedent_.size(); ++i ) {
            int index = (antecedent_[i] > 0 ? antecedent_[i] : (literals.size() >> 1) + -antecedent_[i]) - 1;
            literals[index]->print(os);
            if( i + 1 < antecedent_.size() ) os << " & ";
        }
        os << " => ";
        for( size_t i = 0; i < consequent_.size(); ++i ) {
            int index = (consequent_[i] > 0 ? consequent_[i] : (literals.size() >> 1) + -consequent_[i]) - 1;
            literals[index]->print(os);
            if( i + 1 < consequent_.size() ) os << " v ";
        }
        os << flush;
    }
};

template<typename T>
class Items : public set<T> {
    int num_;
  public:
    Items(int num = 0) : num_(num) { }
    ~Items() { }
    int num() const {
        return num_;
    }
    void dump(ostream &os) const {
        os << set<T>::size();
        for( typename set<T>::const_iterator it = set<T>::begin(); it != set<T>::end(); ++it )
            os << " " << *it;
        os << endl;
    }
    void read(istream &is) {
        for( size_t i = 0; i < num_; ++i ) {
            T item;
            is >> item;
            set<T>::insert(item);
        }
    }
    static const Items<T>* read_dump(istream &is) {
        int num;
        is >> num;
        Items<T> *items = new Items<T>(num);
        items->read(is);
        cout << "Items::read_dump: #items=" << items->num() << endl;
        return items;
    }
};

class Matrix {
  protected:
    const int num_states_;
    const int num_features_;
    const int last_numerical_feature_;
    const int first_boolean_feature_;
    int num_non_zero_entries_;
    vector<pair<int, pair<int, int>*> > rows_;
    vector<string> features_;
    map<string, int> feature_map_;

  public:
    Matrix(int num_states, int num_features, int last_numerical_feature, int first_boolean_feature)
      : num_states_(num_states),
        num_features_(num_features),
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
    int last_numerical_feature() const {
        return last_numerical_feature_;
    }
    int first_boolean_feature() const {
        return first_boolean_feature_;
    }
    int num_non_zero_entries() const {
        return num_non_zero_entries_;
    }
    const string& feature(int i) const {
        assert((0 <= i) && (i < num_features_));
        return features_[i];
    }
    int feature(const string &f) const {
        map<string, int>::const_iterator it = feature_map_.find(f);
        return it != feature_map_.end() ? it->second : -1;
    }
    bool numerical(int i) const {
        return i < last_numerical_feature_;
    }

    // accessors
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
    void get_row(int s, vector<int> &row) const {
        assert((0 <= s) && (s < rows_.size()));
        row = vector<int>(num_features_, 0);
        for( int i = 0; i < rows_[s].first; ++i )
            row[rows_[s].second[i].first] = rows_[s].second[i].second;
    }
    void get_column(int f, vector<int> &col) const {
        col = vector<int>(num_states_, 0);
        for( int s = 0; s < num_states_; ++s )
            col[f] = entry(s, f);;
    }
    int operator()(int s, int f) const {
        return entry(s, f);
    }

    // output
    void dump(ostream &os) const {
        os << num_states_ << " " << num_features_ << " " << last_numerical_feature_ << " " << first_boolean_feature_ << endl;
        for( size_t s = 0; s < rows_.size(); ++s ) {
            int count = rows_[s].first;;
            os << count;
            for( int i = 0; i < count; ++i )
                os << " " << rows_[s].second[i].first << ":" << rows_[s].second[i].second;
            os << endl;
        }
    }
    void print(ostream &os) const {
        os << "Matrix stats: #states=" << num_states_
           << ", #features=" << num_features_
           << ", last-numerical-feature=" << last_numerical_feature_
           << ", first-boolean-feature=" << first_boolean_feature_
           << ", #non-zero-entries=" << num_non_zero_entries_
           << endl;
        for( int s = 0; s < num_states_; ++s ) {
            os << "state " << s << ":";
            for( int f = 0; f < num_features_; ++f ) {
                int value = entry(s, f);
                if( value > 0 )
                    os << " " << f << ":" << value;
            }
            os << endl;
        }
    }

    // readers
    void read(istream &is) {
        // read features
        int n;
        is >> n;
        assert(n == num_features_);
        for( int i = 0; i < num_features_; ++i ) {
            string feature;
            is >> feature;
            feature_map_.insert(make_pair(feature, features_.size())); // there may be duplicate entries (e.g. 'dummy()')
            features_.push_back(feature);
        }

        // read rows
        rows_ = vector<pair<int, pair<int, int>*> >(num_states_, make_pair(0, nullptr));
        for( int i = 0; i < num_states_; ++i ) {
            int s, count, f, value;
            is >> s >> count;
            rows_[s] = make_pair(count, new pair<int, int>[count]);
            for( int j = 0; j < count; ++j ) {
                char filler;
                is >> f >> filler >> value;
                assert(filler == ':');
                assert((0 <= f) && (f < num_features_));
                assert(value > 0);
                assert((j == 0) || (rows_[s].second[j - 1].first < f)); // verify unique and ordered entries
                rows_[s].second[j] = make_pair(f, value);
                ++num_non_zero_entries_;
            }
        }
    }
    static const Matrix* read_dump(istream &is) {
        int num_states, num_features, last_numerical_feature, first_boolean_feature;
        is >> num_states >> num_features >> last_numerical_feature >> first_boolean_feature;
        Matrix *M = new Matrix(num_states, num_features, last_numerical_feature, first_boolean_feature);
        M->read(is);
        cout << "Matrix::read_dump: "
             << "#states=" << M->num_states()
             << ", #features=" << M->num_features()
             << ", #last-numerical-feature=" << M->last_numerical_feature()
             << ", #first-boolean-feature=" << M->first_boolean_feature()
             << ", #non-zero-entries=" << M->num_non_zero_entries()
             << endl;
        return M;
    }
};

class Transitions {
  protected:
    const int num_states_;
    const int num_transitions_;
    vector<pair<int, int*> > tr_;
    vector<int> offset_;

  public:
    Transitions(int num_states, int num_transitions)
      : num_states_(num_states), num_transitions_(num_transitions) {
        tr_ = vector<pair<int, int*> >(num_states_, make_pair(0, nullptr));
        offset_ = vector<int>(num_states_, 0);
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
    int successors(int s, vector<int> &succ) const {
        assert((0 <= s) && (s < num_states_));
        int count = tr_[s].first;
        if( count == 0 ) {
            succ.clear();
        } else {
            succ = vector<int>(count, 0);
            for( int i = 0; i < count; ++i )
                succ[i] = tr_[s].second[i];
        }
        return count;
    }

    // output
    void dump(ostream &os) const {
        os << num_states_ << " " << num_transitions_ << endl;

        int count = 0;
        for( size_t s = 0; s < num_states_; ++s )
            count += tr_[s].first > 0 ? 1 : 0;
        os << count << endl;

        for( size_t s = 0; s < num_states_; ++s ) {
            int count = tr_[s].first;
            if( count > 0 ) {
                os << s << " " << count;
                for( int i = 0; i < count; ++i )
                    os << " " << tr_[s].second[i];
                os << endl;
            }
        }
    }
    void print(ostream &os) const {
        os << "Transitions stats: #states=" << num_states_
           << ", #transitions=" << num_transitions_
           << endl;
        for( size_t s = 0; s < num_states_; ++s ) {
            int count = tr_[s].first;;
            if( count > 0 ) {
                os << "state " << s << ":";
                for( int i = 0; i < count; ++i )
                    os << " " << tr_[s].second[i];
                os << endl;
            }
        }
    }

    // readers
    void read(istream &is) {
        int num_records = 0, last_src = -1;
        is >> num_records;
        for( int i = 0; i < num_records; ++i ) {
            int src, count, dst;
            is >> src >> count;

            while( last_src < src ) {
                offset_[last_src + 1] = last_src == -1 ? 0 : offset_[last_src];
                ++last_src;
            }

            tr_[src] = make_pair(count, new int[count]);
            for( int j = 0; j < count; ++j ) {
                is >> dst;
                // std::cout << "src: " << src << ", j: " << j << " tr_[src].second[j - 1]: " << tr_[src].second[j - 1] << " dst: " << dst << " num_states: " << num_states_ << std::endl;
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
    static const Transitions* read_dump(istream &is) {
        int num_states, num_transitions;
        is >> num_states >> num_transitions;
        Transitions *T = new Transitions(num_states, num_transitions);
        T->read(is);
        cout << "Transitions::read_dump: #states=" << T->num_states() << ", #transitions=" << T->num_transitions() << endl;
        return T;
    }
};

class Theory {
  public:
    struct AbstractAction {
        string name_;
        int selected_precondition_;
        int precondition_;
        int selected_effect_;
        int effect_;

        static map<string, int> *feature_map_;

        AbstractAction() : selected_precondition_(0), precondition_(0), selected_effect_(0), effect_(0) { }
        AbstractAction(const string &name) : name_(name), selected_precondition_(0), precondition_(0), selected_effect_(0), effect_(0) { }
        void add_precondition(const string &feature, int prec) {
            if( (feature_map_ == nullptr) || (feature_map_->find(feature) == feature_map_->end()) ) {
                cout << "error: add_precondition: feature '" << feature << "' not found in feature map; did you forget --add-features to add features in forced actions?" << endl;
                exit(0);
            }
            assert((prec == 0) || (prec == 1));
            int k = feature_map_->find(feature)->second;
            selected_precondition_ = selected_precondition_ | (1 << k);
            precondition_ = precondition_ & ~(1 << k);
            precondition_ += prec << k;
        }
        void add_effect(const string &feature, int effect) {
            if( (feature_map_ == nullptr) || (feature_map_->find(feature) == feature_map_->end()) ) {
                cout << "error: add_precondition: feature '" << feature << "' not found in feature map; did you forget --add-features to add features in forced actions?" << endl;
                exit(0);
            }
            assert((effect == 0) || (effect == 1));
            int k = feature_map_->find(feature)->second;
            selected_effect_ = selected_effect_ | (1 << k);
            effect_ = effect_ & ~(1 << k);
            effect_ += effect << k;
        }
        bool operator<(const AbstractAction &a) const {
            return
              (selected_precondition_ < a.selected_precondition_) ||
              ((selected_precondition_ == a.selected_precondition_) && (precondition_ < a.precondition_)) ||
              ((selected_precondition_ == a.selected_precondition_) && (precondition_ == a.precondition_) && (selected_effect_ < a.selected_effect_)) ||
              ((selected_precondition_ == a.selected_precondition_) && (precondition_ == a.precondition_) && (selected_effect_ == a.selected_effect_) && (effect_ < a.effect_));
        }
        void print(ostream &os, const Theory &T, const vector<bool> &model, bool qnp_format) const {
            if( qnp_format ) {
                os << name_ << endl;

                int num_preconditions = 0;
                for( int k = 0; k < T.K_; ++k )
                    num_preconditions += selected_precondition_ & (1 << k) ? 1 : 0;
                os << num_preconditions;
                for( int k = 0; k < T.K_; ++k ) {
                    if( selected_precondition_ & (1 << k) ) {
                        os << " " << T.condition_on_selected_feature(model, k, (precondition_ & (1 << k)) != 0, qnp_format);
                    }
                }
                os << endl;

                int num_effects = 0;
                for( int k = 0; k < T.K_; ++k )
                    num_effects += selected_effect_ & (1 << k) ? 1 : 0;
                os << num_effects;
                for( int k = 0; k < T.K_; ++k ) {
                    if( selected_effect_ & (1 << k) ) {
                        os << " " << T.effect_on_selected_feature(model, k, (effect_ & (1 << k)) != 0, qnp_format);
                    }
                }
                os << endl;
            } else {
                os << "Prec={";
                for( int k = 0; k < T.K_; ++k ) {
                    if( selected_precondition_ & (1 << k) )
                        os << " " << T.condition_on_selected_feature(model, k, (precondition_ & (1 << k)) != 0, qnp_format) << ",";
                }
                os << " }; Effect={";
                for( int k = 0; k < T.K_; ++k ) {
                    if( selected_effect_ & (1 << k) )
                        os << " " << T.effect_on_selected_feature(model, k, (effect_ & (1 << k)) != 0) << ",";
                }
                os << " }" << endl;
            }
        }
        void dump(ostream &os) const {
            os << bit_string(selected_precondition_)
               << " " << bit_string(precondition_)
               << " " << bit_string(selected_effect_)
               << " " << bit_string(effect_)
               << endl;
        }
    };

  public:
    struct Options {
        string prefix_;
        bool ordering_;
        bool decode_;
        bool explicit_goals_;
        bool add_actions_;
        bool add_features_;
        bool only_soundness_;
        Options()
          : ordering_(true),
            decode_(false),
            explicit_goals_(false),
            add_actions_(false),
            add_features_(false),
            only_soundness_(false) {
        }
    };

  protected:
    const Matrix &M_;
    const Transitions &T_;
    const int K_;           // number of features in abstraction
    const int N_;           // number of actions in abstraction

    const Options options_;
    Items<AbstractAction> actions_;
    Items<string> features_;
    Items<int> goals_;

    const int num_states_;
    const int num_features_;
    int log_num_features_;
    int log_last_numerical_feature_;

    vector<pair<int, string> > var_offsets_;
    vector<const Var*> variables_;
    vector<const Literal*> literals_;

    vector<pair<int, const string> > comments_;
    vector<const Implication*> implications_;
    vector<pair<int, string> > imp_offsets_;

    mutable bool satisfiable_;
    mutable vector<bool> model_;

  public:
    Theory(const Matrix &M, const Transitions &T, int K, int N, const Options &options)
      : M_(M), T_(T), K_(K), N_(N), options_(options),
        num_states_(M_.num_states()),
        num_features_(M_.num_features()),
        log_num_features_(0),
        log_last_numerical_feature_(0),
        satisfiable_(false) {

        while( (1 << log_num_features_) < num_features_ )
            ++log_num_features_;
        while( (1 << log_last_numerical_feature_) < M_.last_numerical_feature() )
            ++log_last_numerical_feature_;

        cout << "Theory: parameters:"
             << " K=" << K_
             << ", N=" << N_
             << ", log2(#features)=" << log_num_features_
             << ", log2(last_numerical_feature)=" << log_last_numerical_feature_
             << endl;

        if( options_.add_features_ )
            read_features(options_.prefix_);
        if( options_.explicit_goals_ )
            read_goals(options_.prefix_);
        if( options_.add_actions_ )
            read_actions(options_.prefix_);

        cout << "-------- variables + literals --------" << endl;
        build_variables();
        build_literals();
        cout << "--------- base implications ----------" << endl;
        build_base();

        add_implications_for_features();
        add_implications_for_goals();
        add_implications_for_actions();

        cout << "------ soundness + completeness ------" << endl;
        enforce_soundness();
        if( !options_.only_soundness_ ) enforce_completeness();
        cout << "--------------------------------------" << endl;
    }
    virtual ~Theory() {
        for( size_t i = 0; i < implications_.size(); ++i )
            delete implications_[i];
        for( size_t i = 0; i < literals_.size(); ++i )
            delete literals_[i];
        for( size_t i = 0; i < variables_.size(); ++i )
            delete variables_[i];
    }

    const Var& variable(int index) const {
        assert((0 <= index) && (index < variables_.size()));
        return *variables_[index];
    }
    const Literal& literal(int index) const {
        assert(index != 0);
        assert((-int(literals_.size()) <= index) && (index <= int(literals_.size())));
        return index > 0 ? *literals_[index - 1] : *literals_[variables_.size() + -index - 1];
    }
    int num_variables() const {
        return variables_.size();
    }

    void add_implication(const Implication *IP) {
        implications_.push_back(IP);
    }
    void add_comment(const string &comment) {
        comments_.push_back(make_pair(implications_.size(), comment));
    }
    int num_implications() const {
        return implications_.size();
    }

    bool satisfiable() const {
        return satisfiable_;
    }
    const vector<bool>& model() const {
        return model_;
    }

    int index_for_selected_feature(const vector<bool> &model, int j) const {
        int index = 0;
        for( int k = 0; k < log_num_features_; ++k ) {
            assert(F(j, k) < model.size());
            index += model[F(j, k)] ? (1 << k) : 0;
        }
        return index;
    }
    string condition_on_selected_feature(const vector<bool> &model, int k, bool prec, bool qnp_format = false) const {
        // get feature index from model
        int f = index_for_selected_feature(model, k);

        // generate string
        string feature;
        if( qnp_format ) {
            feature += M_.feature(f);
            feature += prec ? " 1" : " 0";
        } else if( M_.numerical(f) ) {
            feature += to_string(f) + "." + M_.feature(f);
            feature += prec ? ">0" : "=0";
        } else {
            feature += prec ? "" : "-";
            feature += to_string(f) + "." + M_.feature(f);
        }
        return feature;
    }
    string effect_on_selected_feature(const vector<bool> &model, int k, bool eff, bool qnp_format = false) const {
        // get feature index from model
        int f = index_for_selected_feature(model, k);

        // generate string
        string feature;
        if( qnp_format ) {
            feature += M_.feature(f);
            feature += eff ? " 1" : " 0";
        } else if( M_.numerical(f) ) {
            feature += eff ? "INC(" : "DEC(";
            feature += to_string(f) + "." + M_.feature(f) + ")";
        } else {
            feature += eff ? "" : "-";
            feature += to_string(f) + "." + M_.feature(f);
        }
        return feature;
    }

    void build_variables() {
        // A variables: 4 x K x N
        var_offsets_.push_back(make_pair(0, "A"));
        for( int t = 1; t < 5; ++t ) {
            for( int j = 0; j < N_; ++j ) {
                for( int k = 0; k < K_; ++k ) {
                    int index = variables_.size();
                    string name = string("A(") + to_string(t) + "," + to_string(j) + "," + to_string(k) + ")";
                    variables_.push_back(new Var(index, name));
                    assert(index == A(t, j, k));
                }
            }
        }
        cout << "A: #variables=" << variables_.size() - var_offsets_.back().first << " (4KN), offset=" << var_offsets_.back().first << endl;

        // F variables: K x log2(F)
        var_offsets_.push_back(make_pair(variables_.size(), "F"));
        for( int j = 0; j < K_; ++j ) {
            for( int k = 0; k < log_num_features_; ++k ) {
                int index = variables_.size();
                string name = string("F(") + to_string(j) + "," + to_string(k) + ")";
                variables_.push_back(new Var(index, name));
                assert(index == F(j, k));
            }
        }
        cout << "F: #variables=" << variables_.size() - var_offsets_.back().first << " (KL) where L=log2(F), offset=" << var_offsets_.back().first << endl;

        // phi and phi* variables: 2 x S x K
        var_offsets_.push_back(make_pair(variables_.size(), "phi"));
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                int index = variables_.size();
                string name = string("phi(") + to_string(k) + "," + to_string(i) + ")";
                variables_.push_back(new Var(index, name));
                assert(index == phi(k, i, false));
            }
        }
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                int index = variables_.size();
                string name = string("phi*(") + to_string(k) + "," + to_string(i) + ")";
                variables_.push_back(new Var(index, name));
                assert(index == phi(k, i, true));
            }
        }
        cout << "phi: #variables=" << variables_.size() - var_offsets_.back().first << " (2SK), offset=" << var_offsets_.back().first << endl;

        // B variables: 2 x S x K x N
        var_offsets_.push_back(make_pair(variables_.size(), "B"));
        for( int t = 1; t < 3; ++t ) {
            for( int j = 0; j < N_; ++j ) {
                for( int i = 0; i < num_states_; ++i ) {
                    for( int k = 0; k < K_; ++k ) {
                        int index = variables_.size();
                        string name = string("B(") + to_string(t) + "," + to_string(j) + "," + to_string(i) + "," + to_string(k) + ")";
                        variables_.push_back(new Var(index, name));
                        assert(index == B(t, j, i, k));
                    }
                }
            }
        }
        cout << "B: #variables=" << variables_.size() - var_offsets_.back().first << " (2SKN), offset=" << var_offsets_.back().first << endl;

        // RES variables: T x N
        var_offsets_.push_back(make_pair(variables_.size(), "RES"));
        for( int j = 0; j < N_; ++j ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < T_.num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("RES(") + to_string(i) + "," + to_string(l) + "," + to_string(j) + ")";
                    variables_.push_back(new Var(index, name));
                    assert(index == RES(i, l, j));
                }
            }
        }
        cout << "RES: #variables=" << variables_.size() - var_offsets_.back().first << " (TN), offset=" << var_offsets_.back().first << endl;

        // INC variables: T x K
        var_offsets_.push_back(make_pair(variables_.size(), "INC"));
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < T_.num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("INC(") + to_string(k) + "," + to_string(i) + "," + to_string(l) + ")";
                    variables_.push_back(new Var(index, name));
                    assert(index == INC(k, i, l));
                }
            }
        }
        cout << "INC: #variables=" << variables_.size() - var_offsets_.back().first << " (TK), offset=" << var_offsets_.back().first << endl;

        // DEC variables: T x K
        var_offsets_.push_back(make_pair(variables_.size(), "DEC"));
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < T_.num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("DEC(") + to_string(k) + "," + to_string(i) + "," + to_string(l) + ")";
                    variables_.push_back(new Var(index, name));
                    assert(index == DEC(k, i, l));
                }
            }
        }
        cout << "DEC: #variables=" << variables_.size() - var_offsets_.back().first << " (TK), offset=" << var_offsets_.back().first << endl;

        // MATCH variables: T x N
        var_offsets_.push_back(make_pair(variables_.size(), "MATCH"));
        for( int j = 0; j < N_; ++j ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < T_.num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("MATCH(") + to_string(i) + "," + to_string(l) + "," + to_string(j) + ")";
                    variables_.push_back(new Var(index, name));
                    assert(index == MATCH(i, l, j));
                }
            }
        }
        cout << "MATCH: #variables=" << variables_.size() - var_offsets_.back().first << " (TN), offset=" << var_offsets_.back().first << endl;

        // APPMATCH variables: T x N
        var_offsets_.push_back(make_pair(variables_.size(), "APPMATCH"));
        for( int j = 0; j < N_; ++j ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < T_.num_transitions(i); ++l ) {
                    int index = variables_.size();
                    string name = string("APPMATCH(") + to_string(i) + "," + to_string(l) + "," + to_string(j) + ")";
                    variables_.push_back(new Var(index, name));
                    assert(index == APPMATCH(i, l, j));
                }
            }
        }
        cout << "APPMATCH: #variables=" << variables_.size() - var_offsets_.back().first << " (TN), offset=" << var_offsets_.back().first << endl;

        bool disable_ordering = (features_.size() == K_) || ((features_.size() < K_) && !actions_.empty());
        if( options_.ordering_ && !disable_ordering ) {
            // df variables:
            var_offsets_.push_back(make_pair(variables_.size(), "df"));
            for( int j = 0; j + 1 < K_; ++j ) {
                for( int k = 0; k < log_num_features_; ++k ) {
                    int index = variables_.size();
                    string name = string("df(") + to_string(j) + "," + to_string(k) + ")";
                    variables_.push_back(new Var(index, name));
                    assert(index == df(j, k));
                }
            }
            cout << "df: #variables=" << variables_.size() - var_offsets_.back().first << " ((K-1) x log2(F)), offset=" << var_offsets_.back().first << endl;
        }

        add_variables_for_features();
        add_variables_for_goals();
    }
    void build_literals() {
        for( size_t i = 0; i < variables_.size(); ++i )
            literals_.push_back(new Literal(*variables_[i], false));
        for( size_t i = 0; i < variables_.size(); ++i )
            literals_.push_back(new Literal(*variables_[i], true));
    }

    // references to variables
    // A(1, j, k) iff j-th action *affects* k-th feature
    // A(2, j, k) iff type of effect of j-th action on k-th feature (0=false/decrement, 1=true/increment)
    // A(3, j, k) iff k-th feature is *precondition* of j-th action
    // A(4, j, k) iff type of precondition of k-th feature in j-th action (0=false/=0, 1=true/>0)
    int A(int t, int j, int k) const { // t in {1, ... ,4}, j in {0, ..., N-1}, k in {0, ... ,K-1}
        assert((0 < t) && (t < 5));
        assert((0 <= j) && (j < N_));
        assert((0 <= k) && (k < K_));
        int index = k + K_ * (j + N_ * (t - 1));
        assert(var_offsets_[0].second == "A");
        return var_offsets_[0].first + index;
    }

    // F(j, k) iff j-th feature encodes feature whose index has k-th set to 1 (logarithmic encoding)
    int F(int j, int k) const { // j in {0, ..., K-1}, k in {0, ..., log2(F)-1}
        assert((0 <= j) && (j < K_));
        assert((0 <= k) && (k < log_num_features_));
        int index = k + log_num_features_ * j;
        assert(var_offsets_[1].second == "F");
        return var_offsets_[1].first + index;
    }
    int df(int j, int k) const { // j in {0, ..., K-2}, k in {0, ..., log2(F)-1}
        assert((0 <= j) && (j + 1 < K_));
        assert((0 <= k) && (k < log_num_features_));
        int index = k + log_num_features_ * j;
        assert(var_offsets_[9].second == "df");
        return var_offsets_[9].first + index;
    }
    int MappedBy(int j, int k) const { // j in {0, ..., #forced-features-1}, k in {0, ..., K-1}
        assert((0 <= j) && (j < int(features_.size())));
        assert((0 <= k) && (k < K_));
        int index = k + K_ * j;
        assert(var_offsets_[10].second == "MappedBy");
        return var_offsets_[10].first + index;
    }

    // phi(k, i, F) iff k-th feature *holds* in state s_i in sample
    // phi(k, i, T) iff k-th feature *does not hold* in state s_i in sample
    int phi(int k, int i, bool star) const { // k in {0, ..., K-1}, i in {0, ..., states-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        int index = i + num_states_ * (k + K_ * int(star));
        assert(var_offsets_[2].second == "phi");
        return var_offsets_[2].first + index;
    }

    // used for soundness
    int B(int t, int j, int i, int k) const { // t in {1,...,2}, j in {0, ..., N-1}, i in {0, ..., states-1}, k in {0, ..., K-1}
        assert((1 <= t) && (t < 3));
        assert((0 <= j) && (j < N_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= k) && (k < K_));
        int index = k + K_ * (i + num_states_ * (j + N_ * (t - 1)));
        assert(var_offsets_[3].second == "B");
        return var_offsets_[3].first + index;
    }

    // used for soundness/completeness
    int RES(int i, int l, int j) const { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int index = (T_.offset(i) + l) + T_.num_transitions() * j;
        assert(var_offsets_[4].second == "RES");
        return var_offsets_[4].first + index;
    }

    // used for soundness/completeness
    int INC(int k, int i, int l) const { // k in {0, ..., K}, (i,l) in {0, ..., transitions-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        int index = (T_.offset(i) + l) + T_.num_transitions() * k;
        assert(var_offsets_[5].second == "INC");
        return var_offsets_[5].first + index;
    }

    // used for soundness/completeness
    int DEC(int k, int i, int l) const { // k in {0, ..., K}, (i,l) in {0, ..., transitions-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        int index = (T_.offset(i) + l) + T_.num_transitions() * k;
        assert(var_offsets_[6].second == "DEC");
        return var_offsets_[6].first + index;
    }

    // used for soundness/completeness
    int MATCH(int i, int l, int j) const { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int index = (T_.offset(i) + l) + T_.num_transitions() * j;
        assert(var_offsets_[7].second == "MATCH");
        return var_offsets_[7].first + index;
    }

    // used only for completeness
    int APPMATCH(int i, int l, int j) const { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int index = (T_.offset(i) + l) + T_.num_transitions() * j;
        assert(var_offsets_[8].second == "APPMATCH");
        return var_offsets_[8].first + index;
    }

    // construction of base formulas
    void build_formulas_A() {
        // constraints for A formulas
        for( int j = 0; j < N_; ++j ) {
            for( int k = 0; k < K_; ++k ) {
                // if effect on F^k, it must be selected
                Implication *IP1 = new Implication;
                IP1->add_antecedent(1 + A(2, j, k));
                IP1->add_consequent(1 + A(1, j, k));
                add_implication(IP1);

                // if prec on F^k, it must be selected
                Implication *IP2 = new Implication;
                IP2->add_antecedent(1 + A(4, j, k));
                IP2->add_consequent(1 + A(3, j, k));
                add_implication(IP2);

                // if F^k is numeric and it is decremented, precondition must be value is > 0
                // A(1,j,k) & -A(2,j,k) & [F^k is numeric] => A(4,j,k)
                // [F^k is numeric] = -F(k,t) for all t >= index_for_numerical
                Implication *IP3 = new Implication;
                IP3->add_antecedent(1 + A(1, j, k));
                IP3->add_antecedent(-(1 + A(2, j, k)));
                for( int t = log_last_numerical_feature_; t < log_num_features_; ++t )
                    IP3->add_antecedent(-(1 + F(k, t)));
                IP3->add_consequent(1 + A(4, j, k));
                add_implication(IP3);
            }
        }

        // some feature must be affected (selected)
        for( int j = 0; j < N_; ++j ) {
            Implication *IP = new Implication;
            for( int k = 0; k < K_; ++k )
                IP->add_consequent(1 + A(1, j, k));
            add_implication(IP);
        }

#if 0 // CHECK
        // selected action must be applicable at some state
        for( int j = 0; j < N_; ++j ) {
            Implication *IP = new Implication;
            for( int i = 0; i < num_states_; ++i )
                IP->add_consequent(1 + APP(i, j));
            add_implication(IP);
        }
#endif

        // ordering: if first affected feature by A^j is F^k, then A^{j-1} must affect some feature in {F^1, ..., F^k}
        if( options_.ordering_ && !options_.add_actions_ ) {
            for( int j = 1; j < N_; ++j ) {
                for( int k = 0; k < K_; ++k ) {
                    Implication *IP = new Implication;
                    for( int t = 0; t < k; ++t )
                        IP->add_antecedent(-(1 + A(1, j, t)));
                    IP->add_antecedent(1 + A(1, j, k));
                    for( int t = 0; t <= k; ++t )
                        IP->add_consequent(1 + A(1, j - 1, t));
                    add_implication(IP);
                }
            }
        }
    }

    void build_formulas_F() {
        // constraints: dummy features are not selected
        for( int j = 0; j < K_; ++j ) {
            for( int index = M_.last_numerical_feature(); index < M_.first_boolean_feature(); ++index ) {
                Implication *IP = new Implication;
                for( int k = 0; k < log_num_features_; ++k )
                    IP->add_antecedent(bit(k, index) ? 1 + F(j, k) : -(1 + F(j, k)));
                add_implication(IP);
            }
        }

        // constraints: inexistent features are not selected
        for( int j = 0; j < K_; ++j ) {
            for( int index = num_features_; index < (1 << log_num_features_); ++index ) {
                Implication *IP = new Implication;
                for( int k = 0; k < log_num_features_; ++k )
                    IP->add_antecedent(bit(k, index) ? 1 + F(j, k) : -(1 + F(j, k)));
                add_implication(IP);
            }
        }

        // constraints: ordering of selected features F^0 < F^0 < ... < F^{K-1}
        bool disable_ordering = (features_.size() == K_) || ((features_.size() < K_) && !actions_.empty());
        if( options_.ordering_ && !disable_ordering ) {
            for( int k = 0; k + 1 < K_; ++k ) {
#if 1 // CHECK: for some reason, it is more efficient to generate these implications here
                for( int t = 0; t < log_num_features_; ++t ) {
                    Implication *IP1 = new Implication;
                    IP1->add_antecedent(1 + df(k, t));
                    IP1->add_consequent(-(1 + F(k, t)));
                    add_implication(IP1);

                    Implication *IP2 = new Implication;
                    IP2->add_antecedent(1 + df(k, t));
                    IP2->add_consequent(1 + F(k + 1, t));
                    add_implication(IP2);

                    for( int i = t + 1; i < log_num_features_; ++i ) {
                        Implication *IQ1 = new Implication;
                        IQ1->add_antecedent(1 + df(k, t));
                        IQ1->add_antecedent(1 + F(k, i));
                        IQ1->add_consequent(1 + F(k + 1, i));
                        add_implication(IQ1);

                        Implication *IQ2 = new Implication;
                        IQ2->add_antecedent(1 + df(k, t));
                        IQ2->add_antecedent(-(1 + F(k, i)));
                        IQ2->add_consequent(-(1 + F(k + 1, i)));
                        add_implication(IQ2);
                    }
                }
#endif

                Implication *IP = new Implication;
                for( int t = 0; t < log_num_features_; ++t )
                    IP->add_consequent(1 + df(k, t));
                add_implication(IP);
            }
        }
    }

    void build_phi(int k, int i, bool star) {
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));

        // forward implications
        for( int l = 0; l < num_features_; ++l ) {
            if( (!star && !M_(i, l)) || (star && M_(i, l)) ) {
                Implication *IP = new Implication;
                IP->add_antecedent(1 + phi(k, i, star));
                for( int t = 0; t < log_num_features_; ++t )
                    IP->add_consequent(bit(t, l) ? -(1 + F(k, t)) : (1 + F(k, t)));
                add_implication(IP);
            }
        }

        // backward implications
        for( int l = 0; l < num_features_; ++l ) {
            if( (!star && M_(i, l)) || (star && !M_(i, l)) ) {
                Implication *IP = new Implication;
                IP->add_consequent(1 + phi(k, i, star));
                for( int t = 0; t < log_num_features_; ++t )
                    IP->add_antecedent(bit(t, l) ? (1 + F(k, t)) : -(1 + F(k, t)));
                add_implication(IP);
            }
        }
    }
    void build_formulas_phi() {
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                build_phi(k, i, false);
                build_phi(k, i, true);
            }
        }
    }

    void build_B(int t, int j, int i, int k) { // t in {1,...,2}, j in {0, ..., N-1}, i in {0, ..., states-1}, k in {0, ..., K-1}
        assert((1 <= t) && (t < 3));
        assert((0 <= j) && (j < N_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= k) && (k < K_));
        // B(1,j,i,k) = A(4,j,k) -> phi(k,i)
        // B(2,j,i,k) = A(3,j,k) & -A(4,j,k) -> phi*(k,i)

        bool star = t == 1;

        // forward implication: one
        Implication *IP1 = new Implication;
        IP1->add_antecedent(1 + B(t, j, i, k));
        if( star ) IP1->add_antecedent(1 + A(3, j, k));
        IP1->add_antecedent(star ? -(1 + A(4, j, k)) : 1 + A(4, j, k));
        IP1->add_consequent(1 + phi(k, i, star));
        add_implication(IP1);

        // backward implications: two or three
        Implication *IP2 = new Implication;
        IP2->add_antecedent(star ? 1 + A(4, j, k) : -(1 + A(4, j, k)));
        IP2->add_consequent(1 + B(t, j, i, k));
        add_implication(IP2);

        Implication *IP3 = new Implication;
        IP3->add_antecedent(1 + phi(k, i, star));
        IP3->add_consequent(1 + B(t, j, i, k));
        add_implication(IP3);

        if( star ) {
            Implication *IP4 = new Implication;
            IP4->add_antecedent(-(1 + A(3, j, k)));
            IP4->add_consequent(1 + B(t, j, i, k));
            add_implication(IP4);
        }
    }
    void build_formulas_B() {
        for( int t = 1; t < 3; ++t ) {
            for( int j = 0; j < N_; ++j ) {
                for( int i = 0; i < num_states_; ++i ) {
                    for( int k = 0; k < K_; ++k )
                        build_B(t, j, i, k);
                }
            }
        }
    }

    void build_RES(int i, int l, int j) { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int state_l = T_.transition_state(i, l);

        // forward implications: 4 of them
        for( int k = 0; k < K_; ++k ) {
            // RES(i,l,j) & A(2,j,k) -> phi(k,l)
            Implication *IP1 = new Implication;
            IP1->add_antecedent(1 + RES(i, l, j));
            IP1->add_antecedent(1 + A(2, j, k));
            IP1->add_consequent(1 + phi(k, state_l, false));
            add_implication(IP1);

            // RES(i,l,j) & A(1,j,k) & -A(2,j,k) & -[F^k is numeric] -> phi*(k,l)
            for( int t = log_last_numerical_feature_; t < log_num_features_; ++t ) {
                Implication *IP2 = new Implication;
                IP2->add_antecedent(1 + RES(i, l, j));
                IP2->add_antecedent(1 + A(1, j, k));
                IP2->add_antecedent(-(1 + A(2, j, k)));
                IP2->add_antecedent(1 + F(k, t)); // -[F^k is numeric]
                IP2->add_consequent(1 + phi(k, state_l, true));
                add_implication(IP2);
            }

            // RES(i,l,j) & -A(1,j,k) & phi(k,i) -> phi(k,l)
            Implication *IP3 = new Implication;
            IP3->add_antecedent(1 + RES(i, l, j));
            IP3->add_antecedent(-(1 + A(1, j, k)));
            IP3->add_antecedent(1 + phi(k, i, false));
            IP3->add_consequent(1 + phi(k, state_l, false));
            add_implication(IP3);

            // RES(i,l,j) & -A(1,j,k) & phi*(k,i) -> phi*(k,l)
            Implication *IP4 = new Implication;
            IP4->add_antecedent(1 + RES(i, l, j));
            IP4->add_antecedent(-(1 + A(1, j, k)));
            IP4->add_antecedent(1 + phi(k, i, true));
            IP4->add_consequent(1 + phi(k, state_l, true));
            add_implication(IP4);
        }
    }
    void build_formulas_RES() {
        for( int j = 0; j < N_; ++j ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < T_.num_transitions(i); ++l )
                    build_RES(i, l, j);
            }
        }
    }

    void build_INC(int k, int i, int l) { // k in {0, ..., K-1}, (i,l) in {0, ..., transitions-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        int state_l = T_.transition_state(i, l);

        // forward implications
        for( int f = 0; f < num_features_; ++f ) {
            if( M_(i, f) >= M_(state_l, f) ) {
                Implication *IP = new Implication;
                IP->add_antecedent(1 + INC(k, i, l));
                for( int z = 0; z < log_num_features_; ++z )
                    IP->add_consequent(bit(z, f) ? -(1 + F(k, z)) : (1 + F(k, z)));
                add_implication(IP);
            }
        }

        // backward implications
        for( int f = 0; f < num_features_; ++f ) {
            if( M_(i, f) < M_(state_l, f) ) {
                Implication *IP = new Implication;
                for( int z = 0; z < log_num_features_; ++z )
                    IP->add_antecedent(bit(z, f) ? 1 + F(k, z) : -(1 + F(k, z)));
                IP->add_consequent(1 + INC(k, i, l));
                add_implication(IP);
            }
        }
    }
    void build_formulas_INC() {
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < T_.num_transitions(i); ++l )
                    build_INC(k, i, l);
            }
        }
    }

    void build_DEC(int k, int i, int l) { // k in {0, ..., K-1}, (i,l) in {0, ..., transitions-1}
        assert((0 <= k) && (k < K_));
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        int state_l = T_.transition_state(i, l);

        // forward implications
        for( int f = 0; f < num_features_; ++f ) {
            if( M_(i, f) <= M_(state_l, f) ) {
                Implication *IP = new Implication;
                IP->add_antecedent(1 + DEC(k, i, l));
                for( int z = 0; z < log_num_features_; ++z )
                    IP->add_consequent(bit(z, f) ? -(1 + F(k, z)) : (1 + F(k, z)));
                add_implication(IP);
            }
        }

        // backward implications
        for( int f = 0; f < num_features_; ++f ) {
            if( M_(i, f) > M_(state_l, f) ) {
                Implication *IP = new Implication;
                for( int z = 0; z < log_num_features_; ++z )
                    IP->add_antecedent(bit(z, f) ? 1 + F(k, z) : -(1 + F(k, z)));
                IP->add_consequent(1 + DEC(k, i, l));
                add_implication(IP);
            }
        }
    }
    void build_formulas_DEC() {
        for( int k = 0; k < K_; ++k ) {
            for( int i = 0; i < num_states_; ++i ) {
                for( int l = 0; l < T_.num_transitions(i); ++l )
                    build_DEC(k, i, l);
            }
        }
    }

    void build_MATCH(int i, int l, int j) { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int state_l = T_.transition_state(i, l);

        // forward: implication for RES
        Implication *IP = new Implication;
        IP->add_antecedent(1 + MATCH(i, l, j));
        IP->add_consequent(1 + RES(i, l, j));
        add_implication(IP);

        // forward: increments for INCs
        for( int k = 0; k < K_; ++k ) {
            Implication *IP1 = new Implication;
            IP1->add_antecedent(1 + MATCH(i, l, j));
            // F^k is numeric
            for( int z = log_last_numerical_feature_; z < log_num_features_; ++z )
                IP1->add_antecedent(-(1 + F(k, z)));
            IP1->add_antecedent(1 + A(2, j, k));
            IP1->add_consequent(1 + INC(k, i, l));
            add_implication(IP1);

            Implication *IP2 = new Implication;
            IP2->add_antecedent(1 + MATCH(i, l, j));
            // F^k is numeric
            for( int z = log_last_numerical_feature_; z < log_num_features_; ++z )
                IP2->add_antecedent(-(1 + F(k, z)));
            IP2->add_antecedent(1 + INC(k, i, l));
            IP2->add_consequent(1 + A(2, j, k));
            add_implication(IP2);
        }

        // forward: increments for DECs
        for( int k = 0; k < K_; ++k ) {
            Implication *IP1 = new Implication;
            IP1->add_antecedent(1 + MATCH(i, l, j));
            // F^k is numeric
            for( int z = log_last_numerical_feature_; z < log_num_features_; ++z )
                IP1->add_antecedent(-(1 + F(k, z)));
            IP1->add_antecedent(1 + A(1, j, k));
            IP1->add_antecedent(-(1 + A(2, j, k)));
            IP1->add_consequent(1 + DEC(k, i, l));
            add_implication(IP1);

            Implication *IP2 = new Implication;
            IP2->add_antecedent(1 + MATCH(i, l, j));
            // F^k is numeric
            for( int z = log_last_numerical_feature_; z < log_num_features_; ++z )
                IP2->add_antecedent(-(1 + F(k, z)));
            IP2->add_antecedent(1 + DEC(k, i, l));
            IP2->add_consequent(1 + A(1, j, k));
            add_implication(IP2);

            Implication *IP3 = new Implication;
            IP3->add_antecedent(1 + MATCH(i, l, j));
            // F^k is numeric
            for( int z = log_last_numerical_feature_; z < log_num_features_; ++z )
                IP3->add_antecedent(-(1 + F(k, z)));
            IP3->add_antecedent(1 + DEC(k, i, l));
            IP3->add_consequent(-(1 + A(2, j, k)));
            add_implication(IP3);
        }
    }
    void build_formulas_MATCH() {
        for( int i = 0; i < num_states_; ++i ) {
            for( int l = 0; l < T_.num_transitions(i); ++l ) {
                for( int j = 0; j < N_; ++j )
                    build_MATCH(i, l, j);
            }
        }
    }

    void build_APPMATCH(int i, int l, int j) { // (i,l) in {0, ..., transitions-1}, j in {0, ..., N-1}
        assert((0 <= i) && (i < num_states_));
        assert((0 <= l) && (l < T_.num_transitions(i)));
        assert((0 <= j) && (j < N_));
        int state_l = T_.transition_state(i, l);

        // forward implications
        Implication *IP = new Implication;
        IP->add_antecedent(1 + APPMATCH(i, l, j));
        IP->add_consequent(1 + MATCH(i, l, j));
        add_implication(IP);

        for( int k = 0; k < K_; ++k ) {
            for( int t = 1; t < 3; ++t ) {
                Implication *IP = new Implication;
                IP->add_antecedent(1 + APPMATCH(i, l, j));
                IP->add_consequent(1 + B(t, j, i, k));
                add_implication(IP);
            }
        }
    }
    void build_formulas_APPMATCH() {
        for( int i = 0; i < num_states_; ++i ) {
            for( int l = 0; l < T_.num_transitions(i); ++l ) {
                for( int j = 0; j < N_; ++j )
                    build_APPMATCH(i, l, j);
            }
        }
    }

    void build_df(int j, int k) { // j in {0, ..., K-2}, k in {0, ..., log(F)}
        assert((0 <= j) && (j + 1 < K_));
        assert((0 <= k) && (k < log_num_features_));

#if 0 // CHECK: generated with formulas for F as it is more efficient (for unknown reasons... see above)
        Implication *IP1 = new Implication;
        IP1->add_antecedent(1 + df(j, k));
        IP1->add_consequent(-(1 + F(j, k)));
        add_implication(IP1);

        Implication *IP2 = new Implication;
        IP2->add_antecedent(1 + df(j, k));
        IP2->add_consequent(1 + F(j + 1, k));
        add_implication(IP2);

        for( int i = k + 1; i < log_num_features_; ++i ) {
            Implication *IP1 = new Implication;
            IP1->add_antecedent(1 + df(j, k));
            IP1->add_antecedent(1 + F(j, i));
            IP1->add_consequent(1 + F(j + 1, i));
            add_implication(IP1);

            Implication *IP2 = new Implication;
            IP2->add_antecedent(1 + df(j, k));
            IP2->add_antecedent(1 + F(j + 1, i));
            IP2->add_consequent(1 + F(j, i));
            //IP2->add_antecedent(-(1 + F(j, i)));
            //IP2->add_consequent(-(1 + F(j + 1, i)));
            add_implication(IP2);
        }
#endif
    }
    void build_formulas_df() {
        for( int j = 0; j + 1 < K_; ++j ) {
            for( int k = 0; k < log_num_features_; ++k )
                build_df(j, k);
        }
    }

    void build_base() {
        imp_offsets_.push_back(make_pair(0, "A"));
        add_comment("Theory for formulas A(t,j,k)");
        build_formulas_A();
        cout << "A: #implications=" << implications_.size() - imp_offsets_.back().first << " (3NK + N + (N-1) x K = N(1 + 4K) - K)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "F"));
        add_comment("Theory for formulas F(j,k)");
        build_formulas_F();
        cout << "F: #implications=" << implications_.size() - imp_offsets_.back().first << " (K x (D+U) + (K-1)) where D=#dummy-features and U=#unused-features" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "phi"));
        add_comment("Theory for formulas phi(k,i)");
        build_formulas_phi();
        cout << "phi: #implications=" << implications_.size() - imp_offsets_.back().first << " (2SK x F)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "B"));
        add_comment("Theory for formulas B(t,j,i,k)");
        build_formulas_B();
        cout << "B: #implications=" << implications_.size() - imp_offsets_.back().first << " (7SKN)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "RES"));
        add_comment("Theory for formulas RES(i,l,j)");
        build_formulas_RES();
        cout << "RES: #implications=" << implications_.size() - imp_offsets_.back().first << " (3TNK + TNK x Delta(L)) where Delta(L)=log2(#features)-log2(last_numerical_feature)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "INC"));
        add_comment("Theory for formulas INC(k,i,l)");
        build_formulas_INC();
        cout << "INC: #implications=" << implications_.size() - imp_offsets_.back().first << " (TKF)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "DEC"));
        add_comment("Theory for formulas DEC(k,i,l)");
        build_formulas_DEC();
        cout << "DEC: #implications=" << implications_.size() - imp_offsets_.back().first << " (TKF)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "MATCH"));
        add_comment("Theory for formulas MATCH(i,l,j)");
        build_formulas_MATCH();
        cout << "MATCH: #implications=" << implications_.size() - imp_offsets_.back().first << " (TN + 5TNK)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "APPMATCH"));
        add_comment("Theory for formulas APPMATCH(i,l,j)");
        build_formulas_APPMATCH();
        cout << "APPMATCH: #implications=" << implications_.size() - imp_offsets_.back().first << " (TN + 2TNK)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());

        imp_offsets_.push_back(make_pair(implications_.size(), "df"));
        add_comment("Theory for formulas df(j,k)");
        build_formulas_df();
        cout << "df: #implications=" << implications_.size() - imp_offsets_.back().first << " ((K-1) x (2L + L(L-1)))" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }

    // features added (forced)
    void force_feature(int j, int feature_index) {
        bool disable_ordering = (features_.size() == K_) || ((features_.size() < K_) && !actions_.empty());
        if( options_.ordering_ && !disable_ordering ) {
            Implication *IP = new Implication;
            for( int k = 0; k < K_; ++k )
                IP->add_consequent(1 + MappedBy(j, k));
            add_implication(IP);

            for( int k = 0; k < K_; ++k ) {
                for( int t = 0; t < log_num_features_; ++t ) {
                    Implication *IP = new Implication;
                    IP->add_antecedent(1 + MappedBy(j, k));
                    IP->add_consequent(bit(t, feature_index) ? (1 + F(k, t)) : -(1 + F(k, t)));
                    add_implication(IP);
                }
            }
        } else {
            for( int t = 0; t < log_num_features_; ++t ) {
                Implication *IP = new Implication;
                IP->add_consequent(bit(t, feature_index) ? 1 + F(j, t) : -(1 + F(j, t)));
                add_implication(IP);
            }
        }
    }
    void force_feature(int j, const string &feature) {
        int feature_index = M_.feature(feature);
        if( feature_index == -1 )
            cout << "error: inexistent feature '" << feature << "'" << endl;
        else
            force_feature(j, feature_index);
    }
    void add_implications_for_features() {
        assert(features_.size() <= K_);
        if( !features_.empty() ) {
            imp_offsets_.push_back(make_pair(implications_.size(), "forced-features"));
            add_comment("Implications for added (forced) features");
        }

        int j = 0;
        for( Items<string>::const_iterator it = features_.begin(); it != features_.end(); ++it, ++j )
            force_feature(j, *it);

        bool disable_ordering = (features_.size() == K_) || ((features_.size() < K_) && !actions_.empty());
        if( options_.ordering_ && !disable_ordering )
            cout << "features: #implications=" << implications_.size() - imp_offsets_.back().first << " (M x (1 + K x log2(#features)))" << endl;
        else
            cout << "features: #implications=" << implications_.size() - imp_offsets_.back().first << " (M x log2(#features))" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }
    void add_variables_for_features() {
        bool disable_ordering = (features_.size() == K_) || ((features_.size() < K_) && !actions_.empty());
        if( options_.ordering_ && !disable_ordering ) {
            var_offsets_.push_back(make_pair(variables_.size(), "MappedBy"));
            for( int j = 0; j < int(features_.size()); ++j ) {
                for( int k = 0; k < K_; ++k ) {
                    int index = variables_.size();
                    string name = string("MappedBy(") + to_string(j) + "," + to_string(k) + ")";
                    variables_.push_back(new Var(index, name));
                    assert(index == MappedBy(j, k));
                }
            }
            cout << "MappedBy: #variables=" << variables_.size() - var_offsets_.back().first << " (MK), offset=" << var_offsets_.back().first << endl;
        }
    }
    void read_features(const string &prefix) {
        string features_filename = filename(prefix, K_, N_, "_features.dat", false);
        cout << "reading '" << features_filename << "' ... " << flush;
        ifstream ifs(features_filename.c_str());
        if( !ifs.fail() ) {
            const Items<string> *features = Items<string>::read_dump(ifs);
            ifs.close();
            features_ = *features;
            delete features;
        } else {
            cout << "error: opening file '" << features_filename << "'" << endl;
        }
    }

    // CHECK: clauses to separate explicit goals from non-goals
    // Suggestion:
    //   1. Define variable not-sim(g,n,k) for goal state g and non-goal state n using feature k
    //   2. Require for each g, n, that there is f such that not-sim(g,n,k). Each requirement is
    //      a disjunction over k for not-sim(g,n,k)
    //   3. not-sim(g,n,k) need two implications:
    //      phi(g,f) & phi*(g,f) -> not-sim(g,n,k)
    //      phi*(g,f) & phi(g,f) -> not-sim(g,n,k)
    //   In total need G x N x L vars and G x N x (1 + 2K) implications where G=#goals and N=#non-goals
    void add_implications_for_goals() {
        assert(goals_.empty());
    }
    void add_variables_for_goals() {
        //CHECK assert(0);
    }
    void read_goals(const string &prefix) {
        string goals_filename = filename(prefix, K_, N_, "_goals.dat", false);
        cout << "reading '" << goals_filename << "' ... " << flush;
        ifstream ifs(goals_filename.c_str());
        if( !ifs.fail() ) {
            const Items<int> *goals = Items<int>::read_dump(ifs);
            ifs.close();
            goals_ = *goals;
            delete goals;
        } else {
            cout << "error: opening file '" << goals_filename << "'" << endl;
        }
    }

    // actions added (forced)
    void force_action(int j, const AbstractAction &action) {
        add_comment(string("Forced action ") + action.name_);
        for( int k = 0; k < K_; ++k ) {
            Implication *IP1 = new Implication;
            IP1->add_consequent(action.selected_effect_ & (1 << k) ? 1 + A(1, j, k) : -(1 + A(1, j, k)));
            add_implication(IP1);
            Implication *IP2 = new Implication;
            IP2->add_consequent(action.effect_ & (1 << k) ? 1 + A(2, j, k) : -(1 + A(2, j, k)));
            add_implication(IP2);
            Implication *IP3 = new Implication;
            IP3->add_consequent(action.selected_precondition_ & (1 << k) ? 1 + A(3, j, k) : -(1 + A(3, j, k)));
            add_implication(IP3);
            Implication *IP4 = new Implication;
            IP4->add_consequent(action.precondition_ & (1 << k) ? 1 + A(4, j, k) : -(1 + A(4, j, k)));
            add_implication(IP4);
        }
    }
    void add_implications_for_actions() {
        assert(actions_.size() <= N_);

        imp_offsets_.push_back(make_pair(implications_.size(), "force-actions"));
        add_comment("Theory for forced actions");

        int j = 0;
        for( Items<AbstractAction>::const_iterator it = actions_.begin(); it != actions_.end(); ++it, ++j )
            force_action(j, *it);

        cout << "actions: #implications=" << implications_.size() - imp_offsets_.back().first << " (?)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }
    void read_actions(const string &prefix); // defined below

    void enforce_soundness() {
        imp_offsets_.push_back(make_pair(implications_.size(), "soundness"));
        add_comment("Theory for enforcing soundness");
        for( int i = 0; i < num_states_; ++i ) {
            for( int j = 0; j < N_; ++j ) {
                Implication *IP = new Implication;

                // antecedent
                for( int k = 0; k < K_; ++k ) {
                    for( int t = 1; t < 3; ++t )
                        IP->add_antecedent(1 + B(t, j, i, k));
                }

                // consequent
                for( int l = 0; l < T_.num_transitions(i); ++l )
                    IP->add_consequent(1 + MATCH(i, l, j));

                // add implications
                add_implication(IP);
            }
        }
        cout << "soundness: #implications=" << implications_.size() - imp_offsets_.back().first << " (SN)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }

    void enforce_completeness() {
        imp_offsets_.push_back(make_pair(implications_.size(), "completeness"));
        add_comment("Theory for enforcing completeness");
        for( int i = 0; i < num_states_; ++i ) {
            for( int l = 0; l < T_.num_transitions(i); ++l ) {
                Implication *IP = new Implication;
                for( int j = 0; j < N_; ++j )
                    IP->add_consequent(1 + APPMATCH(i, l, j));
                add_implication(IP);
            }
        }
        cout << "completeness: #implications=" << implications_.size() - imp_offsets_.back().first << " (T)" << endl;
        //print_implications(cout, imp_offsets_.back().first, implications_.size());
    }

    // readers
    void read_minisat_output(ifstream &is) const {
        string status;
        is >> status;
        satisfiable_ = status == "SAT";
        if( satisfiable_ ) {
            int var, lit;
            model_ = vector<bool>(variables_.size(), false);
            for( size_t i = 0; i < variables_.size(); ++i ) {
                is >> lit;
                var = lit > 0 ? lit - 1 : -lit - 1;
                assert(var == int(i));
                model_[var] = lit > 0;
            }
            is >> lit;
            assert(lit == 0);
        } else {
            model_.clear();
        }
    }

    // output
    void dump(ostream &os) const {
        os << "p cnf " << variables_.size() << " " << implications_.size() << endl;
        size_t i = 0;
        for( size_t j = 0; j < implications_.size(); ++j ) {
            while( (i < comments_.size()) && (comments_[i].first == j) ) {
                os << "c " << comments_[i].second << endl;
                ++i;
            }
            implications_[j]->dump(os);
        }
        while( i < comments_.size() ) {
            os << "c " << comments_[i].second << endl;
            ++i;
        }
    }
    void print(ostream &os) const {
        size_t i = 0;
        for( size_t j = 0; j < implications_.size(); ++j ) {
            while( (i < comments_.size()) && (comments_[i].first == j) ) {
                os << "% " << comments_[i].second << endl;
                ++i;
            }
            implications_[j]->print(os, literals_);
            os << endl;
        }
        while( i < comments_.size() ) {
            os << "% " << comments_[i].second << endl;
            ++i;
        }
    }
    void print_implications(ostream &os, int start, int end) const {
        while( start < end ) {
            implications_[start++]->print(os, literals_);
            os << endl;
        }
    }

    void dump_model(ostream &os) const {
        for( size_t var = 0; var < variables_.size(); ++var ) {
            os << variables_[var] << " ";
        }
        os << "0" << endl;
    }
    void print_model(ostream &os) const {
        for( size_t var = 0; var < model_.size(); ++var ) {
            bool sign = model_[var];
            if( sign ) {
                assert(var < variables_.size());
                variables_[var]->print(os);
                os << endl;
            }
        }
    }

    // print abstraction coded in model
    void decode_model(const string &name, ostream &os, bool qnp_format) const {
        assert(satisfiable_ && (model_.size() == variables_.size()));

        // name
        os << name << endl;
        if( qnp_format )
            os << K_;
        else
            os << "Features:" << endl;

        // features
        for( int j = 0; j < K_; ++j ) {
            int feature_index = index_for_selected_feature(model_, j);
            assert(feature_index < num_features_);
            if( qnp_format )
                os << " " << M_.feature(feature_index) << " " << M_.numerical(feature_index);
            else
                os << "F" << j << ": " << M_.feature(feature_index) << endl;
        }
        if( qnp_format ) os << endl;

        // abstract actions
        set<AbstractAction> unique;
        for( int j = 0; j < N_; ++j ) {
            AbstractAction act(string("decoded-") + to_string(j));
            for( int k = 0; k < K_; ++k ) {
                act.selected_precondition_ += model_[A(3, j, k)] << k;
                act.precondition_ += model_[A(4, j, k)] << k;
                act.selected_effect_ += model_[A(1, j, k)] << k;
                act.effect_ += model_[A(2, j, k)] << k;
            }
            unique.insert(act);
        }

        if( qnp_format ) {
            os << unique.size() << endl;
            for( set<AbstractAction>::const_iterator it = unique.begin(); it != unique.end(); ++it )
                it->print(os, *this, model_, qnp_format);
        } else {
            os << "Abstract actions (#total=" << N_ << ", #unique=" << unique.size() << "):" << endl;
            for( set<AbstractAction>::const_iterator it = unique.begin(); it != unique.end(); ++it )
                it->print(os, *this, model_, qnp_format);
        }
    }
};

// static member
map<string, int> *Theory::AbstractAction::feature_map_ = nullptr;

istream& operator>>(istream &is, Theory::AbstractAction &act) {
    string name;
    is >> name;
    act.name_ = name;

    int nprec, neffect;
    for( is >> nprec; nprec > 0; --nprec ) {
        string feature;
        int value;
        is >> feature >> value;
        act.add_precondition(feature, value);
    }
    for( is >> neffect; neffect > 0; --neffect ) {
        string feature;
        int value;
        is >> feature >> value;
        act.add_effect(feature, value);
    }
  
    return is;
}

void Theory::read_actions(const string &prefix) {
    // construct and set feature map
    int j = 0;
    map<string, int> feature_map;
    for( Items<string>::const_iterator it = features_.begin(); it != features_.end(); ++it, ++j )
        feature_map.insert(make_pair(*it, j));
    Theory::AbstractAction::feature_map_ = &feature_map;

    // read file
    string actions_filename = filename(prefix, K_, N_, "_actions.dat", false);
    cout << "reading '" << actions_filename << "' ... " << flush;
    ifstream ifs(actions_filename.c_str());
    if( !ifs.fail() ) {
        const Items<AbstractAction> *actions = Items<AbstractAction>::read_dump(ifs);
        ifs.close();
        actions_ = *actions;
        delete actions;
    } else {
        cout << "error: opening file '" << actions_filename << "'" << endl;
    }

    // reset feature map
    Theory::AbstractAction::feature_map_ = nullptr;
}

pair<const Matrix*, const Transitions*> read_data(const string &matrix_filename, const string &transitions_filename) {
    pair<const Matrix*, const Transitions*> p(nullptr, nullptr);

    cout << "reading '" << matrix_filename << "' ... " << flush;
    ifstream ifs_matrix(matrix_filename.c_str());
    if( !ifs_matrix.fail() ) {
        p.first = Matrix::read_dump(ifs_matrix);
        ifs_matrix.close();
    } else {
            throw std::runtime_error("error: opening file '" + matrix_filename + "'");
        cout << "error: opening file '" << matrix_filename << "'" << endl;
    }

    cout << "reading '" << transitions_filename << "' ... " << flush;
    ifstream ifs_transitions(transitions_filename.c_str());
    if( !ifs_transitions.fail() ) {
        p.second = Transitions::read_dump(ifs_transitions);
        ifs_transitions.close();
    } else {
        cout << "error: opening file '" << transitions_filename << "'" << endl;
    }

    return p;
}

void usage(ostream &os, const string &name) {
    cout << endl
         << "usage: " << name << " [--add-actions] [--add-features] [--decode] [--explicit-goals] [--only-soundness] <prefix> <K> <N>" << endl
         << endl
         << "where" << endl
         << "    <prefix> is prefix for all files" << endl
         << "    K is numer of features to select" << endl
         << "    N is number of abstract actions." << endl
         << endl
         << "For the options," << endl
         << "    --add-actions to force inclusion of actions in '<prefix>_actions.dat'" << endl
         << "    --add-features to force inclusion of features in '<prefix>_features.dat'" << endl
         << "    --decode to decode model in '<prefix>_<K>_<N>_model.cnf' found by minisat" << endl
         << "    --explicit-goals to separate goals from non-goals using goals in '<prefix>_goals.dat'" << endl
         << "    --only-soundness to generate clauses that only enforce soundness instead of soundness+completeness" << endl
         << endl
         << "A feature file consists of number <n> of features in file, followed by space" << endl
         << "followed by space-separated list of features names (the names need to match those" << endl
         << "in matrix file." << endl
         << endl
         ;
}

int main(int argc, const char **argv) {
    // print call
    cout << "call:";
    for( int i = 0; i < argc; ++i )
        cout << " " << argv[i];
    cout << endl;

    // read executable name
    string name(*argv);

    // read options
    Theory::Options options;
    for( ++argv, --argc; (argc > 0) && (**argv == '-'); ++argv, --argc ) {
        if( string(*argv) == "--decode" ) {
            options.decode_ = true;
        } else if( string(*argv) == "--explicit-goals" ) {
            options.explicit_goals_ = true;
        } else if( string(*argv) == "--add-actions" ) {
            options.add_actions_ = true;
        } else if( string(*argv) == "--add-features" ) {
            options.add_features_ = true;
        } else if( string(*argv) == "--only-soundness" ) {
            options.only_soundness_ = true;
        } else {
            cout << "error: unrecognized option '" << *argv << "'" << endl;
            usage(cout, name);
            exit(0);
        }
    }

    // check we have enough arguments
    if( argc < 3 ) {
        usage(cout, name);
        exit(0);
    }

    // read arguments
    options.prefix_ = argv[0];
    int K = atoi(argv[1]);
    int N = atoi(argv[2]);
    cout << "arguments: prefix=" << options.prefix_ << ", K=" << K << ", N=" << N << endl;

    // input filenames
    string matrix_filename = filename(options.prefix_, K, N, "_matrix.dat", false);
    string transitions_filename = filename(options.prefix_, K, N, "_transitions.dat", false);

    // read data
    pair<const Matrix*, const Transitions*> p;
    p = read_data(matrix_filename, transitions_filename);
    const Matrix *M = p.first;
    const Transitions *Tr = p.second;

    // create theory
    Theory Th(*M, *Tr, K, N, options);
    cout << "#variables=" << Th.num_variables() << endl;
    cout << "#implications=" << Th.num_implications() << endl;

    // decode
    if( options.decode_ ) {
        string model_filename = filename(options.prefix_, K, N, "_model.cnf");
        cout << "reading file '" << model_filename << "' ..." << flush;
        ifstream ifs(model_filename.c_str());
        if( !ifs.fail() ) {
            Th.read_minisat_output(ifs);
            ifs.close();
            cout << " done!" << endl;

            // output model in qnp format
            string qnp_filename = filename(options.prefix_, K, N, ".qnp");
            cout << "writing file '" << qnp_filename << "' ..." << flush;
            ofstream os(qnp_filename.c_str());
            Th.decode_model(model_filename, os, true);
            os.close();
            cout << " done!" << endl;
        } else {
            cout << "error: opening file '" << model_filename << "'" << endl;
        }
    } else {
        // output theory in human readable format
        //Th.print(cout);

        // output theory in SATLIB format
        string theory_filename = filename(options.prefix_, K, N, "_theory.cnf");
        cout << "writing file '" << theory_filename << "' ..." << flush;
        ofstream os(theory_filename.c_str());
        Th.dump(os);
        os.close();
        cout << " done!" << endl;
    }

    delete M;
    delete Tr;
    return 0;
}

