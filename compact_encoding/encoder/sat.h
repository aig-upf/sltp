#ifndef SAT_H
#define SAT_H

#include <bitset>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

namespace SAT {

// returns k-th bit in integer l
inline bool bit(int k, int l) {
    return l & (1 << k);
}

// returns a string with binary representation of integer k (16 bits)
std::string bit_string(int k) {
    return std::bitset<8>(k).to_string();
}

class Var {
  protected:
    int index_;
    std::string str_;

  public:
    Var(int index, const std::string &str) : index_(index), str_(str) { }
    virtual ~Var() { }
    int index() const { return index_; }
    std::string str() const { return str_; }
    void print(std::ostream &os) const {
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
    std::string as_str() const {
        return !sign_ ? var_.str() : std::string("-") + var_.str();
    }
    void print(std::ostream &os) const {
        os << as_str();
    }
};

class Implication {
  protected:
    std::vector<int> antecedent_; // joined by AND
    std::vector<int> consequent_; // joined by OR

  public:
    Implication() { }
    ~Implication() { }

    void add_antecedent(int L) { antecedent_.push_back(L); }
    void add_consequent(int L) { consequent_.push_back(L); }

    void dump(std::ostream &os) const {
        for( size_t i = 0; i < antecedent_.size(); ++i ) {
            os << -antecedent_[i];
            if( i + 1 < antecedent_.size() ) os << " ";
        }
        if( !consequent_.empty() ) os << " ";
        for( size_t i = 0; i < consequent_.size(); ++i ) {
            os << consequent_[i];
            if( i + 1 < consequent_.size() ) os << " ";
        }
        os << " 0" << std::endl;
    }
    void print(std::ostream &os, const std::vector<const Literal*> &pos_literals, const std::vector<const Literal*> &neg_literals) const {
        for( size_t i = 0; i < antecedent_.size(); ++i ) {
            if( antecedent_[i] > 0 )
                pos_literals[antecedent_[i] - 1]->print(os);
            else
                neg_literals[-antecedent_[i] - 1]->print(os);
            if( i + 1 < antecedent_.size() ) os << " & ";
        }
        os << " => ";
        for( size_t i = 0; i < consequent_.size(); ++i ) {
            if( consequent_[i] > 0 )
                pos_literals[consequent_[i] - 1]->print(os);
            else
                neg_literals[-consequent_[i] - 1]->print(os);
            if( i + 1 < consequent_.size() ) os << " v ";
        }
        os << std::flush;
    }
};

}; // namespace SAT

#endif

