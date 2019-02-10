#ifndef FEATURES_HXX
#define FEATURES_HXX

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cassert>

//#include <utility>
//#include <cstdint>
//#include <vector>
//#include <stdexcept>

namespace SLTP {

namespace DL {

class GroundedPredicate;

// A state denotation for concept C is a subset of objects
// implemented as a bitmap (ie. vector<bool>). These will
// be cached so that at most one copy of such a bitmap
// exists.
struct state_denotation_t : public std::vector<bool> {
    state_denotation_t(size_t n, bool value) : std::vector<bool>(n, value) { }
};

// A (full) sample denotation for concept C is a vector of
// state denotations, one per each state in the sample.
// Since we will cache state denotation, a sample denotation
// is implemented as a vector of pointers to bitmaps.
struct sample_denotation_t : public std::vector<const state_denotation_t*> { };

// We cache sample and state denotations. The latter are cached
// into a simple hash (i.e. unordered_set). The former are cached
// using two hash maps (i.e. unordered_map): one that maps sample
// denotations to concept names, and the other that maps concept
// names to sample denotations.
//
// We also cache grounded predicates that are used to 
// represent states in the sample.
class Cache {
  public:
    struct cache_support_t {
        bool operator()(const sample_denotation_t *d1, const sample_denotation_t *d2) const {
            assert(d1->size() == d2->size());
            for( int i = 0; i < int(d1->size()); ++i ) {
                if( d1[i] != d2[i] )
                    return false;
            }
            return true;
        }
        size_t operator()(const sample_denotation_t *obj) const {
            return 0;
        }

        bool operator()(const state_denotation_t *sd1, const state_denotation_t *sd2) const {
            assert(sd1->size() == sd2->size());
            return *sd1 == *sd2;
        }
        size_t operator()(const state_denotation_t *obj) const {
            return 0;
        }
    };

    typedef typename std::unordered_map<const sample_denotation_t*, std::string, cache_support_t, cache_support_t> cache1_t;
    typedef typename std::unordered_map<std::string, const sample_denotation_t*> cache2_t;
    typedef typename std::unordered_set<const state_denotation_t*, cache_support_t, cache_support_t> cache3_t;
    typedef typename std::unordered_map<std::string, const GroundedPredicate*> cache4_t;

  protected:
    cache1_t cache1_;
    cache2_t cache2_;
    cache3_t cache3_;
    cache4_t cache4_;

  public:
    Cache() { }
    ~Cache() { }

    // cache1: (full) sample denotations for concepts
    //
    // We maintain a hash of (full) sample denotations so that one copy
    // exists for each such denotations. This hash is implemented with
    // two tables, one that provides mappings from sample denotations
    // to concept names, and the other is the inverse
    //
    // sample denotation -> concept name
    const sample_denotation_t* find_sample_denotation(const sample_denotation_t &d) const {
        cache1_t::const_iterator it = cache1_.find(&d);
        return it == cache1_.end() ? nullptr : it->first;
    }
    const sample_denotation_t* find_or_insert_sample_denotation(const sample_denotation_t &d);

    // cache2: (full) sample denotations for concepts
    //
    // concept name -> sample denotation
    const sample_denotation_t* find_sample_denotation(const std::string name) const {
        cache2_t::const_iterator it = cache2_.find(name);
        return it == cache2_.end() ? nullptr : it->second;
    }

    // cache3: state denotations (bit vectors)
    //
    // We maintain a hash of these denotations such as to keep just one
    // copy of each such bitmap. A sample denotation is a vector of
    // state denotations. The same bitmap may appear associated with
    // different states in the same sample, or across sample denotations
    // for different concepts
    //
    // state denotation -> pointer
    const state_denotation_t* find_state_denotation(const state_denotation_t &sd) const {
        cache3_t::const_iterator it = cache3_.find(&sd);
        return it == cache3_.end() ? nullptr : *it;
    }

    const state_denotation_t* find_or_insert_state_denotation(const state_denotation_t &sd) {
        cache3_t::const_iterator it = cache3_.find(&sd);
        if( it == cache3_.end() ) {
            const state_denotation_t *nsd = new state_denotation_t(sd);
            cache3_.insert(nsd);
            return nsd;
        } else {
            return *it;
        }
    }

    // cache4: grounded predicates
    //
    // name -> grounded predicate
    const GroundedPredicate* find_grounded_predicate(const std::string &name) const {
        cache4_t::const_iterator it = cache4_.find(name);
        return it == cache4_.end() ? nullptr : it->second;
    }
    const GroundedPredicate* find_or_insert_grounded_predicate(const GroundedPredicate &gp);
};

// We represent states as subets of grounded predicates.
// A grounded predicate is a predicate and a vector of
// objects to ground the predicates.
struct Object {
    const std::string name_;
    Object(const std::string name) : name_(name) { }
    const std::string& as_str() const {
        return name_;
    }
    friend std::ostream& operator<<(std::ostream &os, const Object &obj) {
        return os << obj.as_str() << std::flush;
    }
};

struct Predicate {
    const std::string name_;
    int arity_;
    Predicate(const std::string &name, int arity) : name_(name), arity_(arity) { }
    std::string as_str(const std::vector<const Object*> *objects) const {
        std::string str = name_ + "(";
        if( objects == nullptr ) {
            for( int i = 0; i < arity_; ++i ) {
                str += std::string("x") + std::to_string(1 + i);
                if( 1 + i < arity_ ) str += ",";
            }
        } else {
            assert(objects->size() == arity_);
            for( int i = 0; i < arity_; ++i ) {
                str += (*objects)[i]->as_str();
                if( 1 + i < arity_ ) str += ",";
            }
        }
        return str + ")";
    }
    friend std::ostream& operator<<(std::ostream &os, const Predicate &pred) {
        return os << pred.as_str(nullptr) << std::flush;
    }
};

struct GroundedPredicate {
    const Predicate &predicate_;
    const std::vector<const Object*> objects_;
    GroundedPredicate(const Predicate &predicate, const std::vector<const Object*> &objects)
      : predicate_(predicate), objects_(objects) {
    }
    GroundedPredicate(const Predicate &predicate, std::vector<const Object*> &&objects)
      : predicate_(predicate), objects_(std::move(objects)) {
    }
    std::string as_str() const {
        return predicate_.as_str(&objects_);
    }
    friend std::ostream& operator<<(std::ostream &os, const GroundedPredicate &pred) {
        return os << pred.as_str() << std::flush;
    }
};

struct State {
    const int id_;
    std::vector<const GroundedPredicate*> atoms_;
    State(int id) : id_(id) { }
};

class Sample {
  protected:
    const std::string name_;
    const std::vector<Object> objects_;
    const std::vector<const Predicate*> predicates_;
    const std::vector<const GroundedPredicate*> grounded_predicates_;
    const std::vector<State> states_;

    Sample(std::string name, std::vector<Object> objects,
           std::vector<const Predicate*> predicates,
           std::vector<const GroundedPredicate*> grounded_predicates,
           const std::vector<State>& states) :
        name_(std::move(name)),
        objects_(std::move(objects)),
        predicates_(std::move(predicates)),
        grounded_predicates_(std::move(grounded_predicates)),
        states_(std::move(states))
    { }

public:
    ~Sample() {
      for (auto grounded_predicate : grounded_predicates_) delete grounded_predicate;
      for (auto predicate : predicates_) delete predicate;
    }

    const std::string& name() const {
        return name_;
    }
    std::size_t num_objects() const {
        return objects_.size();
    }
    std::size_t num_states() const {
        return states_.size();
    }
    const std::vector<State>& states() const {
        return states_;
    }

    const State& state(unsigned id) const { return states_.at(id); }

    static Sample read(std::istream &is);
};

struct Denotation {
    enum class denotation_t : bool { concept_denotation_t, role_denotation_t };
    denotation_t type_;
    std::vector<bool> values_;
    Denotation(denotation_t type, size_t dimension)
      : type_(type),
        values_(std::vector<bool>(dimension, false)) {
    }
};

// Denotation matrix is just a matrix of DL Denotations, each state being a row, each concept / role a column
class DenotationMatrix {
  protected:
    std::vector<std::vector<Denotation> > data_;
};

class Model {
};

class Base {
  protected:
    // unsigned id_;
    unsigned complexity_;

  public:
    explicit Base(unsigned complexity) : complexity_(complexity) { }
    // const unsigned id() const { return id_; }
    const unsigned complexity() const { return complexity_; }

    //virtual void denotation(Denotation &d) const = 0;
    virtual const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const = 0;
    virtual const state_denotation_t* denotation(Cache &cache, const State &state) const = 0;
    virtual std::string as_str() const = 0;

    friend std::ostream& operator<<(std::ostream &os, const Base &base) {
        return os << base.as_str() << std::flush;
    }
};

class Concept : public Base {
  public:
    explicit Concept(unsigned complexity) : Base(complexity) { }
    virtual ~Concept() { }
};

class Role : public Base {
  public:
    explicit Role(unsigned complexity) : Base(complexity) { }
    virtual ~Role() { }
};

class PrimitiveConcept : public Concept {
  protected:
    const std::string name_;

  public:
    PrimitiveConcept(const std::string &name) : Concept(1), name_(name) { }
    virtual ~PrimitiveConcept() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            sample_denotation_t d;
            for( unsigned i = 0; i < sample.num_states(); ++i ) {
                const State &s = sample.state(i);
                d.emplace_back(denotation(cache, s));
            }
            return cache.find_or_insert_sample_denotation(d);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("TODO: UNIMPLEMENTED: PrimitiveConcept::denotation()");
    }
    std::string as_str() const override {
        return name_;
    }
};

class UniversalConcept : public Concept {
  public:
    UniversalConcept() : Concept(1) { }
    virtual ~UniversalConcept() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            int num_objects = sample.num_objects();
            state_denotation_t sd(num_objects, true);
            const state_denotation_t *cached_sd = cache.find_or_insert_state_denotation(sd);

            sample_denotation_t d;
            for( int i = 0; i < int(sample.num_states()); ++i )
                d.emplace_back(cached_sd);
            return cache.find_or_insert_sample_denotation(d);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("Unexpected call: UniversalConcept::denotation(Cache&, const State&)");
        return nullptr;
    }
    std::string as_str() const override {
        return "<universal>";
    }
};

class AndConcept : public Concept {
  protected:
    const Concept *concept1_;
    const Concept *concept2_;

  public:
    AndConcept(const Concept *concept1, const Concept *concept2)
      : Concept(1 + concept1->complexity() + concept2->complexity()),
        concept1_(concept1),
        concept2_(concept2) {
    }
    virtual ~AndConcept() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            const sample_denotation_t *d1 = cache.find_sample_denotation(concept1_->as_str());
            assert(d1 != nullptr);
            const sample_denotation_t *d2 = cache.find_sample_denotation(concept2_->as_str());
            assert(d2 != nullptr);

            sample_denotation_t d;
            assert(d1->size() == d2->size());
            for( int i = 0; i < int(d1->size()); ++i ) {
                const state_denotation_t &sd1 = *(*d1)[i];
                const state_denotation_t &sd2 = *(*d2)[i];
                assert(sd1.size() == sd2.size());

                state_denotation_t sd(sd1.size(), false);
                for( int j = 0; j < int(sd1.size()); ++j ) {
                    if( sd1[j] && sd2[j] )
                        sd[j] = true;
                }

                const state_denotation_t *cached_sd = cache.find_or_insert_state_denotation(sd);
                d.emplace_back(cached_sd);
            }
            return cache.find_or_insert_sample_denotation(d);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("Unexpected call: AndConcept::denotation(Cache&, const State&)");
        return nullptr;
    }

    std::string as_str() const override {
        return std::string("And[") +  concept1_->as_str() + "," + concept2_->as_str() + "]";
    }
};

class NotConcept : public Concept {
  protected:
    const Concept *concept_;

  public:
    NotConcept(const Concept *concept)
      : Concept(1 + concept->complexity()),
        concept_(concept) {
    }
    virtual ~NotConcept() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            const sample_denotation_t *d = cache.find_sample_denotation(concept_->as_str());
            assert(d != nullptr);

            sample_denotation_t nd;
            for( int i = 0; i < int(d->size()); ++i ) {
                const state_denotation_t &sd = *(*d)[i];
                state_denotation_t nsd(sd.size(), false);
                for( int j = 0; j < int(sd.size()); ++j )
                    nsd[j] = !sd[j];
                const state_denotation_t *cached_sd = cache.find_or_insert_state_denotation(nsd);
                nd.emplace_back(cached_sd);
            }
            return cache.find_or_insert_sample_denotation(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("Unexpected call: NotConcept::denotation(Cache&, const State&)");
        return nullptr;
    }

    std::string as_str() const override {
        return std::string("Not[") + concept_->as_str() + "]";
    }
};

class ExistsConcept : public Concept {
  protected:
    const Concept *concept_;
    const Role *role_;

  public:
    ExistsConcept(const Concept *concept, const Role *role)
      : Concept(1 + concept->complexity() + role->complexity()),
        concept_(concept),
        role_(role) {
    }
    virtual ~ExistsConcept() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            const sample_denotation_t *d = cache.find_sample_denotation(concept_->as_str());
            assert((d != nullptr) && (d->size() == sample.num_states()));
            const sample_denotation_t *r = cache.find_sample_denotation(role_->as_str());
            assert((r != nullptr) && (r->size() == sample.num_states() * sample.num_states()));

            sample_denotation_t nd;
            for( unsigned i = 0; i < sample.num_states(); ++i ) {
                const State &s = sample.state(i);
                state_denotation_t sd(sample.num_objects(), false);
                for( int j = 0; j < sample.num_objects(); ++j ) {
                    for( int k = 0; !sd[j] && (k < sample.num_objects()); ++k ) {
                        int index = j * sample.num_objects() + k;
                        sd[j] = (*(*r)[i])[index];
                    }
                }
                nd.emplace_back(cache.find_or_insert_state_denotation(sd));
            }
            return cache.find_or_insert_sample_denotation(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("Unexpected call: ExistsConcept::denotation(Cache&, const State&)");
        return nullptr;
    }

    std::string as_str() const override {
        return std::string("Exists[") + concept_->as_str() + "," + role_->as_str() + "]";
    }
};

class ForallConcept : public Concept {
  protected:
    const Concept *concept_;
    const Role *role_;

  public:
    ForallConcept(const Concept *concept, const Role *role)
      : Concept(1 + concept->complexity() + role->complexity()),
        concept_(concept),
        role_(role) {
    }
    virtual ~ForallConcept() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            const sample_denotation_t *d = cache.find_sample_denotation(concept_->as_str());
            assert((d != nullptr) && (d->size() == sample.num_states()));
            const sample_denotation_t *r = cache.find_sample_denotation(role_->as_str());
            assert((r != nullptr) && (r->size() == sample.num_states() * sample.num_states()));

            sample_denotation_t nd;
            for( unsigned i = 0; i < sample.num_states(); ++i ) {
                const State &s = sample.state(i);
                state_denotation_t sd(sample.num_objects(), true);
                for( int j = 0; j < sample.num_objects(); ++j ) {
                    for( int k = 0; sd[j] && (k < sample.num_objects()); ++k ) {
                        int index = j * sample.num_objects() + k;
                        sd[j] = (*(*r)[i])[index];
                    }
                }
                nd.emplace_back(cache.find_or_insert_state_denotation(sd));
            }
            return cache.find_or_insert_sample_denotation(nd);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("Unexpected call: ForallConcept::denotation(Cache&, const State&)");
        return nullptr;
    }

    std::string as_str() const override {
        return std::string("Forall[") + concept_->as_str() + "," + role_->as_str() + "]";
    }
};

class PrimitiveRole : public Role {
  protected:
    const std::string name_;

  public:
    PrimitiveRole(const std::string &name) : Role(1), name_(name) { }
    virtual ~PrimitiveRole() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            sample_denotation_t d;
            for( unsigned i = 0; i < sample.num_states(); ++i ) {
                const State &s = sample.state(i);
                d.emplace_back(denotation(cache, s));
            }
            return cache.find_or_insert_sample_denotation(d);
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("TODO: UNIMPLEMENTED: PrimitiveRole::denotation()");
    }

    std::string as_str() const override {
        return name_;
    }
};

class PlusRole : public Role {
  protected:
    const Role *role_;

  public:
    PlusRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }
    virtual ~PlusRole() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            throw std::runtime_error("TODO: UNIMPLEMENTED: PlusRole::denotation()");
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("TODO: UNIMPLEMENTED: PlusRole::denotation()");
    }

    std::string as_str() const override {
        return std::string("Plus[") + role_->as_str() + "]";
    }
};

class StarRole : public Role {
  protected:
    const Role *role_;

  public:
    StarRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }
    virtual ~StarRole() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            throw std::runtime_error("TODO: UNIMPLEMENTED: StarRole::denotation()");
        } else {
            return cached;
        }
    }
    const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("TODO: UNIMPLEMENTED: StarRole::denotation()");
    }

    std::string as_str() const override {
        return std::string("Star[") + role_->as_str() + "]";
    }
};

class InverseRole : public Role {
  protected:
    const Role *role_;

  public:
    InverseRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }
    virtual ~InverseRole() { }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const sample_denotation_t *cached = cache.find_sample_denotation(as_str());
        if( cached == nullptr ) {
            const sample_denotation_t *r = cache.find_sample_denotation(role_->as_str());
            assert((r != nullptr) && (r->size() == sample.num_states() * sample.num_states()));

            sample_denotation_t nd;
            for( unsigned i = 0; i < sample.num_states(); ++i ) {
                const State &s = sample.state(i);
                state_denotation_t sd(sample.num_objects() * sample.num_objects(), false);
                for( int j = 0; j < sample.num_objects(); ++j ) {
                    for( int k = 0; k < sample.num_objects(); ++k ) {
                        int index = j * sample.num_objects() + k;
                        if( (*(*r)[i])[index] ) {
                            int inv_index = k * sample.num_objects() + j;
                            sd[inv_index] = true;
                        }
                    }
                }
                nd.emplace_back(cache.find_or_insert_state_denotation(sd));
            }
            return cache.find_or_insert_sample_denotation(nd);
        } else {
            return cached;
        }
    }
    virtual const state_denotation_t* denotation(Cache &cache, const State &state) const override {
        throw std::runtime_error("Unexpected call: InverseRole::denotation(Cache&, const State&)");
        return nullptr;
    }

    std::string as_str() const override {
        return std::string("Inverse[") + role_->as_str() + "]";
    }
};

class Factory {
  protected:
    const std::string name_;
    std::vector<const Role*> basis_roles_;
    std::vector<const Concept*> basis_concepts_;

    int complexity_bound_;

    mutable int current_complexity_;
    mutable std::vector<const Role*> roles_;
    mutable std::vector<std::vector<const Concept*> > concepts_;

  public:
    Factory(const std::string &name, int complexity_bound)
      : name_(name),
        complexity_bound_(complexity_bound),
        current_complexity_(-1) {
    }
    virtual ~Factory() { }

    const std::string& name() const {
        return name_;
    }

    void insert_basis(const Role *role) {
        basis_roles_.push_back(role);
    }
    void insert_basis(const Concept *concept) {
        basis_concepts_.push_back(concept);
    }

    int current_complexity() const {
        return current_complexity_;
    }
    void set_complexity_bound(int complexity_bound) {
        complexity_bound_ = complexity_bound;
    }

    int reset(bool remove) {
        while( roles_.size() > 0 ) {
            if( remove && (roles_.size() >= basis_roles_.size()) )
                delete roles_.back();
            roles_.pop_back();
        }

        while( concepts_.size() > 1 ) {
            if( remove ) {
                for( int k = 0; k < int(concepts_.back().size()); ++k )
                    delete concepts_.back()[k];
            }
            concepts_.pop_back();
        }
        concepts_.pop_back();

        current_complexity_ = -1;
        return current_complexity_;
    }

    int advance_step() const {
        if( current_complexity_ == -1 ) {
            assert(concepts_.empty());
            concepts_.emplace_back(basis_concepts_);
        } else {
            std::vector<const Concept*> concepts;
            for( int i = 0; i < int(concepts_.back().size()); ++i ) {
                const Concept *c = concepts_.back()[i];
                concepts.push_back(new NotConcept(c));
                for( int j = 0; j < int(roles_.size()); ++j ) {
                    const Role *r = roles_[j];
                    concepts.push_back(new ExistsConcept(c, r));
                    concepts.push_back(new ForallConcept(c, r));
                }
                for( int j = 1 + i; j < int(concepts_.back().size()); ++j ) {
                    const Concept *oc = concepts_.back()[j];
                    concepts.push_back(new AndConcept(c, oc));
                }
            }
            concepts_.emplace_back(std::move(concepts));
        }
        ++current_complexity_;
        return current_complexity_;
    }

    bool is_superfluous(Cache &cache, const sample_denotation_t *d) const {
        return cache.find_sample_denotation(*d);
    }
    void insert_new_denotation(Cache &cache, const Concept &concept, const sample_denotation_t *d) const {
        //cache.emplace(denotation, concept.as_str());
    }

    int prune_superfluous_concepts_in_last_layer(Cache &cache, const Sample &sample) const {
        if( concepts_.empty() ) return 0;

        // for each concept, check whether there is another concept with
        // same denotation in all states. For this, create a vector of
        // denotations, one for each state, and check whether there is
        // another concept with same vector. We store vectors in a hash
        // table to make the check more efficient

        int num_pruned_concepts = 0;
        for( int i = 0; i < int(concepts_.back().size()); ++i ) {
            const Concept &c = *concepts_.back()[i];
            const sample_denotation_t *d = c.denotation(cache, sample);
            if( !is_superfluous(cache, d) ) {
                insert_new_denotation(cache, c, d);
            } else {
                // make sure we do not eliminate a basis
                if( concepts_.size() > 1 ) delete concepts_.back()[i];
                concepts_.back()[i] = concepts_.back().back();
                concepts_.back().pop_back();
                --i;
                ++num_pruned_concepts;
            }
        }
        return num_pruned_concepts;
    }

    int generate_roles() const {
        if( roles_.empty() ) {
            roles_.insert(roles_.end(), basis_roles_.begin(), basis_roles_.end());
            for( int i = 0; i < int(basis_roles_.size()); ++i ) {
                const Role *r = basis_roles_[i];
                roles_.push_back(new PlusRole(r));
                roles_.push_back(new StarRole(r));
                roles_.push_back(new InverseRole(r));
                const Role *ir = roles_.back();
                roles_.push_back(new PlusRole(ir));
                roles_.push_back(new StarRole(ir));
            }
        }
        return roles_.size();
    }

    int generate_concepts(Cache &cache, const Sample *sample = nullptr, bool prune = false) const {
        while( current_complexity_ < complexity_bound_ ) {
            std::cout << "DL::Factory: name="
                      << name_
                      << ", complexity="
                      << 1 + current_complexity_
                      << std::flush;

            advance_step();
            if( prune && (sample != nullptr) ) {
                int number_pruned_concepts = prune_superfluous_concepts_in_last_layer(cache, *sample);
                std::cout << ", #pruned-concepts="
                          << number_pruned_concepts
                          << std::flush;
            }

            std::cout << ", #concepts="
                      << concepts_.back().size()
                      << std::endl;
        }
        return current_complexity_;
    }
};

}; // DL namespace

}; // SLTP namespace

#endif

