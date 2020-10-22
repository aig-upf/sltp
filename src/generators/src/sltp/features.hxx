
#pragma once

#include <cassert>
#include <deque>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <limits>
#include <algorithm>
#include <fstream>
#include <boost/bimap.hpp>

#include <sltp/base.hxx>
#include <sltp/utils.hxx>
#include <sltp/algorithms.hxx>
#include <common/base.h>

#include <ctime>

namespace sltp {
    class TransitionSample;
}

namespace sltp {
    class Atom;
    class Sample;
    class State;
}

namespace sltp::dl {

const unsigned PRIMITIVE_COMPLEXITY = 1;

class DLBaseElement;
class Concept;
class Role;
class Feature;

//! A cache from features to their sample denotations
using feature_cache_t = std::unordered_map<
        feature_sample_denotation_t, const Feature*, utils::container_hash<feature_sample_denotation_t> >;


//! Command-line option processing
struct Options {
    std::string workspace;
    int timeout;
    unsigned complexity_bound;
    unsigned dist_complexity_bound;
    unsigned cond_complexity_bound;
    bool comparison_features;
    bool generate_goal_concepts;
    bool print_denotations;
};

// We cache sample and state denotations. The latter are cached
// into a simple hash (i.e. unordered_set). The former are cached
// using two hash maps (i.e. unordered_map): one that maps sample
// denotations to concept names, and the other that maps concept
// names to sample denotations.
//
// We also cache atoms that are used to represent states in the sample.
class Cache {
public:
    struct cache_support_t {
        // cache for sample denotations
        bool operator()(const sample_denotation_t *d1, const sample_denotation_t *d2) const {
            assert((d1 != nullptr) && (d2 != nullptr));
            assert(d1->size() == d2->size()); // number of states is fixed in sample
            return *d1 == *d2;
        }

        size_t operator()(const sample_denotation_t *obj) const {
            assert(obj != nullptr);
            size_t hash = (*this)((*obj)[0]);
            for (auto elem:*obj)
                hash = hash ^ (*this)(elem);
            return hash;
        }

        // cache for state denotations
        bool operator()(const state_denotation_t *sd1, const state_denotation_t *sd2) const {
            assert((sd1 != nullptr) && (sd2 != nullptr));
            // dimension of sd1 and sd2 may not be equal since they may
            // be for states with different number of objects, or for
            // denotation of concept/roles
            return *sd1 == *sd2;
        }

        size_t operator()(const state_denotation_t *obj) const {
            assert(obj != nullptr);
            std::hash<std::vector<bool> > hasher;
            return hasher(*static_cast<const std::vector<bool>*>(obj));
        }
    };

    //! A map from sample denotations to the concept/role ID with that denotation
    using cache1_t = std::unordered_map<const sample_denotation_t*, unsigned long, cache_support_t, cache_support_t>;

    //! A map from concept/role IDs to their denotation
    using cache2_t = std::unordered_map<unsigned long, const sample_denotation_t*>;

    //! A set of state denotations that we'll keep as a register, to avoid having duplicate state_denotation_t
    //! objects around.
    using cache3_t = std::unordered_set<const state_denotation_t*, cache_support_t, cache_support_t>;

protected:
    cache1_t cache1_;
    cache2_t cache2_;
    cache3_t cache3_;

public:
    Cache() = default;
    ~Cache() = default;

    // cache1: (full) sample denotations for concepts
    //
    // We maintain a hash of (full) sample denotations so that one copy
    // exists for each such denotations. This hash is implemented with
    // two tables, one that provides mappings from sample denotations
    // to concept names, and the other is the inverse
    //
    // sample denotation -> concept name
    const cache1_t& cache1() const {
        return cache1_;
    }

    bool contains(const sample_denotation_t* d) const {
        return cache1_.find(d) != cache1_.end();
    }

    const sample_denotation_t* find_or_insert_sample_denotation(const sample_denotation_t &d, unsigned long id) {
        auto it = cache1_.find(&d);
        if( it == cache1_.end() ) {
            assert(cache2_.find(id) == cache2_.end());
            const sample_denotation_t *nd = new sample_denotation_t(d);
            cache1_.emplace(nd, id);
            cache2_.emplace(id, nd);
            return nd;
        } else {
            return it->first;
        }
    }

    // Return the sample denotation corresponding to the DL element with given ID
    const sample_denotation_t& find_sample_denotation(const DLBaseElement& element, std::size_t expected_size) const;

    // cache3: state denotations (bit vectors)
    //
    // We maintain a hash of these denotations such as to keep just one
    // copy of each such bitmap. A sample denotation is a vector of
    // state denotations. The same bitmap may appear associated with
    // different states in the same sample, or across sample denotations
    // for different concepts
    //
    // state denotation -> pointer
    const state_denotation_t* find_or_insert_state_denotation(const state_denotation_t &sd) {
        auto it = cache3_.find(&sd);
        if( it == cache3_.end() ) {
            const state_denotation_t *nsd = new state_denotation_t(sd);
            cache3_.insert(nsd);
            return nsd;
        } else {
            return *it;
        }
    }

    const state_denotation_t& retrieveDLDenotation(
            const DLBaseElement& element, const State &state, std::size_t expected_size) const;
};

//! A common base class for concepts and roles
class DLBaseElement {
private:
    static unsigned long global_id;

protected:
    const unsigned long id_;
    int complexity_;

public:
    explicit DLBaseElement(int complexity) : id_(global_id++), complexity_(complexity) { }

    [[nodiscard]] int complexity() const { return complexity_; }

    [[nodiscard]] unsigned long id() const { return id_; }

    //! By default we raise an exception here, as we want users to use the method that returns the denotation
    //! over the whole sample, which will be more efficient. Some subclasses though will override this, mostly
    //! Primitive concepts and roles.
    virtual const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const {
        throw std::runtime_error("Unexpected call to DLBaseElement::denotation(Cache&, const Sample&, const State&)");
    }

    //! Compute the full sample denotation for the current DL element
    virtual const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const = 0;

    //! Return a string representation of the concept or role
    [[nodiscard]] virtual std::string as_str() const = 0;

    //! Return a string representation of the concept or role that includes its complexity.
    [[nodiscard]] std::string as_str_with_complexity() const {
        return std::to_string(complexity_) + "." + as_str();
    }

    void force_complexity(int c) {
        complexity_ = c;
    }

    [[nodiscard]] virtual const DLBaseElement* clone() const = 0;

    friend std::ostream& operator<<(std::ostream &os, const DLBaseElement &base) {
        return os << base.as_str_with_complexity() << std::flush;
    }
};

class Concept : public DLBaseElement {
public:
    explicit Concept(int complexity) : DLBaseElement(complexity) { }
    virtual ~Concept() = default;
    [[nodiscard]] const Concept* clone() const override = 0;
};

class Role : public DLBaseElement {
public:
    explicit Role(int complexity) : DLBaseElement(complexity) { }
    virtual ~Role() = default;
    [[nodiscard]] virtual const Role* clone() const override = 0;
};

class PrimitiveConcept : public Concept {
protected:
    const Predicate *predicate_;

public:
    explicit PrimitiveConcept(const Predicate *predicate) : Concept(PRIMITIVE_COMPLEXITY), predicate_(predicate) { }
    ~PrimitiveConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new PrimitiveConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();

        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            res->emplace_back(denotation(cache, sample, sample.state(i)));
        }
        return res;
    }

    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        state_denotation_t sd(sample.num_objects(state.id()), false);
        for( int i = 0; i < int(state.atoms().size()); ++i ) {
            atom_id_t id = state.atoms()[i];
            const Atom &atom = sample.atom(state.id(), id);
            if( atom.is_instance(*predicate_) ) {
                assert(atom.objects().size() == 1);
                object_id_t index = atom.object(0);
                assert(index < sample.num_objects(state.id()));
                sd[index] = true;
            }
        }
        return cache.find_or_insert_state_denotation(sd);
    }
    [[nodiscard]] std::string as_str() const override {
        return predicate_->name_;
    }
};


class NominalConcept : public Concept {
protected:
    const std::string& name_;

public:
    explicit NominalConcept(const std::string &name) : Concept(1), name_(name) {}
    ~NominalConcept() override = default;

    [[nodiscard]] const Concept* clone() const override {
        return new NominalConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();

        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            res->emplace_back(denotation(cache, sample, sample.state(i)));
        }
        return res;
    }


    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        const auto& oidx = sample.instance(state.id()).object_index();
        object_id_t id = oidx.left.at(name_);
        assert(id < sample.num_objects(state.id()));

        state_denotation_t sd(sample.num_objects(state.id()), false);
        sd[id] = true;

        return cache.find_or_insert_state_denotation(sd);
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Nominal(") + name_ + ")";
    }
};

class UniversalConcept : public Concept {
public:
    UniversalConcept() : Concept(0) {}

    ~UniversalConcept() override = default;

    [[nodiscard]] const Concept* clone() const override {
        return new UniversalConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            const auto n = sample.num_objects(i);

            state_denotation_t sd(n, true);
            res->emplace_back(cache.find_or_insert_state_denotation(sd));
        }
        return res;
    }

    [[nodiscard]] std::string as_str() const override {
        return "<universe>";
    }
};


class EmptyConcept : public Concept {
public:
    EmptyConcept() : Concept(0) { }
    ~EmptyConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new EmptyConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();

        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            state_denotation_t std(sample.num_objects(i), false);
            res->emplace_back(cache.find_or_insert_state_denotation(std));
        }
        return res;
    }

    [[nodiscard]] std::string as_str() const override {
        return "<empty>";
    }
};

class AndConcept : public Concept {
protected:
    const Concept *concept1_;
    const Concept *concept2_;

public:
    AndConcept(const Concept *concept1, const Concept *concept2) :
            Concept(1 + concept1->complexity() + concept2->complexity()),
//        Concept(1 + concept1->complexity() * concept2->complexity()),
            concept1_(concept1),
            concept2_(concept2) {
    }
    ~AndConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new AndConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*concept1_, m);
        const sample_denotation_t& sd_sub2 = cache.find_sample_denotation(*concept2_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for( int i = 0; i < m; ++i ) {
            const auto n = sample.num_objects(i);
            const auto& sd1 = sd_sub1.get(i, n);
            const auto& sd2 = sd_sub2.get(i, n);

            state_denotation_t nsd(n, false);
            for (int j = 0; j < n; ++j) nsd[j] = sd1[j] && sd2[j];
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("And(") +  concept1_->as_str() + "," + concept2_->as_str() + ")";
    }
};

class NotConcept : public Concept {
protected:
    const Concept *concept_;

public:
    explicit NotConcept(const Concept *concept)
        : Concept(1 + concept->complexity()),
        concept_(concept)
    {}
    ~NotConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new NotConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*concept_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            const auto n = sample.num_objects(i);
            const auto& std_sub1 = sd_sub1.get(i, n);

            state_denotation_t nsd(n, false);
            for (int j = 0; j < n; ++j) nsd[j] = !std_sub1[j];
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Not(") + concept_->as_str() + ")";
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
    ~ExistsConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new ExistsConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*concept_, m);
        const sample_denotation_t& sd_sub2 = cache.find_sample_denotation(*role_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for( int i = 0; i < m; ++i ) {
            const auto n = sample.num_objects(i);
            const auto& c_den = sd_sub1.get(i, n);
            const auto& r_den = sd_sub2.get(i, n*n);

            state_denotation_t nsd(n, false);
            for (int x = 0; x < n; ++x) {
                // x makes it into the denotation if there is an y such that y in c_den and (x,y) in r_den
                for (unsigned y = 0; y < n; ++y) {
                    if(c_den[y]) {
                        auto x_y = x * m + y;
                        if (r_den[x_y]) {
                            nsd[x] = true;
                            break;
                        }
                    }
                }
            }
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Exists(") + role_->as_str() + "," + concept_->as_str() + ")";
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
    ~ForallConcept() override = default;
    [[nodiscard]] const Concept* clone() const override {
        return new ForallConcept(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*concept_, m);
        const sample_denotation_t& sd_sub2 = cache.find_sample_denotation(*role_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for( int i = 0; i < m; ++i ) {
            const auto n = sample.num_objects(i);
            const auto& c_den = sd_sub1.get(i, n);
            const auto& r_den = sd_sub2.get(i, n*n);

            state_denotation_t nsd(n, true);
            for (int x = 0; x < n; ++x) {
                // x does *not* make it into the denotation if there is an y
                // such that y not in c_den and (x,y) in r_den
                for (unsigned y = 0; y < n; ++y) {
                    if(!c_den[y]) {
                        auto x_y = x * m + y;
                        if (r_den[x_y]) {
                            nsd[x] = false;
                            break;
                        }
                    }
                }
            }
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Forall(") + role_->as_str() + "," + concept_->as_str() + ")";
    }
};


class EqualConcept : public Concept {
protected:
    const Role *r1_;
    const Role *r2_;

public:
    EqualConcept(const Role *r1, const Role *r2)
            : Concept(1 + r1->complexity() + r2->complexity()),
              r1_(r1),
              r2_(r2) {
    }
    ~EqualConcept() override = default;
    [[nodiscard]] const Concept* clone() const override { return new EqualConcept(*this); }


    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*r1_, m);
        const sample_denotation_t& sd_sub2 = cache.find_sample_denotation(*r2_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for( int i = 0; i < m; ++i ) {
            const auto n = sample.num_objects(i);
            const auto& sd1 = sd_sub1.get(i, n*n);
            const auto& sd2 = sd_sub2.get(i, n*n);

            state_denotation_t nsd(n, false);
            for (int x = 0; x < n; ++x) {
                // If the set of y such that (x, y) in sd1 is equal to the set of z such that (x, z) in sd2,
                // then x makes it into the denotation of this concept
                bool in_denotation = true;
                for (int z = 0; z < n; ++z) {
                    auto idx = x * n + z;
                    if (sd1[idx] != sd2[idx]) {
                        in_denotation = false;
                        break;
                    }
                }

                if (in_denotation) {
                    nsd[x] = true;
                }
            }
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Equal(") + r1_->as_str() + "," + r2_->as_str() + ")";
    }
};


class PrimitiveRole : public Role {
protected:
    const Predicate *predicate_;

public:
    explicit PrimitiveRole(const Predicate *predicate) : Role(PRIMITIVE_COMPLEXITY), predicate_(predicate) { }
    ~PrimitiveRole() override = default;
    [[nodiscard]] const Role* clone() const override {
        return new PrimitiveRole(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();

        auto res = new sample_denotation_t();
        res->reserve(m);
        for( int i = 0; i < m; ++i ) {
            res->emplace_back(denotation(cache, sample, sample.state(i)));
        }
        return res;
    }

    const state_denotation_t* denotation(Cache &cache, const Sample &sample, const State &state) const override {
        const unsigned m = sample.num_objects(state.id());
        state_denotation_t sr(m*m, false);
        for (const auto atom_id : state.atoms()) {
            const Atom &atom = sample.atom(state.id(), atom_id);
            if( atom.is_instance(*predicate_) ) {
                assert(atom.objects().size() == 2);
                unsigned index = atom.object(0) * m + atom.object(1);
                assert(index < m * m);
                sr[index] = true;
            }
        }
        return cache.find_or_insert_state_denotation(sr);
    }

    [[nodiscard]] const Predicate* predicate() const { return predicate_; }

    [[nodiscard]] std::string as_str() const override {
        return predicate_->name_;
    }
};

class PlusRole : public Role {
protected:
    const Role *role_;

public:
    explicit PlusRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }
    ~PlusRole() override = default;
    [[nodiscard]] const Role* clone() const override {
        return new PlusRole(*this);
    }

    // apply Johnson's algorithm for transitive closure
    static void transitive_closure(int num_objects, state_denotation_t &sd) {
        // create adjacency lists
        std::vector<std::vector<int> > adj(num_objects);
        for( int i = 0; i < num_objects; ++i ) {
            for( int j = 0; j < num_objects; ++j ) {
                if( sd[i * num_objects + j] )
                    adj[i].emplace_back(j);
            }
        }

        // apply dfs starting from each vertex
        for( int r = 0; r < num_objects; ++r ) {
            std::vector<bool> visited(num_objects, false);
            std::vector<int> q = adj[r];
            while( !q.empty() ) {
                int i = q.back();
                q.pop_back();
                sd[r * num_objects + i] = true;
                if( !visited[i] ) {
                    visited[i] = true;
                    q.insert(q.end(), adj[i].begin(), adj[i].end());
                }
            }
        }
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*role_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            const auto n = sample.num_objects(i);
            const auto& sd1 = sd_sub1.get(i, n*n);

            state_denotation_t nsd(sd1);
            transitive_closure(n, nsd);
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }

    [[nodiscard]] const Role* role() const { return role_; }

    [[nodiscard]] std::string as_str() const override {
        // ATM let us call these Star(X) to get the same output than the Python module
        return std::string("Star(") + role_->as_str() + ")";
    }
};

class StarRole : public Role {
protected:
    const Role *role_;
    const PlusRole *plus_role_;

public:
    explicit StarRole(const Role *role)
            : Role(1 + role->complexity()),
              role_(role),
              plus_role_(new PlusRole(role)) {
    }
    ~StarRole() override {
        delete plus_role_;
    }

    [[nodiscard]] const Role* clone() const override {
        return new StarRole(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        throw std::runtime_error("UNIMPLEMENTED");
    }

    [[nodiscard]] const Role* role() const { return role_; }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Star(") + role_->as_str() + ")";
    }
};

class InverseRole : public Role {
protected:
    const Role *role_;

public:
    explicit InverseRole(const Role *role) : Role(1 + role->complexity()), role_(role) { }

    ~InverseRole() override = default;

    [[nodiscard]] const Role* clone() const override {
        return new InverseRole(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*role_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            const auto n = sample.num_objects(i);
            const auto& sr = sd_sub1.get(i, n*n);

            state_denotation_t nsd(n*n, false);
            for (int j = 0; j < n; ++j) {
                for(int k = 0; k < n; ++k ) {
                    int index = j * n + k;
                    if (sr[index]) {
                        int inv_index = k * n + j;
                        nsd[inv_index] = true;
                    }
                }
            }
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }


    [[nodiscard]] const Role* role() const { return role_; }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Inverse(") + role_->as_str() + ")";
    }
};

// RoleRestriction are only used for distance features
// and thus they are generated when generating such features
class RoleRestriction : public Role {
protected:
    const Role *role_;
    const Concept *restriction_;

public:
    RoleRestriction(const Role *role, const Concept *restriction)
            : Role(1 + role->complexity() + restriction->complexity()),
              role_(role),
              restriction_(restriction) {
    }

    ~RoleRestriction() override = default;

    [[nodiscard]] const Role* clone() const override {
        return new RoleRestriction(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*role_, m);
        const sample_denotation_t& sd_sub2 = cache.find_sample_denotation(*restriction_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            const auto n = sample.num_objects(i);
            const auto& sr = sd_sub1.get(i, n*n);
            const auto& sd = sd_sub2.get(i, n);

            state_denotation_t nsd(sr);
            for (int j = 0; j < n*n; ++j) {
                if( nsd[j] ) {
                    //int src = j / n;
                    int dst = j % n;
                    nsd[j] = sd[dst];
                }
            }
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }

    [[nodiscard]] const Role* role() const { return role_; }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Restrict(") + role_->as_str() + "," + restriction_->as_str() + ")";
    }
};

//! R - R' represents the set of pairs (a,b) such that R(a,b) and not R'(a,b)
//! Makes sense e.g. when R is a goal predicate
class RoleDifference : public Role {
protected:
    const Role *r1_;
    const Role *r2_;

public:
    RoleDifference(const Role *r1, const Role *r2)
            : Role(1 + r1->complexity() + r2->complexity()),
              r1_(r1),
              r2_(r2) {
    }
    ~RoleDifference() override = default;
    [[nodiscard]] const Role* clone() const override {
        return new RoleDifference(*this);
    }

    const sample_denotation_t* denotation(Cache &cache, const Sample &sample) const override {
        const auto m = sample.num_states();
        const sample_denotation_t& sd_sub1 = cache.find_sample_denotation(*r1_, m);
        const sample_denotation_t& sd_sub2 = cache.find_sample_denotation(*r2_, m);

        auto res = new sample_denotation_t();
        res->reserve(m);
        for (int i = 0; i < m; ++i) {
            const auto n = sample.num_objects(i);
            const auto& sd1 = sd_sub1.get(i, n*n);
            const auto& sd2 = sd_sub2.get(i, n*n);

            state_denotation_t nsd(n*n, false);
            for (int x = 0; x < n*n; ++x) {
                if (sd1[x] && !sd2[x]) {
                    nsd[x] = true;
                }
            }
            res->emplace_back(cache.find_or_insert_state_denotation(nsd));
        }
        return res;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("RoleDifference(") + r1_->as_str() + "," + r2_->as_str() + ")";
    }
};

class Feature {
public:
    Feature() = default;
    virtual ~Feature() = default;
    virtual const Feature* clone() const = 0;

    [[nodiscard]] virtual int complexity() const = 0;
    [[nodiscard]] virtual int value(const Cache &cache, const Sample &sample, const State &state) const = 0;
    [[nodiscard]] virtual std::string as_str() const = 0;

    [[nodiscard]] std::string as_str_with_complexity() const {
        return std::to_string(complexity()) + "." + as_str();
    }

    friend std::ostream& operator<<(std::ostream &os, const Feature &f) {
        return os << f.as_str() << std::flush;
    }

    [[nodiscard]] virtual bool is_boolean() const = 0;
};

class NullaryAtomFeature : public Feature {
protected:
    const Predicate* predicate_;

public:
    explicit NullaryAtomFeature(const Predicate* predicate) : predicate_(predicate) { }
    ~NullaryAtomFeature() override = default;

    [[nodiscard]] const Feature* clone() const override {
        return new NullaryAtomFeature(*this);
    }

    [[nodiscard]] int complexity() const override { // Nullary atoms have complexity 0 by definition
        return 1;
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        // Return true (1) iff some atom in the state has the coindiding predicate ID with the predicate of the
        // nullary atom, since there can only be one single atom with this predicate
        // TODO: This is unnecessarily expensive, and it is not cached anywhere
        for (const auto& atom_id:state.atoms()) {
            const Atom &atom = sample.atom(state.id(), atom_id);
            if( atom.is_instance(*predicate_)) {
                return 1;
            }
        }
        return 0;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Atom[") + predicate_->name_ + "]";
    }

    [[nodiscard]] bool is_boolean() const override { return true; }
};

class BooleanFeature : public Feature {
protected:
    const Concept *concept_;

public:
    explicit BooleanFeature(const Concept *concept) : Feature(), concept_(concept) { }

    ~BooleanFeature() override = default;

    [[nodiscard]] const Feature* clone() const override {
        return new BooleanFeature(concept_);
    }

    [[nodiscard]] int complexity() const override {
        return concept_->complexity();
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        // we retrieve the sample denotation from the cache, then the state denotation from the sample denotation,
        // and compute the cardinality (this assumes that state id is index of state into sample.states())
        const auto m = sample.num_states();
        const sample_denotation_t& d = cache.find_sample_denotation(*concept_, m);
        const state_denotation_t& std = d.get(state.id(), sample.num_objects(state.id()));
        assert(std.cardinality() < 2);
        return std.cardinality();
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Bool[") + concept_->as_str() + "]";
    }

    [[nodiscard]] bool is_boolean() const override { return true; }
};

class NumericalFeature : public Feature {
protected:
    const Concept *concept_;

public:
    explicit NumericalFeature(const Concept *concept) : Feature(), concept_(concept) { }

    ~NumericalFeature() override = default;

    [[nodiscard]] const Feature* clone() const override {
        return new NumericalFeature(concept_);
    }

    [[nodiscard]] int complexity() const override {
        return concept_->complexity();
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        // we retrieve the sample denotation from the cache, then the state denotation from the sample denotation,
        // and compute the cardinality (this assumes that state id is index of state into sample.states())
        const auto m = sample.num_states();
        const sample_denotation_t& d = cache.find_sample_denotation(*concept_, m);
        const state_denotation_t& std = d.get(state.id(), sample.num_objects(state.id()));
        return std.cardinality();
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("Num[") + concept_->as_str() + "]";
    }

    [[nodiscard]] bool is_boolean() const override { return false; }
};

class DistanceFeature : public Feature {
protected:
    const Concept *start_;
    const Concept *end_;
    const Role *role_;

    mutable bool valid_cache_;
    mutable std::vector<int> cached_distances_;

public:
    DistanceFeature(const Concept *start, const Concept *end, const Role *role)
            : Feature(),
              start_(start),
              end_(end),
              role_(role),
              valid_cache_(false)
    {}
    ~DistanceFeature() override = default;
    const Feature* clone() const override {
        auto *f = new DistanceFeature(start_, end_, role_);
        f->valid_cache_ = valid_cache_;
        f->cached_distances_ = cached_distances_;
        return f;
    }

    int complexity() const override {
        return 1 + start_->complexity() + end_->complexity() + role_->complexity();
    }
    int value(const Cache &cache, const Sample &sample, const State &state) const override {
        assert(sample.state(state.id()).id() == state.id());
        const auto k = sample.num_states();
        if( !valid_cache_ ) {
            const sample_denotation_t& start_d = cache.find_sample_denotation(*start_, k);
            const sample_denotation_t& end_d = cache.find_sample_denotation(*end_, k);
            const sample_denotation_t& role_d = cache.find_sample_denotation(*role_, k);

            cached_distances_ = std::vector<int>(k, std::numeric_limits<int>::max());
            for(int i = 0; i < k; ++i ) {
                const auto m = sample.num_objects(i);

                const state_denotation_t& start_sd = start_d.get(i, m);
                const state_denotation_t& end_sd = end_d.get(i, m);
                const state_denotation_t& role_sd = role_d.get(i, m*m);
                int distance = compute_distance(m, start_sd, end_sd, role_sd);
                cached_distances_[i] = distance;
            }
            valid_cache_ = true;
        }
        return cached_distances_[state.id()];
    }
    [[nodiscard]] std::string as_str() const override {
        return std::string("Dist[") + start_->as_str() + ";" + role_->as_str() + ";" + end_->as_str() + "]";
    }

    bool is_boolean() const override { return false; }
};


class ConditionalFeature : public Feature {
protected:
    const Feature* condition_;
    const Feature* body_;

public:
    ConditionalFeature(const Feature* condition, const Feature* body) :
            Feature(), condition_(condition), body_(body) { }

    ~ConditionalFeature() override = default;
    [[nodiscard]] const Feature* clone() const override {
        return new ConditionalFeature(condition_, body_);
    }

    [[nodiscard]] int complexity() const override {
        return 1 + condition_->complexity() + body_->complexity();
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        return value_from_components(cache, sample, state, condition_, body_);
    }

    static int value_from_components(const Cache &cache, const Sample &sample, const State &state,
                                     const Feature* condition, const Feature* body) {
        const auto condval = condition->value(cache, sample, state);
        return (condval) ? body->value(cache, sample, state) : std::numeric_limits<int>::max();
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("If{") + condition_->as_str() + "}{" + body_->as_str() + "}{Infty}";
    }

    [[nodiscard]] bool is_boolean() const override { return false; }
};


//! A feature with value f1 < f2
class DifferenceFeature : public Feature {
protected:
    const Feature* f1;
    const Feature* f2;

public:
    DifferenceFeature(const Feature* f1, const Feature* f2) :
            Feature(), f1(f1), f2(f2)
    {}

    ~DifferenceFeature() override = default;
    [[nodiscard]] const Feature* clone() const override {
        return new DifferenceFeature(f1, f2);
    }

    [[nodiscard]] int complexity() const override {
        return 1 + f1->complexity() + f2->complexity();
    }

    [[nodiscard]] int value(const Cache &cache, const Sample &sample, const State &state) const override {
        const auto f1val = f1->value(cache, sample, state);
        const auto f2val = f2->value(cache, sample, state);
        return f1val < f2val;
    }

    [[nodiscard]] std::string as_str() const override {
        return std::string("LessThan{") + f1->as_str() + "}{" + f2->as_str() + "}";
    }

    [[nodiscard]] bool is_boolean() const override { return true; }
};

//! We use this to store a number of properties of the denotations of concepts
//! and features over entire samples; properties which we might be interested in analyzing
//! for diverse ends such as pruning redundant features, etc.
struct SampleDenotationProperties {
    bool denotation_is_bool = false;
    bool denotation_is_constant = false;
};


class Factory {
protected:
    const std::vector<std::string> nominals_;
    std::vector<const Role*> basis_roles_;
    std::vector<const Concept*> basis_concepts_;
    
    Options options;

    mutable std::vector<const Role*> roles_;

    // A layered set of concepts, concepts_[k] contains all concepts generated in the k-th application of
    // the concept grammar
    mutable std::vector<std::vector<const Concept*> > concepts_;
    std::vector<const Feature*> features_;

    //! Indices of features generated from goal concepts
    std::unordered_set<unsigned> goal_features_;

public:
    Factory(std::vector<std::string>  nominals, Options options) :
            nominals_(std::move(nominals)),
            options(std::move(options))
    {}
    virtual ~Factory() = default;

    void insert_basis(const Role *role) {
        basis_roles_.push_back(role);
    }

    void insert_basis(const Concept *concept) {
//        std::cout << "insert_basis: " << *concept << std::endl;
        basis_concepts_.push_back(concept);
    }


    //! Insert the given concept/role as long as it is not redundant with some previous role and it is
    //! below the complexity bound. Return whether the role was effectively inserted.
    template <typename T1, typename T2>
    bool attempt_insertion(const T1& elem, Cache &cache, const Sample &sample, std::vector<const T2*>& container) const {
        if (elem.complexity() > options.complexity_bound) {
//            std::cout << elem.as_str() << " superfluous because complexity " << base.complexity() << ">" << options.complexity_bound << std::endl;
            return false;
        }

        const sample_denotation_t *d = elem.denotation(cache, sample);

        const auto& index = cache.cache1();
        auto it = index.find(d);
        if (it != index.end()) {
            // There is in the index some other concept/role with same sample denotation,
            // hence we consider this one redundant
//              std::cout << elem.as_str() << " superfluous because complexity "
//                        << elem.complexity() << ">" << options.complexity_bound << std::endl;
            delete d;
            return false;

        } else {

            container.push_back(elem.clone());
            assert (container.back()->id() == elem.id());
            cache.find_or_insert_sample_denotation(*d, elem.id());
            delete d;
            return true;
        }
    }

    bool check_timeout(const std::clock_t& start_time) const {
        if (options.timeout <= 0) return false;
        double elapsed = (std::clock() - start_time) / (double) CLOCKS_PER_SEC;
        if (elapsed > options.timeout) {
            std::cout << "\tTimeout of " << options.timeout << " sec. reached while generating concepts" << std::endl;
            return true;
        }
        return false;
    }

    // apply one iteration of the concept generation grammar
    // new concepts are left on a new last layer in concepts_
    // new concepts are non-redundant if sample != nullptr
    int advance_step(Cache &cache, const Sample &sample, const std::clock_t& start_time) const {
        int num_pruned_concepts = 0;
        if( concepts_.empty() ) { // On the first iteration, we simply process the basis concepts and return
            concepts_.emplace_back();
            for (const auto* concept : basis_concepts_) {
                num_pruned_concepts += !attempt_insertion(*concept, cache, sample, concepts_.back());
            }
            return num_pruned_concepts;
        }


        // Classify concepts and roles by their complexity - we will use this to minimize complexity checks later
        std::vector<std::vector<const Role*>> roles_by_complexity(options.complexity_bound+1);
        std::vector<std::vector<const Concept*>> concepts_in_last_layer_by_complexity(options.complexity_bound+1);
        std::vector<std::vector<const Concept*>> concepts_in_previous_layers_by_complexity(options.complexity_bound+1);

        for (const auto* r:roles_) roles_by_complexity.at(r->complexity()).push_back(r);  // TODO Do this just once
        for (const auto* c:concepts_.back()) concepts_in_last_layer_by_complexity.at(c->complexity()).push_back(c);

        for( unsigned layer = 0; layer < concepts_.size()-1; ++layer ) {
            for (const auto* c:concepts_[layer]) {
                concepts_in_previous_layers_by_complexity.at(c->complexity()).push_back(c);
            }
        }

        // create new concept layer
        bool is_first_non_basis_iteration = (concepts_.size() == 1);
        concepts_.emplace_back();

        if (is_first_non_basis_iteration) {
            // Let's create equal concepts for those pairs of roles (whose number is already fixed at this point)
            // such that both arise from same predicate and one of them is the "goal version" of the other, e.g.
            // on and on_g, in blocksworld.
            for (const auto* r1:roles_) {
                for (const auto* r2:roles_) {
                    if (dynamic_cast<const RoleDifference*>(r1) || dynamic_cast<const RoleDifference*>(r2)) {
                        // R=R' Makes no sense when R or R' are already a role difference
                        continue;
                    }
                    auto name1 = get_role_predicate(r1)->name();
                    auto name2 = get_role_predicate(r2)->name();

                    // A hacky way to allow only equal concepts R=R' where R and R' come from the same predicate
                    // and the first one is a goal role
                    if (name1.substr(name1.size()-2) != "_g") continue;
                    if (name1.substr(0, name1.size()-2) != name2) continue;

                    EqualConcept eq_concept(r1, r2);
                    // Force the complexity of R=R' to be 1, if both roles are of the same type
                    // This is currently disabled, as the concept-pruning strategy doesn't work as expected when
                    // this is enabled - sometimes concepts with lower complexity get pruned in favor of others
                    // with higher complexity
                    // if (typeid(*r1) == typeid(*r2)) eq_concept.force_complexity(1);

                    num_pruned_concepts += !attempt_insertion(eq_concept, cache, sample, concepts_.back());

                    if (check_timeout(start_time)) return -1;
                }
            }
        }

        for (int k = 0; k <= (options.complexity_bound-1); ++k) {
            for (const auto* concept:concepts_in_last_layer_by_complexity[k]) {

                // Negate concepts in the last layer, only if they are not already negations
                if (!dynamic_cast<const NotConcept*>(concept)) {
                    num_pruned_concepts += !attempt_insertion(NotConcept(concept), cache, sample, concepts_.back());
                }


                // generate exist and forall combining a role with a concept in the last layer
                for (int k2 = 0; k2 <= (options.complexity_bound-k-1); ++k2) {
                    for (const auto *role:roles_by_complexity[k2]) {
                        num_pruned_concepts += !attempt_insertion(ExistsConcept(concept, role), cache, sample, concepts_.back());
                        num_pruned_concepts += !attempt_insertion(ForallConcept(concept, role), cache, sample, concepts_.back());

                        if (check_timeout(start_time)) return -1;
                    }
                }
                if (check_timeout(start_time)) return -1;
            }
        }


        // generate conjunctions of (a) a concept of the last layer with a concept in some previous layer, and
        // (b) two concepts of the last layer, avoiding symmetries
        for (int k = 0; k <= (options.complexity_bound-1); ++k) {
            for (unsigned i_k = 0; i_k < concepts_in_last_layer_by_complexity[k].size(); ++i_k) {
                const auto& concept1 = concepts_in_last_layer_by_complexity[k][i_k];

                // (a) a concept of the last layer with a concept in some previous layer
                for (int k2 = 0; k2 <= (options.complexity_bound-k); ++k2) {
                    for (const auto *concept2:concepts_in_previous_layers_by_complexity[k2]) {
                        num_pruned_concepts += !attempt_insertion(AndConcept(concept1, concept2), cache, sample, concepts_.back());

                        if (check_timeout(start_time)) return -1;
                    }
                }


                // (b) two concepts of the last layer, avoiding symmetries
                for (int k2 = k; k2 <= (options.complexity_bound-k); ++k2) {
                    unsigned start_at = (k == k2) ? i_k+1: 0; // Break symmetries within same complexity bucket
                    for (unsigned i_k2 = start_at; i_k2 < concepts_in_last_layer_by_complexity[k2].size(); ++i_k2) {
                        const auto& concept2 = concepts_in_last_layer_by_complexity[k2][i_k2];
                        num_pruned_concepts += !attempt_insertion(AndConcept(concept1, concept2), cache, sample, concepts_.back());

                        if (check_timeout(start_time)) return -1;
                    }
                }
            }
        }
        return num_pruned_concepts;
    }

    //! Retrieve the predicate at the basis of a given role (given the current grammar restrictions, there will be
    //! exactly one such predicate
    static const Predicate* get_role_predicate(const Role* r) {

        if (const auto *c = dynamic_cast<const PrimitiveRole*>(r)) {
            return c->predicate();

        } else if (const auto* c = dynamic_cast<const StarRole*>(r)) {
            return get_role_predicate(c->role());

        } else if (const auto *c = dynamic_cast<const PlusRole*>(r)) {
            return get_role_predicate(c->role());

        } else if (const auto *c = dynamic_cast<const RoleRestriction*>(r)) {
            return get_role_predicate(c->role());

        } else if (const auto *c = dynamic_cast<const InverseRole*>(r)) {
            return get_role_predicate(c->role());

        } else if (const auto *c = dynamic_cast<const RoleDifference*>(r)) {
            throw std::runtime_error("Unimplemented");
        }
        throw std::runtime_error("Unknown role type");
    }

    void generate_basis(const Sample &sample) {
        insert_basis(new UniversalConcept);
        insert_basis(new EmptyConcept);
        for( int i = 0; i < int(sample.num_predicates()); ++i ) {
            const Predicate *predicate = &sample.predicates()[i];
            if( predicate->arity() == 1 ) {
                insert_basis(new PrimitiveConcept(predicate));
            } else if( predicate->arity() == 2 ) {
                insert_basis(new PrimitiveRole(predicate));
            }
        }

        for (const auto& nominal:nominals_) {
            insert_basis(new NominalConcept(nominal));
        }
        std::cout << "BASIS: #concepts=" << basis_concepts_.size() << ", #roles=" << basis_roles_.size() << std::endl;
//        report_dl_data(std::cout);
    }

    std::vector<const Concept*> generate_goal_concepts_and_roles(Cache &cache, const Sample &sample) {
        std::vector<const Concept*> goal_concepts;

        const Sample::PredicateIndex& predicate_idx = sample.predicate_index();
        for (const predicate_id_t& pid:sample.goal_predicates()) {
            const Predicate& pred = sample.predicate(pid);
            unsigned arity = pred.arity();
            const auto& name = pred.name();
            auto name_n = name.size();
            assert(name[name_n-2] == '_' && name[name_n-1] == 'g');  // Could try something more robust, but ATM this is OK
            auto nongoal_name = name.substr(0, name_n-2);
            const Predicate& nongoal_pred = sample.predicate(predicate_idx.at(nongoal_name));
            if (arity == 0) {
                // No need to do anything here, as the corresponding NullaryAtomFeature will always be generated
                // with complexity 1
                assert (0); // This will need some refactoring
//                goal_features.push_back(new NullaryAtomFeature(&pred));

            } else if (arity == 1) {
                // Force into the basis a concept "p_g AND Not p"
                auto *c = new PrimitiveConcept(&nongoal_pred);
                auto *not_c = new NotConcept(c);
                auto *c_g = new PrimitiveConcept(&pred);
                auto *and_c = new AndConcept(c_g, not_c);

                // Insert the denotations into the cache
//                for (const auto elem:std::vector<const Concept*>({and_c})) {
                for (const auto* elem:std::vector<const Concept*>({c, not_c, c_g, and_c})) {
                    const sample_denotation_t *d = elem->denotation(cache, sample);
                    cache.find_or_insert_sample_denotation(*d, elem->id());
                }
//                for (const auto elem:std::vector<const Concept*>({c, not_c, c_g})) {
//                    cache.remove_sample_denotation(elem->as_str());
//                }


                and_c->force_complexity(1);
                for (const auto* elem:std::vector<const Concept*>({c, not_c, c_g, and_c})) {
                    insert_basis(elem);
                }

//                goal_features.push_back(new NumericalFeature(and_c));
                goal_concepts.push_back(and_c);

            } else if (arity == 2) {
                // Force into the basis a new RoleDifference(r_g, r)
                auto *r = new PrimitiveRole(&nongoal_pred);
                auto *r_g = new PrimitiveRole(&pred);
                auto *r_diff = new RoleDifference(r_g, r);

                // Insert the denotations into the cache
//                for (const auto elem:std::vector<const Role*>({r_diff})) {
//                for (const auto elem:std::vector<const Role*>({r, r_g, r_diff})) {
//                    const sample_denotation_t *d = elem->denotation(cache, sample, false);
//                    cache.find_or_insert_sample_denotation(*d, elem->as_str());
//                }
//                for (const auto elem:std::vector<const Role*>({r, r_g})) {
//                    cache.remove_sample_denotation(elem->as_str());
//                }

                r_diff->force_complexity(1);
                for (const auto* elem:std::vector<const Role*>({r, r_g, r_diff})) {
                    insert_basis(elem);
                }

                auto *ex_c = new ExistsConcept(new UniversalConcept, r_diff); // memleak
//                goal_features.push_back(new NumericalFeature(ex_c));
                goal_concepts.push_back(ex_c);
            }
        }
        return goal_concepts;
    }

    int generate_roles(Cache &cache, const Sample &sample) const {
        assert(roles_.empty());

        std::vector<const Role*> non_redundant_base_roles;

        // Insert the basis (i.e. primitive) roles as long as they are not redundant
        for (const auto *role : basis_roles_) {
            if (attempt_insertion(*role, cache, sample, roles_)) {
                non_redundant_base_roles.push_back(role);
            }
        }

        // Now insert a few compounds based on those base roles that are not redundant
        for (const auto *role:non_redundant_base_roles) {

            // Create Inverse(R) role from the primitive role
            attempt_insertion(InverseRole(role), cache, sample, roles_);

            // Create Plus(R) role from the primitive role
            PlusRole p_role(role);
            if (attempt_insertion(p_role, cache, sample, roles_)) {
                // Create Inverse(Plus(R)) only if Plus(R) is NOT redundant
                attempt_insertion(InverseRole(p_role.clone()), cache, sample, roles_);
            }
            // Create Star(R) roles from the primitive roles !!! NOTE ATM we deactivate Star roles
            // attempt_role_insertion(StarRole(role), cache, sample);
        }
        std::cout << "ROLES: #roles=" << roles_.size() << std::endl;
        return roles_.size();
    }

    std::vector<const Concept*> generate_concepts(Cache &cache, const Sample &sample, const std::clock_t& start_time) const {
        std::size_t num_concepts = 0;
        bool some_new_concepts = true;
        bool timeout_reached = false;
        for( int iteration = 0; some_new_concepts && !timeout_reached; ++iteration ) {
            std::cout << "dl::concept-generation: iteration=" << iteration
                      << ", #concepts=" << num_concepts
                      << ", #concepts-in-last-layer=" << (concepts_.empty() ? 0 : concepts_.back().size())
                      << std::endl;

            std::size_t num_concepts_before_step = num_concepts;

            // Let the fun begin:
            auto num_pruned_concepts = advance_step(cache, sample, start_time);
            timeout_reached = num_pruned_concepts < 0;

            num_concepts += concepts_.empty() ? 0 : concepts_.back().size();
            some_new_concepts = num_concepts > num_concepts_before_step;

            std::cout << "\tResult of iteration:"
                      << " #concepts-in-layer=" << concepts_.back().size()
                      << ", #pruned-concepts=" << num_pruned_concepts
                      << std::endl;
//            report_dl_data(std::cout);
        }

        if (concepts_.back().empty()) {  // i.e. we stopped because last iteration was already a fixpoint
            concepts_.pop_back();
        }

        // Important to use this vector from here on, as it is sorted by complexity
        auto all_concepts = consolidate_concepts();
        assert(all_concepts.size() == num_concepts);

        std::cout << "dl::Factory: #concepts-final=" << num_concepts << std::endl;
        return all_concepts;
    }

    void generate_comparison_features(
            const std::vector<const Feature*>& base_features,
            Cache& cache,
            const Sample& sample,
            const TransitionSample& transitions,
            feature_cache_t& seen_denotations)
    {
        if (!options.comparison_features) return;

        // get the max index here, as we'll keep adding elements to the same `features_` vector:
        auto n = features_.size();

        for (std::size_t i = 0; i < n; ++i) {
            const auto* f_i = features_[i];
            if (f_i->is_boolean() || f_i->complexity() + 1 + 1 > options.complexity_bound) continue;

            for (std::size_t j = 0; j < n; ++j) {
                const auto* f_j = features_[j];
                if (i == j || f_j->is_boolean() || f_i->complexity() + f_j->complexity() + 1 > options.complexity_bound)
                    continue;

                const auto *feature = new DifferenceFeature(f_i, f_j);
                if (!attempt_feature_insertion(
                        feature, options.complexity_bound, cache, sample, transitions, seen_denotations, true)) {
                    delete feature;
                }
            }
        }
    }

    void generate_conditional_features(
            const std::vector<const Feature*>& base_features,
            Cache& cache,
            const Sample& sample,
            const TransitionSample& transitions,
            feature_cache_t& seen_denotations)
    {
        if (options.cond_complexity_bound <= 0) return;

        for (std::size_t i = 0, n = features_.size(); i < n; ++i) {
            const auto* cond = features_[i];
            if (!cond->is_boolean() || cond->complexity() + 1 + 1 > options.cond_complexity_bound) continue;

            for (std::size_t j = i + 1; j < n; ++j) {
                const auto* body = features_[j];
                if (body->is_boolean() ||
                    cond->complexity() + body->complexity() + 1 > options.cond_complexity_bound)
                    continue;

                const auto *feature = new ConditionalFeature(cond, body);
                if (!attempt_feature_insertion(
                        feature, options.cond_complexity_bound, cache, sample, transitions, seen_denotations, true)) {
                    delete feature;
                }
            }
        }
    }

    void generate_features(
            const std::vector<const Concept*>& concepts,
            Cache &cache, const Sample &sample,
            const TransitionSample& transitions,
            const std::vector<const Concept*>& forced_goal_features);

    void print_feature_count() const {
        unsigned num_nullary_features = 0, num_boolean_features = 0, num_numeric_features = 0,
                num_distance_features = 0, num_conditional_features = 0, num_comparison_features = 0;
        auto nf = features_.size();
        for (const auto *f:features_) {
            if (dynamic_cast<const NullaryAtomFeature*>(f)) num_nullary_features++;
            else if (dynamic_cast<const BooleanFeature*>(f)) num_boolean_features++;
            else if (dynamic_cast<const NumericalFeature*>(f)) num_numeric_features++;
            else if (dynamic_cast<const DistanceFeature*>(f)) num_distance_features++;
            else if (dynamic_cast<const ConditionalFeature*>(f)) num_conditional_features++;
            else if (dynamic_cast<const DifferenceFeature*>(f)) num_comparison_features++;
            else throw std::runtime_error("Unknown feature type");
        }

        assert(nf == num_nullary_features+num_boolean_features+num_numeric_features
                     +num_distance_features+num_conditional_features+num_comparison_features);
        std::cout << "FEATURES: #features=" << nf
                  << ", #nullary="   << num_nullary_features
                  << ", #boolean="   << num_boolean_features
                  << ", #numerical=" << num_numeric_features
                  << ", #distance="  << num_distance_features
                  << ", #conditional="  << num_conditional_features
                  << ", #comparison="  << num_comparison_features
                  << std::endl;
    }

    inline static feature_sample_denotation_t compute_feature_sample_denotation(
            const Feature& feature, const Sample &sample, const Cache &cache) {
        SampleDenotationProperties _;
        return compute_feature_sample_denotation(feature, sample, cache, _);
    }

    static feature_sample_denotation_t compute_feature_sample_denotation(
            const Feature& feature, const Sample &sample, const Cache &cache, SampleDenotationProperties& properties);

    bool attempt_cardinality_feature_insertion(
            const Concept* c,
            Cache &cache,
            const Sample &sample,
            const TransitionSample& transitions,
            feature_cache_t& seen_denotations,
            bool check_redundancy);

    //! Insert the given feature if its complexity is below the given bound, its denotation is not constant,
    //! and its denotation trail is not redundant with that of some previously-generated feature
    //! Return whether the feature was indeed inserted or not
    bool attempt_feature_insertion(
            const Feature* feature, unsigned bound,
            Cache &cache, const Sample &sample,
            const TransitionSample& transitions, feature_cache_t& seen_denotations,
            bool check_redundancy);

    void generate_distance_features(
            const std::vector<const Concept*>& concepts, Cache &cache,
            const Sample &sample,
            const TransitionSample& transitions,
            feature_cache_t& seen_denotations) {
        if (options.dist_complexity_bound<=0) return;

        const auto m = sample.num_states();

        // Identify concepts with singleton denotation across all states: these are the candidates for start concepts
        std::vector<const Concept*> start_concepts;
        for (const Concept* c:concepts) {
            const sample_denotation_t& d = cache.find_sample_denotation(*c, m);
            bool singleton_denotations = true;
            for (int j = 0; singleton_denotations && (j < m); ++j) {
                if (d[j]->cardinality() != 1) {
                    singleton_denotations = false;
                    break;
                }
            }

            if (singleton_denotations) {
                start_concepts.push_back(c);
            }
        }

        // create role restrictions to be used in distance features
        std::vector<const Role*> role_restrictions(roles_);  // Start with all base roles
        for (const Role* r:roles_) {
            for (const Concept* c:concepts) {
                RoleRestriction role_restriction(r, c);
                if( role_restriction.complexity()+3 > options.dist_complexity_bound ) continue;
                const sample_denotation_t *d = role_restriction.denotation(cache, sample);

                if (!cache.contains(d)) {  // The role is not redundant
                    role_restrictions.push_back(role_restriction.clone());
                    cache.find_or_insert_sample_denotation(*d, role_restriction.id());
                    //std::cout << "ACCEPT RR(sz=" << cache_for_role_restrictions.cache1().size() << "): "
                    // + role_restriction.as_str_with_complexity() << std::endl;
                } else {
                    //std::cout << "PRUNE RR: " + role_restriction.as_str() << std::endl;
                }
                delete d;
            }
        }

        // create distance features
        int num_distance_features = 0;
        std::vector<const DistanceFeature*> candidates;

        for (const Concept* start:start_concepts) {
            for (const Concept* end:concepts) {
                if (start == end) continue;

                for (const Role* role:role_restrictions) {
                    const auto* df = new DistanceFeature(start, end, role);
                    if (df->complexity() > options.dist_complexity_bound) {
                        delete df;
                        continue;
                    }

                    candidates.push_back(df);

                    // fill cache with denotations for start and end concepts
//            const sample_denotation_t& ds = cache.find_sample_denotation(*start, m);
//            cache.find_or_insert_sample_denotation(ds, start->id());
//            const sample_denotation_t& de = cache.find_sample_denotation(*end, m);
//            cache.find_or_insert_sample_denotation(de, end->id());
                }
            }
        }

        // Sort the candidate distance features along increasing complexity
        std::sort(std::begin(candidates), std::end(candidates),
                  [](const DistanceFeature* f1, const DistanceFeature* f2) {
            return f1->complexity() < f2->complexity();
        });

        for (const auto* df:candidates) {
            SampleDenotationProperties properties;
            const auto denotation = compute_feature_sample_denotation(*df, sample, cache, properties);

            if (!prune_feature_denotation(
                    *df, denotation, properties, sample, transitions, seen_denotations, true)) {
                ++num_distance_features;
                features_.emplace_back(df);
                seen_denotations.emplace(denotation, features_.back());
            } else {
                delete df;
            }
        }
    }

    static bool prune_feature_denotation(
            const Feature& f,
            const feature_sample_denotation_t& fd,
            const SampleDenotationProperties& properties,
            const Sample &sample,
            const TransitionSample& transitions,
            feature_cache_t& seen_denotations,
            bool check_redundancy);

    std::ostream& report_dl_data(std::ostream &os) const {

        unsigned nconcepts = 0;
        for (const auto& layer:concepts_) nconcepts += layer.size();
        auto nroles = roles_.size();
        os << "Total number of concepts / roles: " << nconcepts << "/" << nroles << std::endl;

        os << "Base concepts (sz=" << basis_concepts_.size() << "): ";
        for( int i = 0; i < int(basis_concepts_.size()); ++i ) {
            os << basis_concepts_[i]->as_str();
            if( 1 + i < int(basis_concepts_.size()) ) os << ", ";
        }
        os << std::endl;

        os << "Base roles (sz=" << basis_roles_.size() << "): ";
        for( int i = 0; i < int(basis_roles_.size()); ++i ) {
            os << basis_roles_[i]->as_str();
            if( 1 + i < int(basis_roles_.size()) ) os << ", ";
        }
        os << std::endl;

        os << "All Roles under complexity " << options.complexity_bound << " (sz=" << roles_.size() << "): ";
        for( int i = 0; i < int(roles_.size()); ++i ) {
            os << roles_[i]->as_str_with_complexity();
            if( 1 + i < int(roles_.size()) ) os << ", ";
        }
        os << std::endl;

        os << "All concepts (by layer) under complexity " << options.complexity_bound << ": " << std::endl;
        for( int layer = 0; layer < concepts_.size(); ++layer ) {
            os << "    Layer " << layer << " (sz=" << concepts_[layer].size() << "): ";
            for( int i = 0; i < int(concepts_[layer].size()); ++i ) {
                os << concepts_[layer][i]->as_str_with_complexity();
                if( 1 + i < int(concepts_[layer].size()) ) os << ", ";
            }
            os << std::endl;
        }

        os << "Nullary-atom features: ";
        bool need_comma = false;
        for (const auto &feature : features_) {
            if( dynamic_cast<const NullaryAtomFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Boolean features: ";
        need_comma = false;
        for (const auto &feature : features_) {
            if( dynamic_cast<const BooleanFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Numerical features: ";
        need_comma = false;
        for (const auto &feature : features_) {
            if( dynamic_cast<const NumericalFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        os << "Distance features: ";
        need_comma = false;
        for (const auto &feature : features_) {
            if( dynamic_cast<const DistanceFeature*>(feature) ) {
                if( need_comma ) os << ", ";
                os << feature->as_str();
                need_comma = true;
            }
        }
        os << std::endl;

        return os;
    }

    void output_feature_matrix(std::ostream &os, const Cache &cache, const Sample &sample) const {
        auto num_features = features_.size();
        for( int i = 0; i < sample.num_states(); ++i ) {
            const State &state = sample.state(i);

            // one line per state with the numeric denotation of all features
            for( int j = 0; j < num_features; ++j ) {
                const Feature &feature = *features_[j];
                os << feature.value(cache, sample, state);
                if( 1 + j < num_features ) os << " ";
            }
            os << std::endl;
        }
    }

    void output_feature_info(std::ostream &os, const Cache &cache, const Sample &sample) const {
        auto num_features = features_.size();

        // Line #1: feature names
        for( int i = 0; i < num_features; ++i ) {
            const Feature &feature = *features_[i];
            os << feature.as_str();
            if( 1 + i < num_features ) os << "\t";
        }
        os << std::endl;

        // Line #2: feature complexities
        for( int i = 0; i < num_features; ++i ) {
            const Feature &feature = *features_[i];
            os << feature.complexity();
            if( 1 + i < num_features ) os << "\t";
        }
        os << std::endl;

        // Line #3: feature types (0: boolean; 1: numeric)
        for( int i = 0; i < num_features; ++i ) {
            const Feature* feature = features_[i];
            os << (feature->is_boolean() ? 0 : 1);
            if( 1 + i < num_features ) os << "\t";
        }
        os << std::endl;

        // Line #4: whether feature is goal feature (0: No; 1: yes)
        for( int i = 0; i < num_features; ++i ) {
            auto it = goal_features_.find(i);
            os << (it == goal_features_.end() ? 0 : 1);
            if( 1 + i < num_features ) os << "\t";
        }
        os << std::endl;
    }

    void log_all_concepts_and_features(const std::vector<const Concept*>& concepts,
                                       const sltp::dl::Cache &cache, const Sample &sample,
                                       const std::string& workspace, bool print_denotations);

    //! Return all generated concepts in a single, *unlayered* vector, and sorted by increasing complexity
    std::vector<const Concept*> consolidate_concepts() const {
        std::vector<const Concept*> all;
        for (const auto& layer:concepts_) all.insert(all.end(), layer.begin(), layer.end());

        std::sort(std::begin(all), std::end(all), [](const Concept* c1, const Concept* c2) {
            return c1->complexity() < c2->complexity();
        });

        return all;
    }

    static bool check_some_transition_pair_distinguished(
            const feature_sample_denotation_t &fsd, const Sample &sample, const TransitionSample &transitions) ;
};

} // namespaces
