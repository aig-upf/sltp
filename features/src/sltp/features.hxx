

#pragma once

#include <utility>
#include <cstdint>
#include <vector>
#include <stdexcept>

namespace sltp {

    class DLDenotation {
    public:
        enum class denotation_t : bool {
            concept_denotation_t,
            role_denotation_t
        };

        std::vector<bool> _values;
        denotation_t _type;
    };

    //! A DL denotation matrix is just a matrix of DL Denotations, each state being a row, each concept / role a column
    class DLDenotationMatrix {

        std::vector<std::vector<DLDenotation>> data;
    };

    class DLModel {

    };


    class DLBase {
    protected:
//        unsigned _id;

        unsigned _complexity;

    public:
        explicit DLBase(unsigned complexity) : _complexity(complexity) {}

//        const unsigned id() const { return _id; }
        const unsigned complexity() const { return _complexity; }

        virtual DLDenotation denotation() const = 0;

        //! Print a representation of the object to the given stream.
        friend std::ostream& operator<<(std::ostream &os, const DLBase& o) { return o.print(os); }
        virtual std::ostream& print(std::ostream& os) const = 0;
    };

    class Concept : public DLBase {
    public:
        explicit Concept(unsigned complexity) : DLBase(complexity) {}

    };

    class Role : public DLBase {
    public:
        explicit Role(unsigned complexity) : DLBase(complexity) {}

    };

    class PrimitiveConcept : public Concept {
    protected:
        const std::string _name;

    public:
        PrimitiveConcept(std::string name) : Concept(1), _name(std::move(name)) {}

        DLDenotation denotation() const override {
            throw std::runtime_error("TODO: UNIMPLEMENTED");
            return DLDenotation();
        }

        std::ostream& print(std::ostream& os) const override {
            return os << _name;
        }

    };

    class PrimitiveRole : public Role {
    protected:
        const std::string _name;

    public:
        PrimitiveRole(std::string name) : Role(1), _name(std::move(name)) {}

        DLDenotation denotation() const override {
            throw std::runtime_error("TODO: UNIMPLEMENTED");
            return DLDenotation();
        }

        std::ostream& print(std::ostream& os) const override {
            return os << _name;
        }

    };

    class ExistsConcept : public Concept {
    protected:
        const Concept* _c;
        const Role* _r;

    public:
        ExistsConcept(const Concept* c, const Role* r) : Concept(c->complexity() + r->complexity() + 1), _c(c), _r(r) {}

        DLDenotation denotation() const override {
            throw std::runtime_error("TODO: UNIMPLEMENTED");
            return DLDenotation();
        }
    };




}