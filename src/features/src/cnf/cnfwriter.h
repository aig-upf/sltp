#pragma once

#include <utility>
#include <cstdint>
#include <vector>
#include <ostream>

using cnfvar_t = uint32_t;
using cnflit_t = int64_t;  // To make sure we can negate any variable
using cnfclause_t = std::vector<cnflit_t>;


class CNFWriter {
protected:
    uint32_t next_var_id_;

    ulong accumulated_weight_;

    ulong nclauses_;

    std::ostream& os_;

public:
    // Variable IDs must start with 1
    explicit CNFWriter(std::ostream &os) : next_var_id_(1), accumulated_weight_(0), nclauses_(0), os_(os) {}

    cnfvar_t variable() {
        next_var_id_ += 1;
        return next_var_id_-1;
    }

    uint32_t nvars() const { return next_var_id_ - 1; }

    ulong nclauses() const { return nclauses_; }

    static inline cnflit_t literal(cnfvar_t var, bool polarity) {
        auto res = static_cast<cnflit_t>(var);
        return polarity ? res : -1*res;
    }

    ulong top() const { return accumulated_weight_ + 1; }

    //! Print the given clause - if weight is negative, it is assumed to be TOP
    void print_clause(const cnfclause_t& clause, int weight = -1) {
        //  w <literals> 0
        assert(!clause.empty());
        if (weight == 0) throw std::runtime_error("Cannot use weight-0 clauses");
        auto size = clause.size();
        if (weight < 0) os_ << "TOP ";
        else {
            accumulated_weight_ += weight;
            os_ << weight << " ";
        }
        for (unsigned i = 0; i < size; ++i) {
            os_ << clause[i] << " ";
        }
        os_ << "0" << std::endl;

        nclauses_ += 1;
    }
};