
#pragma once

#include <string>
#include "features.hxx"

namespace sltp {

    class io {
    public:
        static DLDenotationMatrix read_denotation_matrix(const std::string &filename);
        static void write_denotation_matrix(const DLDenotationMatrix& matrix, const std::string& filename);
    };

}