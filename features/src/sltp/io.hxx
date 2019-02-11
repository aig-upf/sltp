#ifndef IO_HXX
#define IO_HXX

#include <string>
#include "features.hxx"

namespace SLTP {

class IO {
  public:
    static DL::DenotationMatrix read_denotation_matrix(const std::string &filename);
    static void write_denotation_matrix(const DL::DenotationMatrix& matrix, const std::string& filename);
};

};

#endif

