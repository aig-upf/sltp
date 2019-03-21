
#pragma once

#include <string>
#include <fstream>
#include <blai/utils.h>


std::ofstream open_file(const std::string &filename) {
    std::ofstream ofs(filename.c_str());
    if(ofs.fail()) {
        throw std::runtime_error(Utils::error() + "opening file '" + filename + "'");
    }
    return ofs;
}