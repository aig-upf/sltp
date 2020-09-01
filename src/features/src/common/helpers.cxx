

#include <common/helpers.h>
#include <blai/utils.h>
#include "utils.h"

Sample::TransitionSample read_transition_data(const std::string& workspace, bool verbose) {
    std::string transitions_filename = workspace + "/sat_transitions.dat";
    std::cout << Utils::blue() << "reading" << Utils::normal() << " '" << transitions_filename << std::endl;
    auto ifs_transitions = get_ifstream(transitions_filename);
    auto transitions = Sample::TransitionSample::read_dump(ifs_transitions, verbose);
    ifs_transitions.close();
    return transitions;
}


Sample::FeatureMatrix read_feature_matrix(const std::string& workspace, bool verbose) {
    std::string matrix_filename = workspace + "/feature-matrix.dat";
    std::cout << Utils::blue() << "reading" << Utils::normal() << " '" << matrix_filename << std::endl;
    auto ifs_matrix = get_ifstream(matrix_filename);
    auto matrix = Sample::FeatureMatrix::read_dump(ifs_matrix, verbose);
    ifs_matrix.close();
    return matrix;
}
