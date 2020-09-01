
#pragma once

#include <blai/transitions.h>
#include <blai/matrix.h>

Sample::TransitionSample read_transition_data(const std::string& workspace, bool verbose);

Sample::FeatureMatrix read_feature_matrix(const std::string& workspace, bool verbose);
