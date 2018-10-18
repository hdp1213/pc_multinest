#ifndef PROFILE_H
#define PROFILE_H

#include "Parameter.hpp"

#include <string>
#include <utility>

void profile_likelihood(std::string& root, param_t profiled_param,
                        std::pair<double, double>& range, double fixed_values[],
                        int num_samples, void* context);

#endif
