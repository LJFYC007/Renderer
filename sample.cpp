#include "sample.h"

thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
thread_local std::mt19937 generator(std::random_device{}());