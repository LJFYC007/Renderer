#pragma once
#include <cmath>
#include <random>

#include "math.h"

inline double radians(double x) { return x * pi / 180.0;  }
inline double fract(double x) { return x - floor(x); }
inline double randomDouble()
{
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}
inline double randomDouble(double min, double max) { return min + randomDouble() * (max - min); }
inline int randomInt(int min, int max) { return static_cast<int>(randomDouble(min, max + 1)); }
inline vec3 vec3Random() { return vec3(randomDouble(), randomDouble(), randomDouble()); }
inline vec3 vec3Random(double min, double max) { return vec3(randomDouble(min, max), randomDouble(min, max), randomDouble(min, max)); }

static vec3 randInSphere()
{
    while (true)
    {
        vec3 p = vec3Random(-1.0, 1.0);
        if (p.length() < 1.0) return p;
    }
}

static vec3 randUnitVector() { return normalize(randInSphere()); }

static vec3 randInDisk()
{
    while (true)
    {
        vec3 p = vec3(randomDouble(-1, 1), randomDouble(-1, 1), 0);
        if (p.length() < 1.0) return p;
    }
}

static vec3 randInHemisphere(double seed, vec3 normal)
{
    vec3 p = normalize(randInSphere());
    return dot(p, normal) > 0.0 ? p : -p;
}

static vec3 SampleUniformHemisphere(const vec2& u) {
    double z = u[0];
    double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    double phi = 2 * pi * u[1];
    return vec3(r * std::cos(phi), r * std::sin(phi), z);
}

static double UniformHemispherePDF() { return 1 / (2 * pi); }
