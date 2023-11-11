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
inline vec2 vec2Random() { return vec2(randomDouble(), randomDouble()); }
inline vec3 vec3Random() { return vec3(randomDouble(), randomDouble(), randomDouble()); }
inline vec3 vec3Random(double min, double max) { return vec3(randomDouble(min, max), randomDouble(min, max), randomDouble(min, max)); }

static vec3 randInDisk()
{
    while (true)
    {
        vec3 p = vec3(randomDouble(-1, 1), randomDouble(-1, 1), 0);
        if (p.length() < 1.0) return p;
    }
}

static vec3 SampleUniformHemisphere(const vec2& u) {
    double z = u[0];
    double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    double phi = 2 * pi * u[1];
    return vec3(r * std::cos(phi), r * std::sin(phi), z);
}

inline double UniformHemispherePDF() { return 1 / (2 * pi); }

inline bool SameHemisphere(vec3 x, vec3 y) { return x.z() * y.z() > 0.0; }

static vec2 SampleConcentricDisk(vec2 u)
{
    vec2 _u = u * 2 - vec2(1);
    if (_u.x() == 0 && _u.y() == 0) return vec2(0);
    double theta, r;
    if (std::abs(_u.x()) > std::abs(_u.y()))
    {
        r = _u.x();
        theta = pi / 4 * (_u.y() / _u.x());
    }
    else
    {
        r = _u.y();
        theta = pi / 2 - pi / 4 * (_u.x() / _u.y());
    }
    return r * vec2(std::cos(theta), std::sin(theta));
}

static vec3 SampleCosineHemisphere(vec2 u)
{
    vec2 d = SampleConcentricDisk(u);
    double z = std::sqrt(std::max(0.0, 1 - d.x() * d.x() - d.y() * d.y()));
    return vec3(d.x(), d.y(), z);
}

inline double CosineHemispherePDF(double cosTheta) { return cosTheta / pi; }

