#pragma once
#include "math.h"
#include "declaration.h"

class RGBColor
{
public : 
    double r, g, b;
    RGBColor() : r(0), g(0), b(0) {}
    RGBColor(double _r, double _g, double _b) : r(_r), g(_g), b(_b) {}
    RGBColor(vec3 x) : r(x.x()), g(x.y()), b(x.z()) {}
};

class RGBSigmoidPolynomial
{
public:
    RGBSigmoidPolynomial(double _c0, double _c1, double _c2) : c0(_c0), c1(_c1), c2(_c2) {}

    static double s(double x) {
        if (std::isinf(x)) return x > 0 ? 1 : 0;
        return 0.5 + x / (2.0 * std::sqrt(1 + x * x));
    }

    double operator()(double lambda) const {
        double x = c0 * lambda * lambda + c1 * lambda + c2;
        return s(x);
    }

    double MaxValue() const {
        double result = fmax((*this)(360), (*this)(830));
        double lambda = -c1 / (2 * c0);
        if (lambda >= 360 && lambda <= 830) result = fmax(result, (*this)(lambda));
        return result;
    }

private:
    double c0, c1, c2;
};

class RGBToSpectrumTable
{
public : 
    static constexpr int res = 64;
    using CoefficientArray = float[3][res][res][res][3];

    RGBToSpectrumTable(const float* _zNodes, const CoefficientArray* _coeffs) : zNodes(_zNodes), coeffs(_coeffs) {}

    RGBSigmoidPolynomial operator ()(RGBColor rgb) const;

private : 
    const float* zNodes;
    const CoefficientArray* coeffs;
};

extern const float sRGBToSpectrumTable_Scale[64];
extern const float sRGBToSpectrumTable_Data[3][64][64][64][3];
static const RGBToSpectrumTable sRGBTable(sRGBToSpectrumTable_Scale, &sRGBToSpectrumTable_Data);

