#include "color.h"

#include <algorithm>

template <typename T, typename U, typename V>
inline constexpr T Clamp(T val, U low, V high) {
    if (val < low)
        return T(low);
    else if (val > high)
        return T(high);
    else
        return val;
}

template <typename Predicate>
inline size_t FindInterval(size_t sz, const Predicate& pred) {
    using ssize_t = std::make_signed_t<size_t>;
    ssize_t size = (ssize_t)sz - 2, first = 1;
    while (size > 0) {
        // Evaluate predicate at midpoint and update _first_ and _size_
        size_t half = (size_t)size >> 1, middle = first + half;
        bool predResult = pred(middle);
        first = predResult ? middle + 1 : first;
        size = predResult ? size - (half + 1) : half;
    }
    return (size_t)Clamp((ssize_t)first - 1, 0, sz - 2);
}

RGBSigmoidPolynomial RGBToSpectrumTable::operator()(RGBColor _rgb) const
{
    vec3 rgb(_rgb.r, _rgb.g, _rgb.b);
    // Handle uniform _rgb_ values
    if (rgb[0] == rgb[1] && rgb[1] == rgb[2])
        return RGBSigmoidPolynomial(0, 0,
            (rgb[0] - .5f) / std::sqrt(rgb[0] * (1 - rgb[0])));

    // Find maximum component and compute remapped component values
    int maxc =
        (rgb[0] > rgb[1]) ? ((rgb[0] > rgb[2]) ? 0 : 2) : ((rgb[1] > rgb[2]) ? 1 : 2);
    double z = rgb[maxc];
    double x = rgb[(maxc + 1) % 3] * (res - 1) / z;
    double y = rgb[(maxc + 2) % 3] * (res - 1) / z;

    // Compute integer indices and offsets for coefficient interpolation
    int xi = std::min((int)x, res - 2), yi = std::min((int)y, res - 2),
        zi = FindInterval(res, [&](int i) { return zNodes[i] < z; });
    double dx = x - xi, dy = y - yi, dz = (z - zNodes[zi]) / (zNodes[zi + 1] - zNodes[zi]);

    // Trilinearly interpolate sigmoid polynomial coefficients _c_
    double c[3];
    for (int i = 0; i < 3; ++i) {
        // Define _co_ lambda for looking up sigmoid polynomial coefficients
        auto co = [&](int dx, int dy, int dz) {
            return (*coeffs)[maxc][zi + dz][yi + dy][xi + dx][i];
            };

        c[i] = Lerp(dz,
            Lerp(dy, Lerp(dx, co(0, 0, 0), co(1, 0, 0)),
                Lerp(dx, co(0, 1, 0), co(1, 1, 0))),
            Lerp(dy, Lerp(dx, co(0, 0, 1), co(1, 0, 1)),
                Lerp(dx, co(0, 1, 1), co(1, 1, 1))));
    }

    return RGBSigmoidPolynomial(c[0], c[1], c[2]);
}
