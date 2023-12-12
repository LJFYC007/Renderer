#pragma once
#include "math.h"

vec3 operator *(const SquareMatrix<3>& m, const vec3& v)
{
	return vec3(m[0][0] * v.a[0] + m[0][1] * v.a[1] + m[0][2] * v.a[2],
		m[1][0] * v.a[0] + m[1][1] * v.a[1] + m[1][2] * v.a[2],
		m[2][0] * v.a[0] + m[2][1] * v.a[1] + m[2][2] * v.a[2]
	);
}

vec3 Transform::operator()(const vec3& p) const
{
    double x = p[0], y = p[1], z = p[2];
    double X = mat[0][0] * x + mat[0][1] * y + mat[0][2] * z + mat[0][3];
    double Y = mat[1][0] * x + mat[1][1] * y + mat[1][2] * z + mat[1][3];
    double Z = mat[2][0] * x + mat[2][1] * y + mat[2][2] * z + mat[2][3];
    double W = mat[3][0] * x + mat[3][1] * y + mat[3][2] * z + mat[3][3];

    assert(W != 0);
    if (W == 1) return vec3(X, Y, Z);
    return vec3(X, Y, Z) / W;
}

Vector3fi Transform::operator()(const Vector3fi& p) const
{
    double x = p.x.Midpoint(), y = p.y.Midpoint(), z = p.z.Midpoint();
    double X = mat[0][0] * x + mat[0][1] * y + mat[0][2] * z + mat[0][3];
    double Y = mat[1][0] * x + mat[1][1] * y + mat[1][2] * z + mat[1][3];
    double Z = mat[2][0] * x + mat[2][1] * y + mat[2][2] * z + mat[2][3];
    double W = mat[3][0] * x + mat[3][1] * y + mat[3][2] * z + mat[3][3];
    assert(W != 0);

    vec3 pError;
    if (p.IsExact()) {
        pError.a[0] = gamma(3) * (std::abs(mat[0][0] * x) + std::abs(mat[0][1] * y) + std::abs(mat[0][2] * z) + mat[0][3]);
        pError.a[1] = gamma(3) * (std::abs(mat[1][0] * x) + std::abs(mat[1][1] * y) + std::abs(mat[1][2] * z) + mat[1][3]);
        pError.a[2] = gamma(3) * (std::abs(mat[2][0] * x) + std::abs(mat[2][1] * y) + std::abs(mat[2][2] * z) + mat[2][3]);
    }
    else {
        vec3 pInError = p.Error();
        pError.a[0] = (gamma(3) + 1) * (std::abs(mat[0][0]) * pInError.a[0] + std::abs(mat[0][1]) * pInError.a[1] + std::abs(mat[0][2]) * pInError.a[2]) +
            gamma(3) * (std::abs(mat[0][0] * x) + std::abs(mat[0][1] * y) + std::abs(mat[0][2] * z) + std::abs(mat[0][3]));
        pError.a[1] = (gamma(3) + 1) * (std::abs(mat[1][0]) * pInError.a[0] + std::abs(mat[1][1]) * pInError.a[1] + std::abs(mat[1][2]) * pInError.a[2]) +
            gamma(3) * (std::abs(mat[1][0] * x) + std::abs(mat[1][1] * y) + std::abs(mat[1][2] * z) + std::abs(mat[1][3]));
        pError.a[2] = (gamma(3) + 1) * (std::abs(mat[2][0]) * pInError.a[0] + std::abs(mat[2][1]) * pInError.a[1] + std::abs(mat[2][2]) * pInError.a[2]) +
            gamma(3) * (std::abs(mat[2][0] * x) + std::abs(mat[2][1] * y) + std::abs(mat[2][2] * z) + std::abs(mat[2][3]));
    }

    if (W == 1) return Vector3fi(vec3(X, Y, Z), pError);
    return Vector3fi(vec3(X, Y, Z), pError) / W;
}

Transform Transform::Translate(const vec3& delta)
{
    double matData[4][4] = {
        {1, 0, 0, delta.x()},
        {0, 1, 0, delta.y()},
        {0, 0, 1, delta.z()},
        {0, 0, 0, 1}
    };
    return Transform(matData);
}

Transform Transform::Scale(double x, double y, double z)
{
    double matData[4][4] = {
        {x, 0, 0, 0},
        {0, y, 0, 0},
        {0, 0, z, 0},
        {0, 0, 0, 1}
    };
    return Transform(matData);
}

Transform Transform::RotateX(double theta)
{
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);

    double matData[4][4] = {
        {1, 0, 0, 0},
        {0, cosTheta, -sinTheta, 0},
        {0, sinTheta, cosTheta, 0},
        {0, 0, 0, 1}
    };
    return Transform(matData);
}

Transform Transform::RotateY(double theta)
{
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);

    double matData[4][4] = {
        {cosTheta, 0, sinTheta, 0},
        {0, 1, 0, 0},
        {-sinTheta, 0, cosTheta, 0},
        {0, 0, 0, 1}
    };
    return Transform(matData);
}

Transform Transform::RotateZ(double theta)
{
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);

    double matData[4][4] = {
        {cosTheta, -sinTheta, 0, 0},
        {sinTheta, cosTheta, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    return Transform(matData);
}

Transform Transform::Rotate(double theta, const vec3& axis)
{
    vec3 a = normalize(axis);
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);
    double oneMinusCosTheta = 1.0 - cosTheta;
    double x = a.x(), y = a.y(), z = a.z();
    double x2 = x * x;
    double y2 = y * y;
    double z2 = z * z;
    double xy = x * y;
    double xz = x * z;
    double yz = y * z;

    double matData[4][4] = {
        {x2 + (1 - x2) * cosTheta, xy * oneMinusCosTheta - z * sinTheta, xz * oneMinusCosTheta + y * sinTheta, 0},
        {xy * oneMinusCosTheta + z * sinTheta, y2 + (1 - y2) * cosTheta, yz * oneMinusCosTheta - x * sinTheta, 0},
        {xz * oneMinusCosTheta - y * sinTheta, yz * oneMinusCosTheta + x * sinTheta, z2 + (1 - z2) * cosTheta, 0},
        {0, 0, 0, 1}
    };
    return Transform(matData);
}
