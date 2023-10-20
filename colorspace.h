#pragma once
#include "color.h"
#include "spectrum.h"

class RGBColorSpace
{
public : 
    RGBColorSpace(vec2 _r, vec2 _g, vec2 _b, const Spectrum& _illuminant) :
        r(_r), g(_g), b(_b), illuminant(_illuminant)
    {
        XYZ W = SpectrumToXYZ(illuminant);
        w = W.xy();
        XYZ R = XYZ::FromxyY(r);
        XYZ G = XYZ::FromxyY(g);
        XYZ B = XYZ::FromxyY(b);

        SquareMatrix<3> rgb(R.x, G.x, B.x, R.y, G.y, B.y, R.z, G.z, B.z);
        XYZ C = XYZ(rgb.Invert() * vec3(W.x, W.y, W.z));
        XYZFromRGB = rgb * SquareMatrix<3>(C.x, 0, 0, 0, C.y, 0, 0, 0, C.z);
        RGBFromXYZ = XYZFromRGB.Invert();
    }

    RGBColor ToRGB(XYZ xyz) const {
        vec3 rgb = RGBFromXYZ * vec3(xyz.x, xyz.y, xyz.z);
        return RGBColor(rgb.x(), rgb.y(), rgb.z());
    }

    XYZ ToXYZ(RGBColor rgb) const {
        vec3 xyz = XYZFromRGB * vec3(rgb.r, rgb.g, rgb.b);
        return XYZ(xyz.x(), xyz.y(), xyz.z());
    }

    RGBSigmoidPolynomial ToRGBCoeffs(RGBColor rgb) const
    {
        return sRGBTable(rgb);
    }

private : 
    vec2 r, g, b, w;
    SquareMatrix<3> XYZFromRGB, RGBFromXYZ;
    DenselySampledSpectrum illuminant;
};

const RGBColorSpace sRGB(vec2(.64, .33), vec2(.3, .6), vec2(.15, .06), SpectrasRGB);
