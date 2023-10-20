#pragma once
#include "spectrum.h"
#include "colorspace.h"

XYZ SampledSpectrum::ToXYZ(const SampledWaveLengths& lambda) const
{
    SampledSpectrum X = SpectraX.Sample(lambda);
    SampledSpectrum Y = SpectraY.Sample(lambda);
    SampledSpectrum Z = SpectraZ.Sample(lambda);
    SampledSpectrum pdf = lambda.PDF();
    return XYZ(((X * *this) / pdf).Average(),
        ((Y * *this) / pdf).Average(),
        ((Z * *this) / pdf).Average()) / CIE_Y_INTEGRAL;
}

RGBColor SampledSpectrum::ToRGB(const SampledWaveLengths& lambda, const RGBColorSpace& space) const
{
    XYZ xyz = ToXYZ(lambda);
    return space.ToRGB(xyz);
}

PiecewiseLinearSpectrum::PiecewiseLinearSpectrum(const double* samples, int n, bool normalize)
{
	if (samples[0] >= LambdaMin)
	{
		lambdas.push_back(LambdaMin - 1);
		values.push_back(samples[1]);
	}
	for (int i = 0; i < n; ++i)
	{
		if (i > 0) assert(samples[2 * i] >= lambdas.back());
		lambdas.push_back(samples[2 * i]);
		values.push_back(samples[2 * i + 1]);
	}
	if (lambdas.back() <= LambdaMax)
	{
		lambdas.push_back(LambdaMax + 1);
		values.push_back(values.back());
	}

	if (normalize)
	{
		double x = CIE_Y_INTEGRAL / InnerProduct(*this, DenselySampledSpectrum(CIE_Y));
		for (auto& it : values) it *= x;
	}
}

RGBAlbedoSpectrum::RGBAlbedoSpectrum(const RGBColorSpace& cs, RGBColor rgb) : rsp(cs.ToRGBCoeffs(rgb)) {}
