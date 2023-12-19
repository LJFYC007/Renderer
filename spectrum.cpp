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
		lambdas.emplace_back(LambdaMin - 1);
		values.emplace_back(samples[1]);
	}
	for (int i = 0; i < n; ++i)
	{
		if (i > 0) assert(samples[2 * i] >= lambdas.back());
		lambdas.emplace_back(samples[2 * i]);
		values.emplace_back(samples[2 * i + 1]);
	}
	if (lambdas.back() <= LambdaMax)
	{
		lambdas.emplace_back(LambdaMax + 1);
		values.emplace_back(values.back());
	}

	if (normalize)
	{
		double x = CIE_Y_INTEGRAL / InnerProduct(*this, DenselySampledSpectrum(CIE_Y));
		for (auto& it : values) it *= x;
	}
}

RGBAlbedoSpectrum::RGBAlbedoSpectrum(const RGBColorSpace& cs, const RGBColor& rgb) : rsp(cs.ToRGBCoeffs(rgb)) {}

RGBIlluminantSpectrum::RGBIlluminantSpectrum(const RGBColorSpace& cs, RGBColor rgb)
{
	double m = std::max(std::max(rgb.r, rgb.g), rgb.b);
	scale = 2 * m;
	rsp = cs.ToRGBCoeffs(scale != 0 ? rgb / scale : RGBColor(0.0, 0.0, 0.0));
}
