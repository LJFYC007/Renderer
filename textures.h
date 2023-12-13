#pragma once
#include "texture.h"

class DoubleConstantTexture : public DoubleTexture
{
public:
	DoubleConstantTexture(double _value) : value(_value) {}
	double Evaluate(TextureEvalContext ctx) const override { return value; }
private:
	double value;
};

class SpectrumConstantTexture : public SpectrumTexture
{
public:
	SpectrumConstantTexture(vec3 color) : value(make_shared<RGBAlbedoSpectrum>(sRGB, RGBColor(color))) {}
	SpectrumConstantTexture(shared_ptr<Spectrum> _value) : value(_value) {}
	SampledSpectrum Evaluate(TextureEvalContext ctx, SampledWaveLengths lambda) const override { return value->Sample(lambda); }
private:
	shared_ptr<Spectrum> value;
};
