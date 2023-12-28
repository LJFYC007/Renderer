#pragma once
#include "texture.h"
#include "color.h"
#include "colorspace.h"
#include "image.h"

#include <cmath>
#include <string>
#include <iostream>
#include <tiny_gltf.h>

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
	SpectrumConstantTexture(vec3 _color) : color(_color), value(make_shared<RGBAlbedoSpectrum>(sRGB, RGBColor(color))) {}
	SpectrumConstantTexture(shared_ptr<Spectrum> _value) : value(_value) {}
	vec3 Evaluate(TextureEvalContext ctx) const override { return color; }
	SampledSpectrum Evaluate(TextureEvalContext ctx, SampledWaveLengths lambda) const override { return value->Sample(lambda); }
private:
	vec3 color;
	shared_ptr<Spectrum> value;
};

struct TexCoord2D {
	vec2 st;
	double dsdx, dsdy, dtdx, dtdy;
};

class UVMapping {
public:
	UVMapping(double _su = 1, double _sv = 1, double _du = 0, double _dv = 0) : su(_su), sv(_sv), du(_du), dv(_dv) {}
	TexCoord2D Map(int texCoord, TextureEvalContext ctx) const {
		double dsdx = su * ctx.dudx, dsdy = su * ctx.dudy;
		double dtdx = sv * ctx.dvdx, dtdy = sv * ctx.dvdy;
		vec2 st;
		if (texCoord == 0) st = vec2(su * ctx.uv[0] + du, sv * ctx.uv[1] + dv);
		else st = vec2(su * ctx.UV[0] + du, sv * ctx.UV[1] + dv);

		st[0] = std::fmod(st[0], 1.0);
		st[1] = std::fmod(st[1], 1.0);
		if (st[0] < 0.0) st[0] += 1.0;
		if (st[1] < 0.0) st[1] += 1.0;

		return TexCoord2D{ st, dsdx, dsdy, dtdx, dtdy };
	}

private:
	double su, sv, du, dv;
};

class AlphaImageTexture : public DoubleTexture 
{
public:
	AlphaImageTexture(int _texCoord, UVMapping _mapping, shared_ptr<Image> _image) :
		texCoord(_texCoord), mapping(_mapping), image(_image) {}

	double Evaluate(TextureEvalContext ctx) const override {
		TexCoord2D c = mapping.Map(texCoord, ctx);
		return image->LookUpAlpha(c.st);
	}

private:
	int texCoord;
	UVMapping mapping;
	shared_ptr<Image> image;
};

class ImageTexture : public SpectrumTexture
{
public:
	ImageTexture(int _texCoord, UVMapping _mapping, shared_ptr<Image> _image, bool _gammaCorrection) :
		texCoord(_texCoord), mapping(_mapping), image(_image), gammaCorrection(_gammaCorrection) {}

	vec3 Evaluate(TextureEvalContext ctx) const override {
		TexCoord2D c = mapping.Map(texCoord, ctx);
		return image->LookUp(c.st, gammaCorrection);
	}

	SampledSpectrum Evaluate(TextureEvalContext ctx, SampledWaveLengths lambda) const override {
		TexCoord2D c = mapping.Map(texCoord, ctx);
		return image->LookUp(c.st, lambda, gammaCorrection);
	}

private:
	int texCoord;
	UVMapping mapping;
	shared_ptr<Image> image;
	bool gammaCorrection;
};
