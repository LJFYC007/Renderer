#pragma once
#include "texture.h"
#include "color.h"
#include "colorspace.h"
#include "stb_image.h"

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

class ImageTexture : public SpectrumTexture
{
public:
	ImageTexture(int _texCoord, UVMapping _mapping, shared_ptr<tinygltf::Image> _image) :
		texCoord(_texCoord), mapping(_mapping), image(_image) {}

	vec3 Evaluate(TextureEvalContext ctx) const override {
		TexCoord2D c = mapping.Map(texCoord, ctx);
		int i = static_cast<int>(c.st[0] * (image->width - 1));
		int j = static_cast<int>(c.st[1] * (image->height - 1));
		i = std::max(0, std::min(i, image->width - 1));
		j = std::max(0, std::min(j, image->height - 1));
		int pixelIndex = (j * image->width + i) * image->component;
		double r = image->image.data()[pixelIndex] / 255.0f;
		double g = image->image.data()[pixelIndex + 1] / 255.0f;
		double b = image->image.data()[pixelIndex + 2] / 255.0f;
		return vec3(r, g, b);
	}

	SampledSpectrum Evaluate(TextureEvalContext ctx, SampledWaveLengths lambda) const override {
		vec3 color = Evaluate(ctx);
		return RGBAlbedoSpectrum(sRGB, RGBColor(color)).Sample(lambda);
	}

private:
	int texCoord;
	UVMapping mapping;
	shared_ptr<tinygltf::Image> image;
};
