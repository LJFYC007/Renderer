#pragma once
#include "texture.h"
#include "color.h"
#include "colorspace.h"
#include "stb_image.h"

#include <string>
#include <iostream>

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

struct TexCoord2D {
	vec2 st;
	double dsdx, dsdy, dtdx, dtdy;
};

class UVMapping {
public:
	UVMapping(double _su = 1, double _sv = 1, double _du = 0, double _dv = 0) : su(_su), sv(_sv), du(_du), dv(_dv) {}
	TexCoord2D Map(TextureEvalContext ctx) const {
		double dsdx = su * ctx.dudx, dsdy = su * ctx.dudy;
		double dtdx = sv * ctx.dvdx, dtdy = sv * ctx.dvdy;
		vec2 st(su * ctx.uv[0] + du, sv * ctx.uv[1] + dv);
		return TexCoord2D{ st, dsdx, dsdy, dtdx, dtdy };
	}

private:
	double su, sv, du, dv;
};

class ImageTexture : public SpectrumTexture
{
public:
	ImageTexture(UVMapping _mapping, std::string _filename, double _scale) : mapping(_mapping), filename(_filename), scale(_scale) {
		data = stbi_load(filename.c_str(), &width, &height, &nrChannels, 0);
	}

	std::string GetPath() const { return filename; }

	double DoubleEvaluate(TextureEvalContext ctx) const {
		TexCoord2D c = mapping.Map(ctx);
		int i = static_cast<int>(c.st[0] * (width - 1));
		int j = static_cast<int>(c.st[1] * (height - 1));
		i = std::max(0, std::min(i, width - 1));
		j = std::max(0, std::min(j, height - 1));
		int pixelIndex = (j * width + i) * nrChannels;
		double x = data[pixelIndex] / 255.0f;
		return x * scale;
	}

	vec3 Evaluate(TextureEvalContext ctx) const override {
		TexCoord2D c = mapping.Map(ctx);
		int i = static_cast<int>(c.st[0] * (width - 1));
		int j = static_cast<int>(c.st[1] * (height - 1));
		i = std::max(0, std::min(i, width - 1));
		j = std::max(0, std::min(j, height - 1));
		int pixelIndex = (j * width + i) * nrChannels;
		assert(nrChannels >= 3);
		double r = data[pixelIndex] / 255.0f;
		double g = data[pixelIndex + 1] / 255.0f;
		double b = data[pixelIndex + 2] / 255.0f;
		return vec3(r, g, b) * scale;
	}

	SampledSpectrum Evaluate(TextureEvalContext ctx, SampledWaveLengths lambda) const override {
		vec3 color = Evaluate(ctx);
		return RGBAlbedoSpectrum(sRGB, RGBColor(color)).Sample(lambda);
	}

private:
	UVMapping mapping;
	std::string filename;
	double scale;
	int width, height, nrChannels;
	unsigned char* data;
};
