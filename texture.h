#pragma once
#include "math.h"
#include "spectrum.h"
#include "colorspace.h"
#include "color.h"
#include "interaction.h"

#include <memory>
using std::shared_ptr;
using std::make_shared;

class texture
{
public:
	virtual ~texture() = default;
	virtual RGBAlbedoSpectrum value(double u, double v, const vec3& p) const = 0;
};

class solidColor : public texture
{
private:
	RGBAlbedoSpectrum color;

public:
	solidColor(vec3 _color) : color(sRGB, RGBColor(_color)) {}
	solidColor(double red, double green, double blue) : color(sRGB, RGBColor(vec3(red, green, blue))) {}
	RGBAlbedoSpectrum value() const { return color; }
	RGBAlbedoSpectrum value(double u, double v, const vec3& p) const override { return color; }
};

class checkerBoard : public texture
{
private:
	double scale;
	shared_ptr<texture> even, odd;
public:
	checkerBoard(double _scale, shared_ptr<texture> _even, shared_ptr<texture> _odd) : scale(_scale), even(_even), odd(_odd) {}
	checkerBoard(double _scale, vec3 _even, vec3 _odd) : scale(_scale), even(make_shared<solidColor>(_even)), odd(make_shared<solidColor>(_odd)) {}

	RGBAlbedoSpectrum value(double u, double v, const vec3& p) const override
	{
		int x = static_cast<int>(std::floor(scale * p.x()));
		int y = static_cast<int>(std::floor(scale * p.y()));
		int z = static_cast<int>(std::floor(scale * p.z()));
		bool isEven = (x + y + z) % 2 == 0;
		return isEven ? even->value(u, v, p) : odd->value(u, v, p);
	}
};

struct TextureEvalContext 
{
	TextureEvalContext(const Interaction& intr) : p(intr.pi), uv(intr.uv) {}
	vec3 p, dpdx, dpdy, n;
	vec2 uv;
	double dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
	int faceIndex = 0;
};

class FloatTexture
{
public:
	virtual double Evaluate(TextureEvalContext ctx) const = 0;
};

class SpectrumTexture
{
public:
	virtual SampledSpectrum Evaluate(TextureEvalContext ctx) const = 0;
};
