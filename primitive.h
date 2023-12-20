#pragma once
#include "math.h"
#include "ray.h"
#include "interval.h"
#include "shape.h"
#include "lights.h"
#include "material.h"

#include <memory>
#include <vector>
#include <assert.h>
#include <optional>
using std::shared_ptr;
using std::make_shared;
using std::vector;

class Primitive
{
public: 
	virtual AABB Bounds() const = 0;
	virtual std::optional<ShapeIntersection> Intersect(const Ray& r, interval t) const = 0;
};

class GeometricPrimitive : public Primitive
{
public:
	GeometricPrimitive(shared_ptr<Shape> _shape, shared_ptr<Material> _material, shared_ptr<Light> _areaLight) : shape(_shape), material(_material), areaLight(_areaLight) {}
	AABB Bounds() const override { return shape->Bounds(); }

	std::optional<ShapeIntersection> Intersect(const Ray& r, interval t) const override {
		std::optional<ShapeIntersection> isect = shape->Intersect(r, t);
		if (!isect) return {};
		isect->intr.SetIntersectionProperties(material, areaLight);
		return isect;
	}

private:
	shared_ptr<Shape> shape;
	shared_ptr<Material> material;
	shared_ptr<Light> areaLight;
};

class SimplePrimitive : public Primitive
{
public:
	SimplePrimitive(shared_ptr<Shape> _shape, shared_ptr<Material> _material) : shape(_shape), material(_material) {}
	AABB Bounds() const override { return shape->Bounds(); }

	std::optional<ShapeIntersection> Intersect(const Ray& r, interval t) const override {
		std::optional<ShapeIntersection> isect = shape->Intersect(r, t);
		if (!isect) return {};
		isect->intr.SetIntersectionProperties(material, nullptr);
		return isect;
	}

private:
	shared_ptr<Shape> shape;
	shared_ptr<Material> material;
};

