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
	virtual std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const = 0;
};

class GeometricPrimitive : public Primitive
{
public:
	GeometricPrimitive(shared_ptr<Shape> _shape, shared_ptr<Material> _material, shared_ptr<Light> _areaLight) : shape(_shape), material(_material), areaLight(_areaLight) {}
	AABB Bounds() const override { return shape->Bounds(); }

	std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const override {
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

	std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const override {
		std::optional<ShapeIntersection> isect = shape->Intersect(r, t);
		if (!isect) return {};
		isect->intr.SetIntersectionProperties(material, nullptr);
		return isect;
	}

private:
	shared_ptr<Shape> shape;
	shared_ptr<Material> material;
};

class PrimitiveList : public Primitive
{
public : 
	vector<shared_ptr<Primitive>> objects;
	AABB box;

	PrimitiveList() {}
	PrimitiveList(shared_ptr<Primitive> object) { add(object); }

	void clear() { objects.clear(); }
	void add(shared_ptr<Primitive> object) { objects.push_back(object); box = AABB(box, object->Bounds()); }
	AABB Bounds() const override { return box; }
	std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const override { assert(-1); return {}; }
};

class BvhNode : public Primitive
{
public:
	shared_ptr<Primitive> left, right;
	AABB box;

	BvhNode(shared_ptr<PrimitiveList> list) : BvhNode(list -> objects, 0, list -> objects.size()) {}
	BvhNode(const std::vector<shared_ptr<Primitive>>& oldObjects, size_t start, size_t end) {
		auto objects = oldObjects;
		int axis = randomInt(0, 2);
		auto comparator = (axis == 0) ? box_x_compare : ((axis == 1) ? box_y_compare : box_z_compare);
		size_t n = end - start;

		if (n == 1) left = right = objects[start];
		else if (n == 2) {
			if (comparator(objects[start], objects[start + 1]))
				left = objects[start], right = objects[start + 1];
			else
				left = objects[start + 1], right = objects[start];
		}
		else
		{
			std::sort(objects.begin() + start, objects.begin() + end, comparator);
			auto mid = start + n / 2;
			left = make_shared<BvhNode>(objects, start, mid);
			right = make_shared<BvhNode>(objects, mid, end);
		}
		box = AABB(left->Bounds(), right->Bounds());
	}

	std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const override
	{
		double t0, t1;
		if (!box.Intersect(r, t, t0, t1)) return {};
		std::optional<ShapeIntersection> hitLeft = left->Intersect(r, t);
		std::optional<ShapeIntersection> hitRight = right->Intersect(r, interval(t.Min, hitLeft ? hitLeft->tHit : t.Max));
		return hitRight ? hitRight : hitLeft;
	}

	AABB Bounds() const override { return box; }

private:
	static bool box_compare(const shared_ptr<Primitive> a, const shared_ptr<Primitive> b, int axis_index) {
		return a->Bounds().axis(axis_index).Min < b->Bounds().axis(axis_index).Min;
	}

	static bool box_x_compare(const shared_ptr<Primitive> a, const shared_ptr<Primitive> b) {
		return box_compare(a, b, 0);
	}

	static bool box_y_compare(const shared_ptr<Primitive> a, const shared_ptr<Primitive> b) {
		return box_compare(a, b, 1);
	}

	static bool box_z_compare(const shared_ptr<Primitive> a, const shared_ptr<Primitive> b) {
		return box_compare(a, b, 2);
	}
};
