#pragma once
#include "math.h"
#include "ray.h"
#include "material.h"
#include "interval.h"
#include "bvh.h"

#include <memory>
#include <vector>
#include <assert.h>
using std::shared_ptr;
using std::make_shared;
using std::vector;

class material;

class hitRecord
{
public : 
	point3 p;
	vec3 normal;
	double t, u, v;
	bool frontFace;
	shared_ptr<material> mat;
};

class primitive
{
public : 
	virtual ~primitive() = default;
	virtual bool hit(const ray& r, interval t, hitRecord& rec) const = 0;
	virtual bvh boundingBox() const = 0;
};

class primitiveList : public primitive
{
public : 
	vector<shared_ptr<primitive>> objects;
	bvh box;

	primitiveList() {}
	primitiveList(shared_ptr<primitive> object) { add(object); }

	void clear() { objects.clear(); }
	void add(shared_ptr<primitive> object) { objects.push_back(object); box = bvh(box, object->boundingBox()); }
	bvh boundingBox() const override { return box; }

	bool hit(const ray& r, interval t, hitRecord& rec) const override
	{
		assert(-1);
        hitRecord x;
        bool hitAnything = false;
        for ( const auto& object : objects)
            if ( object -> hit(r, t, x) == true )
            {
                hitAnything = true;
				t.Max = x.t;
                rec = x;
            }
        return hitAnything;
	}
};

class bvhNode : public primitive
{
public : 
	shared_ptr<primitive> left, right;
	bvh box;

	bvhNode(shared_ptr<primitiveList> list) : bvhNode(list -> objects, 0, list -> objects.size()) {}
	bvhNode(const std::vector<shared_ptr<primitive>>& oldObjects, size_t start, size_t end) {
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
			left = make_shared<bvhNode>(objects, start, mid);
			right = make_shared<bvhNode>(objects, mid, end);
		}
		box = bvh(left->boundingBox(), right->boundingBox());
	}

	bool hit(const ray& r, interval t, hitRecord& rec) const override
	{
		if (!box.hit(r, t)) return false;
		bool hitLeft = left->hit(r, t, rec);
		bool hitRight = right->hit(r, interval(t.Min, hitLeft ? rec.t : t.Max), rec);
		return hitLeft || hitRight;
	}

	bvh boundingBox() const override { return box; }

private : 
	static bool box_compare(const shared_ptr<primitive> a, const shared_ptr<primitive> b, int axis_index) {
		return a -> boundingBox().axis(axis_index).Min < b -> boundingBox().axis(axis_index).Min;
	}

	static bool box_x_compare(const shared_ptr<primitive> a, const shared_ptr<primitive> b) {
		return box_compare(a, b, 0);
	}

	static bool box_y_compare(const shared_ptr<primitive> a, const shared_ptr<primitive> b) {
		return box_compare(a, b, 1);
	}

	static bool box_z_compare(const shared_ptr<primitive> a, const shared_ptr<primitive> b) {
		return box_compare(a, b, 2);
	}
};
