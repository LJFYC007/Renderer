#pragma once
#include "aabb.h"
#include "interaction.h"
#include "material.h"

#include <optional>
#include <vector>
#include <memory>
using std::shared_ptr;
using std::make_shared;
using std::vector;

class material;

class hitRecord
{
public : 
	point3 p;
	vec3 normal, dpdu;
	double t, u, v;
	bool frontFace;
	shared_ptr<Material> mat;
};

class QuadricIntersection
{
public:
	double tHit;
	vec3 p;
};

class Shape
{
public : 
	virtual AABB Bounds() const = 0;
	virtual std::optional<hitRecord> Intersect(const ray& r, interval t) const = 0;
	//virtual bool IntersectP(const ray& r, interval t) const = 0;
	// virtual double Area() const = 0;
};

class Sphere : public Shape
{
public:
	Sphere(point3 _center, double _radius, shared_ptr<Material> _mat) : center(_center), radius(_radius), mat(_mat) {}

	AABB Bounds() const override {
		return AABB(center - vec3(radius), center + vec3(radius));
	}

	std::optional<hitRecord> Intersect(const ray& r, interval t) const override {
		std::optional<QuadricIntersection> ints = BasicIntersect(r, t);
		if (!ints) return {};

		hitRecord rec;
		rec.t = ints->tHit; rec.p = ints->p;
		rec.mat = mat;
		rec.normal = (rec.p - center) / radius;

		rec.frontFace = dot(rec.normal, r.rd) < 0.0;
		rec.normal = rec.frontFace ? rec.normal : -rec.normal;
		return rec;
	}

	std::optional<QuadricIntersection> BasicIntersect(const ray& r, interval t) const {
		vec3 oc = r.ro - center;
		double a = dot(r.rd, r.rd);
		double halfb = dot(oc, r.rd);
		double c = dot(oc, oc) - radius * radius;

		double discriminant = halfb * halfb - a * c;
		if (discriminant < 0.0) return {};
		discriminant = sqrt(discriminant);
		double root = (-halfb - discriminant) / a;
		if (!t.surrounds(root))
		{
			root = (-halfb + discriminant) / a;
			if (!t.surrounds(root)) return {};
		}
		return QuadricIntersection{ root, r.at(root) };
	}

private:
	vec3 center;
	double radius;
	shared_ptr<Material> mat;
};

class TriangleIntersection 
{
public:
	double t, u, v;
};

struct Vertex {
	vec3 p, n; vec2 uv;
	Vertex() {}
	Vertex(vec3 _p, vec3 _n) : p(_p), n(_n) {}
	Vertex(vec3 _p, vec3 _n, vec2 _uv) : p(_p), n(_n), uv(_uv) {}
};

class TriangleMesh
{
public:
	int nTriangles, nVertices;
	vector<int> vertexIndices;
	vector<Vertex> vertices;
	shared_ptr<Material> mat;

	TriangleMesh(const Transform& ObjectToWorld, vector<int> _vertexIndices, vector<Vertex> _vertices, shared_ptr<Material> _mat) :
		nTriangles(static_cast<int>(_vertexIndices.size()) / 3), nVertices(static_cast<int>(_vertices.size())), vertexIndices(_vertexIndices), mat(_mat) {
		vertices.resize(nVertices);
		for (int i = 0; i < nVertices; ++i)
		{
			vertices[i].p = ObjectToWorld(_vertices[i].p);
			vertices[i].n = ObjectToWorld(_vertices[i].n);
		}
	}
};
static vector<TriangleMesh> meshes;

class Triangle : public Shape
{
public:
	Triangle(int _meshIndex, int _triIndex) : meshIndex(_meshIndex), triIndex(_triIndex) {}

	const TriangleMesh* GetMesh() const {
		return &meshes[meshIndex];
	}

	AABB Bounds() const override {
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p;
		vec3 p1 = mesh->vertices[v[1]].p;
		vec3 p2 = mesh->vertices[v[2]].p;
		return AABB(p0, p1,p2);
	}

	std::optional<hitRecord> Intersect(const ray& r, interval t) const override {
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p;
		vec3 p1 = mesh->vertices[v[1]].p;
		vec3 p2 = mesh->vertices[v[2]].p;

		std::optional<TriangleIntersection> ints = BasicIntersect(r, t, p0, p1, p2);
		if (!ints) return {};

		hitRecord rec;
		rec.t = ints->t;
		rec.p = p0 + (p1 - p0) * ints->u + (p2 - p0) * ints->v;
		rec.dpdu = p1 - p0;
		rec.normal = mesh->vertices[v[0]].n
			+ (mesh->vertices[v[1]].n - mesh->vertices[v[1]].n) * ints->u
			+ (mesh->vertices[v[2]].n - mesh->vertices[v[1]].n) * ints->v;

		rec.mat = mesh->mat;
		return rec;
	}

	std::optional<TriangleIntersection> BasicIntersect(const ray& r, interval t, vec3 p0, vec3 p1, vec3 p2) const {
		vec3 e1 = p1 - p0, e2 = p2 - p0, s = r.ro - p0;
		vec3 s1 = cross(r.rd, e2), s2 = cross(s, e1);

		double coeff = 1.0 / dot(s1, e1);
		double T = coeff * dot(s2, e2);
		double U = coeff * dot(s1, s);
		double V = coeff * dot(s2, r.rd);

		if (!t.surrounds(T)) return {};

		if (T >= 0 && U >= 0 && V >= 0 && (1 - U - V) >= 0)
			return TriangleIntersection{ T, U, V };
		return {};
	}

private:
	int meshIndex = -1, triIndex = -1;
};
