#pragma once
#include "aabb.h"
#include "interaction.h"
#include "light.h"
#include "material.h"

#include <optional>
#include <vector>
#include <memory>
using std::shared_ptr;
using std::make_shared;
using std::vector;

struct ShapeIntersection
{
	SurfaceInteraction intr;
	double tHit;
};

struct ShapeSample
{
	Interaction intr;
	double pdf;
};

struct ShapeSampleContext
{
	ShapeSampleContext(Vector3fi _pi, vec3 _n, vec3 _ns, double _t) : pi(_pi), n(_n), ns(_ns), t(_t) {}
	vec3 p() const { return vec3(pi); }
	Vector3fi pi;
	vec3 n, ns;
	double t;
};

class Shape
{
public : 
	virtual AABB Bounds() const = 0;
	virtual std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const = 0;
	// virtual bool IntersectP(const ray& r, interval t) const = 0;
	virtual double Area() const = 0;
	virtual std::optional<ShapeSample> Sample(vec2 u) const = 0;
	virtual double PDF(ShapeSampleContext sample) const {
		return 1.0 / Area();
	}
	virtual std::optional<ShapeSample> Sample(ShapeSampleContext sample, vec2 u) const {
		return {};
	}
	virtual double PDF(ShapeSampleContext sample, vec3 wi) const {
		return 0.0;
	}
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

	TriangleMesh(const Transform& ObjectToWorld, vector<int> _vertexIndices, vector<Vertex> _vertices) :
		nTriangles(static_cast<int>(_vertexIndices.size()) / 3), nVertices(static_cast<int>(_vertices.size())), vertexIndices(_vertexIndices) {
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
		vec3 p0 = mesh->vertices[v[0]].p, p1 = mesh->vertices[v[1]].p, p2 = mesh->vertices[v[2]].p;
		return AABB(p0, p1,p2);
	}

	std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const override {
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p, p1 = mesh->vertices[v[1]].p, p2 = mesh->vertices[v[2]].p;

		std::optional<TriangleIntersection> ints = BasicIntersect(r, t, p0, p1, p2);
		if (!ints) return {};

		vec3 p = p0 + (p1 - p0) * ints->u + (p2 - p0) * ints->v;
		vec2 uv = mesh->vertices[v[0]].uv * (1 - ints->u - ints->v) +
			mesh->vertices[v[1]].uv * ints->u +
			mesh->vertices[v[2]].uv * ints->v;

		vec3 n0 = mesh->vertices[v[0]].n;
		vec3 n1 = mesh->vertices[v[1]].n;
		vec3 n2 = mesh->vertices[v[2]].n;

		vec2 uv0 = mesh->vertices[v[0]].uv;
		vec2 uv1 = mesh->vertices[v[1]].uv;
		vec2 uv2 = mesh->vertices[v[2]].uv;

		vec3 edge1 = p1 - p0;
		vec3 edge2 = p2 - p0;
		vec2 deltaUV1 = uv1 - uv0;
		vec2 deltaUV2 = uv2 - uv0;
		float invDet = 1.0f / (deltaUV1.x() * deltaUV2.y() - deltaUV1.y() * deltaUV2.x());
		vec3 tangent = invDet * (deltaUV2.y() * edge1 - deltaUV1.y() * edge2);
		vec3 bitangent = invDet * (-deltaUV2.x() * edge1 + deltaUV1.x() * edge2);
		vec3 interpolatedNormal = normalize((1 - ints->u - ints->v) * n0 + ints->u * n1 + ints->v * n2);
		vec3 dndu = cross(bitangent, interpolatedNormal);
		vec3 dndv = cross(interpolatedNormal, tangent);
		dndu = normalize(dndu);
		dndv = normalize(dndv);

		vec3 ti(1 - ints->u - ints->v, ints->u, ints->v);
		vec3 pAbsSum = Abs(ti.x() * p0) + Abs(ti.y() * p1) + Abs(ti.z() * p2);
		vec3 pError = gamma(7) * pAbsSum;

		SurfaceInteraction intr(Vector3fi(p, pError), uv, -r.rd, p1 - p0, p2 - p0, dndu, dndv, ints->t, false);
		return ShapeIntersection{ intr, ints->t };
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
	
	double Area() const override {
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p, p1 = mesh->vertices[v[1]].p, p2 = mesh->vertices[v[2]].p;
		return 0.5 * cross(p1 - p0, p2 - p0).length();
	}

	std::optional<ShapeSample> Sample(vec2 uv) const override { 
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p, p1 = mesh->vertices[v[1]].p, p2 = mesh->vertices[v[2]].p;

		vec3 b = SampleUniformTriangle(uv);
		vec3 p = b.x() * p0 + b.y() * p1 + b.z() * p2;

		vec3 n = normalize(cross(p1 - p0, p2 - p0)); 
		vec3 ns = b.x() * mesh->vertices[v[0]].n + b.y() * mesh->vertices[v[1]].n + b.z() * mesh->vertices[v[2]].n;
		n = FaceForward(n, ns);
		vec2 uvSample = b.x() * mesh->vertices[v[0]].uv + b.y() * mesh->vertices[v[1]].uv + b.z() * mesh->vertices[v[2]].uv;
		vec3 pAbsSum = Abs(b.x() * p0) + Abs(b.y() * p1) + Abs(b.z() * p2);
		vec3 pError = gamma(6) * pAbsSum;
		return ShapeSample{Interaction(Vector3fi(p, pError), n, uvSample), 1 / Area()};
	}

	std::optional<ShapeSample> Sample(ShapeSampleContext sample, vec2 u) const override {
		std::optional<ShapeSample> ss = Sample(u);
		if (!ss) return {};
		vec3 wi = ss->intr.p() - sample.p();
		if (wi.lengthSquared() == 0.0) return ShapeSample{ ss->intr, 0.0 };
		wi = normalize(wi);
		double pdf = ss->pdf * (sample.p() - ss->intr.p()).lengthSquared() / std::abs(dot(ss->intr.n, -wi));
		if (std::isinf(pdf)) pdf = 0.0;
		return ShapeSample{ ss->intr, pdf };
	}

	double PDF(ShapeSampleContext sample, vec3 wi) const override {
		ray r = SpawnRay(sample.pi, sample.n, wi);
		std::optional<ShapeIntersection> isect = Intersect(r, interval(0, infinity));
		if (!isect) return 0.0;
		SurfaceInteraction intr = isect->intr;
		double pdf = 1.0 * (sample.p() - intr.p()).lengthSquared() / std::abs(dot(intr.n, -wi)) / Area();
		if (std::isinf(pdf)) pdf = 0.0;
		return pdf;
	}

private:
	int meshIndex = -1, triIndex = -1;
};
