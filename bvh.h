#pragma once
#include "math.h"
#include "interval.h"
#include "ray.h"
#include "primitive.h"

#include <vector>
#include <memory>
#include <atomic>
using std::shared_ptr;

enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

struct BVHBuildNode;
struct alignas(32) LinearBVHNode
{
    AABB bounds;
    union { int primitiveOffset, secondChildOffset; };
    uint16_t nPrimitives;
    uint8_t axis;
};


struct BVHPrimitive
{
	BVHPrimitive() {}
	BVHPrimitive(int _primitiveNumber, const AABB& _bounds)
		: primitiveIndex(_primitiveNumber), bounds(_bounds) {}
	int primitiveIndex;
	AABB bounds;

	vec3 Centroid() const {	return 0.5 * (bounds.pMin + bounds.pMax); }
};

class BVHAggregate : public Primitive
{
public:
	BVHAggregate(std::vector<shared_ptr<Primitive>> _primitives, int _maxPrimsInNode, SplitMethod _splitMethod = SplitMethod::SAH)
		: primitives(_primitives), maxPrimesInNode(_maxPrimsInNode), splitMethod(_splitMethod) {
		bvhPrimitives.resize(primitives.size());
		for (int i = 0; i < primitives.size(); ++i)
			bvhPrimitives[i] = BVHPrimitive(i, primitives[i]->Bounds());

		std::vector<shared_ptr<Primitive>> orderedPrims(primitives.size());
		BVHBuildNode* root;
		std::atomic<int> totalNodes{ 0 };
		std::atomic<int> orderedPrimOffset{ 0 };
		root = buildRecursive(0, bvhPrimitives.size(), &totalNodes, &orderedPrimOffset, orderedPrims);
		primitives.swap(orderedPrims);

		bvhPrimitives.resize(0);
		nodes = new LinearBVHNode[totalNodes];
		int offset = 0;
		flattenBVH(root, &offset);
	}

	AABB Bounds() const override {
		return nodes[0].bounds;
	}

	std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const override {
		std::optional<ShapeIntersection> si;
		vec3 invDir(1.0 / r.rd[0], 1.0 / r.rd[1], 1.0 / r.rd[2]);
		int dirIsNeg[3] = { int(invDir[0] < 0), int(invDir[1] < 0), int(invDir[2] < 0) };
		int toVisitOffset = 0, currentNodeIndex = 0;
		int nodesToVisit[64];
		while (true) {
			const LinearBVHNode* node = &nodes[currentNodeIndex];
			if (node->bounds.Intersect(r.ro, r.rd, t.Max, invDir, dirIsNeg)) {
				if (node->nPrimitives > 0) {
					for (int i = 0; i < node->nPrimitives; ++i) {
						std::optional<ShapeIntersection> primSi = primitives[node->primitiveOffset + i]->Intersect(r, t);
						if (primSi) {
							si = primSi;
							t.Max = si->tHit;
						}
					}
					if (toVisitOffset == 0) break;
					currentNodeIndex = nodesToVisit[--toVisitOffset];
				}
				else {
					if (dirIsNeg[node->axis]) {
						nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
						currentNodeIndex = node->secondChildOffset;
					}
					else {
						nodesToVisit[toVisitOffset++] = node->secondChildOffset;
						currentNodeIndex = currentNodeIndex + 1;
					}
				}
			}
			else {
				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
		}
		return si;
	}

private:
	int maxPrimesInNode;
	std::vector<shared_ptr<Primitive>> primitives;
	std::vector<BVHPrimitive> bvhPrimitives;
	SplitMethod splitMethod;
	LinearBVHNode* nodes = nullptr;

	BVHBuildNode* buildRecursive(size_t begin, size_t end, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimOffset, std::vector<shared_ptr<Primitive>>& orderedPrims);
	int flattenBVH(BVHBuildNode* node, int* offset);
};

