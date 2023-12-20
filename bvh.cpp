#include "bvh.h"

struct BVHBuildNode 
{
    AABB bounds;
    BVHBuildNode* children[2];
    int splitAxis, firstPrimOffset, nPrimitives;

    void InitLeaf(int first, int n, const AABB& b)
    {
        firstPrimOffset = first;
        nPrimitives = n;
		bounds = b;
        children[0] = children[1] = nullptr;
    }

	void InitInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1)
	{
		children[0] = c0;
		children[1] = c1;
		bounds = AABB(c0->bounds, c1->bounds);
		splitAxis = axis;
		nPrimitives = 0;
	}
};

struct BVHSplitBucket {
	int count = 0;
	AABB bounds;
};

BVHBuildNode* BVHAggregate::buildRecursive(size_t begin, size_t end, std::atomic<int>* totalNodes, std::atomic<int>* orderedPrimOffset, std::vector<shared_ptr<Primitive>>& orderedPrims)
{
    BVHBuildNode* node = new BVHBuildNode();
    (*totalNodes)++;
    AABB bounds = bvhPrimitives[begin].bounds;
    for (size_t i = begin; i < end; ++i) bounds = AABB(bounds, bvhPrimitives[i].bounds);

    if (bounds.SurfaceArea() == 0 || begin + 1 == end) {
        int firstPrimOffset = orderedPrimOffset->fetch_add(static_cast<int>(end - begin));
        for (size_t i = begin; i < end; ++i) {
            size_t index = bvhPrimitives[i].primitiveIndex;
            orderedPrims[firstPrimOffset + i - begin] = primitives[index];
        }
        node->InitLeaf(firstPrimOffset, static_cast<int>(end - begin), bounds);
        return node;
    }

    AABB centroidBounds = bvhPrimitives[begin].Centroid();
    for (size_t i = begin; i < end; ++i) centroidBounds = AABB(centroidBounds, bvhPrimitives[i].Centroid());
    int dim = centroidBounds.MaxDimension();
    if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
        int firstPrimOffset = orderedPrimOffset->fetch_add(static_cast<int>(end - begin));
        for (size_t i = begin; i < end; ++i) {
            size_t index = bvhPrimitives[i].primitiveIndex;
            orderedPrims[firstPrimOffset + i - begin] = primitives[index];
        }
        node->InitLeaf(firstPrimOffset, static_cast<int>(end - begin), bounds);
        return node;
    }

    size_t mid;
    switch (splitMethod) {
    case SplitMethod::EqualCounts: {
        mid = (begin + end) / 2;
        std::nth_element(bvhPrimitives.begin() + begin, bvhPrimitives.begin() + mid, bvhPrimitives.begin() + end, [dim](const BVHPrimitive& a, const BVHPrimitive& b) {
            return a.Centroid()[dim] < b.Centroid()[dim];
            });
        break;
    }
    case SplitMethod::SAH: {
        if (begin + 2 >= end) {
            mid = (begin + end) / 2;
            std::nth_element(bvhPrimitives.begin() + begin, bvhPrimitives.begin() + mid, bvhPrimitives.begin() + end, [dim](const BVHPrimitive& a, const BVHPrimitive& b) {
                return a.Centroid()[dim] < b.Centroid()[dim];
                });
            break;
        }
        const int nBuckets = 12;
        BVHSplitBucket buckets[nBuckets];
        for (size_t i = begin; i < end; ++i) {
            int b = static_cast<int>(nBuckets * centroidBounds.Offset(bvhPrimitives[i].Centroid())[dim]);
            if (b == nBuckets) b = nBuckets - 1;
            buckets[b].count++;
            buckets[b].bounds = AABB(buckets[b].bounds, bvhPrimitives[i].bounds);
        }

        const int nSplits = nBuckets - 1;
        double costs[nSplits] = {};
        int countBelow = 0; AABB boundBelow;
        for (int i = 0; i < nSplits; ++i) {
            boundBelow = AABB(boundBelow, buckets[i].bounds);
            countBelow += buckets[i].count;
            costs[i] += countBelow * boundBelow.SurfaceArea();
        }
        int countAbove = 0; AABB boundAbove;
        for (int i = nSplits; i >= 1; --i) {
            boundAbove = AABB(boundAbove, buckets[i].bounds);
            countAbove += buckets[i].count;
            costs[i - 1] += countAbove * boundAbove.SurfaceArea();
        }

        int minCostSplitBucket = -1; double minCost = std::numeric_limits<double>::max();
        for (int i = 0; i < nSplits; ++i) if (costs[i] < minCost) {
            minCost = costs[i];
            minCostSplitBucket = i;
        }
        double leafCost = static_cast<double>(end - begin);
        minCost = 0.5 + minCost / bounds.SurfaceArea();

        if (end - begin > maxPrimesInNode || minCost < leafCost) {
            auto midIter = std::partition(bvhPrimitives.begin() + begin, bvhPrimitives.begin() + end, [dim, minCostSplitBucket, nBuckets, centroidBounds](const BVHPrimitive& a) {
                int b = static_cast<int>(nBuckets * centroidBounds.Offset(a.Centroid())[dim]);
                if (b == nBuckets) b = nBuckets - 1;
                return b <= minCostSplitBucket;
                });
            mid = midIter - bvhPrimitives.begin();
        }
        else {
            int firstPrimOffset = orderedPrimOffset->fetch_add(static_cast<int>(end - begin));
            for (size_t i = begin; i < end; ++i) {
                size_t index = bvhPrimitives[i].primitiveIndex;
                orderedPrims[firstPrimOffset + i - begin] = primitives[index];
            }
			node->InitLeaf(firstPrimOffset, static_cast<int>(end - begin), bounds);
            return node;
        }
        break;
    }
    }

    BVHBuildNode* children[2];
    /*
    if (bvhPrimitives.size() > 128 * 1024) {
    }
    else {
    }
    */
    children[0] = buildRecursive(begin, mid, totalNodes, orderedPrimOffset, orderedPrims);
    children[1] = buildRecursive(mid, end, totalNodes, orderedPrimOffset, orderedPrims);
    node->InitInterior(dim, children[0], children[1]);
    return node;
}

int BVHAggregate::flattenBVH(BVHBuildNode* node, int* offset)
{
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    int nodeOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        linearNode->primitiveOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVH(node->children[0], offset);
		linearNode->secondChildOffset = flattenBVH(node->children[1], offset);
    }
    return nodeOffset;
}
