#pragma once
#include <optional>
#include <memory>
#include <vector>
using std::shared_ptr;

class Light;
class LightSampleContext;

struct SampledLight {
	shared_ptr<Light> light;
	double p = 0;
};

class LightSampler
{
public: 
	virtual std::optional<SampledLight> Sample(const LightSampleContext& sample, double u) const = 0;
	virtual double PMF(const LightSampleContext& ctx, shared_ptr<Light> light) const = 0;
	virtual std::optional<SampledLight> Sample(double u) const = 0;
	virtual double PMF(shared_ptr<Light> light) const = 0;
};

class UniformLightSampler : public LightSampler
{
public:
	UniformLightSampler(const std::vector<shared_ptr<Light>>& _lights) : lights(_lights) {}

	std::optional<SampledLight> Sample(double u) const override {
		if (lights.empty()) return {};
		int lightIndex = std::min(static_cast<int>(u * lights.size()), static_cast<int>(lights.size()) - 1);
		return SampledLight { lights[lightIndex], 1.0 / lights.size() };
	}

	double PMF(shared_ptr<Light> light) const override {
		if (lights.empty()) return 0;
		return 1.0 / lights.size();
	}

	std::optional<SampledLight> Sample(const LightSampleContext& sample, double u) const override { return Sample(u); }
	double PMF(const LightSampleContext& ctx, shared_ptr<Light> light) const override { return PMF(light); }
private:
	std::vector<shared_ptr<Light>> lights;
};