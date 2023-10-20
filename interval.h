#pragma once
#include <cmath>
#include <limits>

const double infinity = std::numeric_limits<double>::infinity();

class interval {
public : 
	double Min, Max;
	interval() : Min(+infinity), Max(-infinity) {}
	interval(double _min, double _max) : Min(_min), Max(_max) {}
	interval(const interval& a, const interval& b) : Min(fmin(a.Min, b.Min)), Max(fmax(a.Max, b.Max)) {}

	double size() const { return Max - Min; }
	bool contains(const double& x) const { return Min <= x && x <= Max; }
	bool surrounds(const double& x) const { return Min < x && x < Max; }

	double clamp(const double& x) const {
		if (x < Min) return Min; 
		if (x > Max) return Max; 
		return x;
	}

	interval expand(const double& delta) const {
		auto padding = delta / 2;
		return interval(Min - padding, Max + padding);
	}

	static const interval empty, universe;
};

const static interval empty(+infinity, -infinity);
const static interval universe(-infinity, +infinity);

