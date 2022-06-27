#pragma once
#ifndef _MATHEXTENSION_H_
#define _MATHEXTENSION_H_

#include <tuple>
#include <cmath>
#include <numbers>

inline bool iszero(double f)
{
	return f < 1e-4 && f > -1e-4;
}

inline std::pair<double, double> sincos(double f)
{
	return { std::sin(f), std::cos(f) };
}
#endif // !_MATHEXTENSION_H_