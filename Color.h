#pragma once
#ifndef _COLOR_H_
#define _COLOR_H_
#include "Vector.h"
#include <cmath>
#include<cassert>
/// <summary>
/// NOTE THAT this class is RGB, and doesn't contain alpha.
/// </summary>
struct Color
{
public:
	Color() = default;
	Color(Vector3 init_color) : color(init_color) {};
	double& r() { return color[0]; }
	double& g() { return color[1]; }
	double& b() { return color[2]; }

	/// <summary>
	/// This color should only be used in output.
	/// </summary>
	/// <returns>Real rgb color.</returns>
	static Color GammaCorrection(Vector3& radiance)
	{
		return Color({ _Correct(radiance.x()), _Correct(radiance.y()), _Correct(radiance.z()) });
	}
	Color operator+(const Color& color1)const { return Color({ color[0] + color1.color[0], color[1] + color1.color[1], color[2] + color1.color[2] }); }
	Color operator*(double coeff) const{ return Color({ color[0] * coeff, color[1] * coeff, color[2] * coeff }); }
	Color operator/(double coeff) const { double rep = 1 / coeff; return Color({ color[0] * rep, color[1] * rep, color[2] * rep });}
	Color& operator+=(const Color& color1) { color[0] += color1.color[0], color[1] += color1.color[1], color[2] += color1.color[2]; return *this; };
	Color& operator*=(double coeff) { color[0] *= coeff, color[1] *= coeff, color[2] *= coeff; return *this; }
	Color& operator/=(double coeff) { double rep = 1 / coeff; color[0] *= rep, color[1] *= rep, color[2] *= rep;  return *this; };

	Vector3 AbsorbedRadiance(const Vector3& radiance) const
	{
		Vector3 coeff = color / 255.0;
		return Vector3::ElemwiseMult(radiance, coeff);
	};

private:
	Vector3 color;
	static double _Clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
	static double _Correct(double c) {
		double clamped = _Clamp(c);
		double powed = std::pow(clamped, 1 / 2.2);
		return std::round(powed * 255.0);
	}
};
#endif // !_COLOR_H_
