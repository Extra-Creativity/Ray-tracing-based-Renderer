#pragma once
#ifndef _VECTOR_H_
#define _VECTOR_H_
#include "MathExtension.h"

#include <array>
#include <ostream>
#include <format>


struct Vector3
{
	friend struct Vector4;
public:
	Vector3() = default;
	Vector3(const double init_x, const double init_y, const double init_z) :vec({ init_x, init_y, init_z }) {};
	double& x() { return vec[0]; };
	double& y() { return vec[1]; };
	double& z() { return vec[2]; };

	double x() const { return vec[0]; };
	double y() const { return vec[1]; };
	double z() const { return vec[2]; };

	// calculate with coeff
	Vector3 operator*(double coeff) const { return Vector3(vec[0] * coeff, vec[1] * coeff, vec[2] * coeff); };
	Vector3 operator/(double coeff) const { double rep = 1.0f / coeff; return Vector3(vec[0] * rep, vec[1] * rep, vec[2] * rep); };
	Vector3& operator*=(double coeff) { vec[0] *= coeff; vec[1] *= coeff; vec[2] *= coeff; return *this; }
	Vector3& operator/=(double coeff) { double rep = 1.0f / coeff; vec[0] *= rep; vec[1] *= rep; vec[2] *= rep; return *this; }

	// calculate with Vector3, without ambiguity.
	Vector3 operator+(const Vector3& vec1) const  { return Vector3(vec[0] + vec1.vec[0], vec[1] + vec1.vec[1], vec[2] + vec1.vec[2]); }
	Vector3 operator-(const Vector3& vec1) const { return Vector3(vec[0] - vec1.vec[0], vec[1] - vec1.vec[1], vec[2] - vec1.vec[2]); }
	Vector3& operator+=(const Vector3& vec1) { vec[0] += vec1.vec[0]; vec[1] += vec1.vec[1]; vec[2] += vec1.vec[2]; return *this; }
	Vector3& operator-=(const Vector3& vec1) { vec[0] -= vec1.vec[0]; vec[1] -= vec1.vec[1]; vec[2] -= vec1.vec[2]; return *this; }
	Vector3 operator-() const { return Vector3(-vec[0], -vec[1], -vec[2]); };

	// calculate with Vector3, with ambiguity.
	void ElemwiseMult(const Vector3& vec1) { vec[0] *= vec1.vec[0]; vec[1] *= vec1.vec[1]; vec[2] *= vec1.vec[2]; }

	// use for normal calulate.
	static double Dot(const Vector3& vec1, const Vector3& vec2) { return vec1.vec[0] * vec2.vec[0] + vec1.vec[1] * vec2.vec[1] + vec1.vec[2] * vec2.vec[2]; }
	static Vector3 Cross(const Vector3& vec1, const Vector3& vec2) {
		return Vector3(vec1.vec[1] * vec2.vec[2] - vec2.vec[1] * vec1.vec[2], vec1.vec[2] * vec2.vec[0] - vec2.vec[2] * vec1.vec[0], vec1.vec[0] * vec2.vec[1] - vec2.vec[0] * vec1.vec[1]);
	};
	static Vector3 ElemwiseMult(const Vector3& vec1, const Vector3& vec2) { return Vector3(vec1.vec[0] * vec2.vec[0], vec1.vec[1] * vec2.vec[1], vec1.vec[2] * vec2.vec[2]);}
	// index op
	double& operator[](int index) { return vec[index]; };
	double operator[](int index) const { return vec[index]; }

	static Vector3 zero() { return Vector3(0, 0, 0); }
	static bool isZero(const Vector3& vec1) { return iszero(vec1.vec[0]) && iszero(vec1.vec[1]) && iszero(vec1.vec[2]); }

	void Normalize() { double length = Length(); if (length < 1e-5) return; double lengthRep = 1 / length; vec[0] *= lengthRep, vec[1] *= lengthRep, vec[2] *= lengthRep; };
	double Magnitude() const { return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]; };
	double Length() const { return std::sqrt(Magnitude()); };

	static Vector3 Normalized(const Vector3& vec) { 
		double length = vec.Length(); 
		if (length < 1e-5) 
			return zero(); 
		double lengthRep = 1 / length; 
		
		Vector3 result{ vec.vec[0] * lengthRep,  vec.vec[1] * lengthRep, vec.vec[2] * lengthRep };
		return result;
	};

	bool operator==(const Vector3& vec1)const { return iszero(vec[0] - vec1.vec[0]) && iszero(vec[1] - vec1.vec[1]) && iszero(vec[2] - vec1.vec[2]); }
	bool operator!=(const Vector3& vec1)const { return !(this->operator==(vec1)); }

	void Print(std::ostream& os) const  { os << std::format("[{:6f},{:6f},{:6f}]\n", vec[0], vec[1], vec[2]); }
private:
	std::array<double, 3> vec;
};

inline Vector3 operator*(double coeff, const Vector3& vec) { return vec * coeff; }

struct Vector4
{
public:
	Vector4() = default;
	Vector4(double init_x, double init_y, double init_z, double init_w) :vec({ init_x, init_y, init_z, init_w }) {};
	Vector4(Vector3 init_vec, double init_w) :vec({ init_vec.vec[0], init_vec.vec[1], init_vec.vec[2], init_w }) {};
	double& x() { return vec[0]; };
	double& y() { return vec[1]; };
	double& z() { return vec[2]; };
	double& w() { return vec[3]; };

	double x() const { return vec[0]; };
	double y() const { return vec[1]; };
	double z() const { return vec[2]; };
	double w() const { return vec[3]; };

	// calculate with coeff
	Vector4 operator*(double coeff) const { return Vector4(vec[0] * coeff, vec[1] * coeff, vec[2] * coeff, vec[3] * coeff); };
	Vector4 operator/(double coeff) const { double rep = 1.0f / coeff; return Vector4(vec[0] * rep, vec[1] * rep, vec[2] * rep, vec[3] * rep); };

	// calculate with Vector4, without ambiguity.
	Vector4 operator+(const Vector4& vec1) const { return Vector4(vec[0] + vec1.vec[0], vec[1] + vec1.vec[1], vec[2] + vec1.vec[2], vec[3] + vec1.vec[3]); }
	Vector4 operator-(const Vector4& vec1) const { return Vector4(vec[0] - vec1.vec[0], vec[1] - vec1.vec[1], vec[2] - vec1.vec[2], vec[3] - vec1.vec[3]); }
	Vector4& operator+=(const Vector4& vec1) { vec[0] += vec1.vec[0]; vec[1] += vec1.vec[1]; vec[2] += vec1.vec[2]; vec[3] += vec1.vec[3]; return *this; }
	Vector4& operator-=(const Vector4& vec1) { vec[0] -= vec1.vec[0]; vec[1] -= vec1.vec[1]; vec[2] -= vec1.vec[2]; vec[3] -= vec1.vec[3]; return *this; }

	// calculate with Vector4, with ambiguity.
	void ElemwiseMult(const Vector4& vec1) { vec[0] *= vec1.vec[0]; vec[1] *= vec1.vec[1]; vec[2] *= vec1.vec[2]; vec[3] *= vec1.vec[3]; }

	// use for normal calulate.
	static double Dot(const Vector4& vec1, const Vector4& vec2) { return vec1.vec[0] * vec2.vec[0] + vec1.vec[1] * vec2.vec[1] + vec1.vec[2] * vec2.vec[2] + vec1.vec[3] * vec2.vec[3]; }
	static Vector4 ElemwiseMult(const Vector4& vec1, const Vector4& vec2) { return Vector4(vec1.vec[0] * vec2.vec[0], vec1.vec[1] * vec2.vec[1], vec1.vec[2] * vec2.vec[2], vec1.vec[3] * vec2.vec[3]); }
	// index op
	double& operator[](int index) { return vec[index]; };
	double operator[](int index) const { return vec[index]; }
	static Vector4 zero() { return Vector4(0, 0, 0, 0); }
	
	void Print(std::ostream& os) const { os << std::format("[{:6f},{:6f},{:6f},{:6f}]\n", vec[0], vec[1], vec[2], vec[3]); }
	/// <summary>
	/// If this Vector4 has a w approximate to 0, we will return (x, y, z); else return (x / w, y / w, z / w).
	/// </summary>
	/// <returns>Vector3</returns>
	Vector3 ToVector3() {
		double temp = vec[3];
		if (temp < 1e-5f && temp > 1e-5f)
			return Vector3(vec[0], vec[1], vec[2]);
		else
		{
			double tempp = 1 / temp;
			return Vector3(vec[0] * tempp, vec[1] * tempp, vec[2] * tempp);
		}
	}

	void Normalize() { double length = Length(); if (length < 1e-5) return; double lengthRep = 1 / length; vec[0] *= lengthRep, vec[1] *= lengthRep, vec[2] *= lengthRep ,vec[3] *= lengthRep; };
	double Magnitude() const { return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] + vec[3] * vec[3]; };
	double Length() const { return std::sqrt(Magnitude()); }

	static Vector4 Normalized(const Vector4& vec) {
		double length = vec.Length();
		if (length < 1e-5)
			return zero();
		double lengthRep = 1 / length;

		Vector4 result{ vec.vec[0] * lengthRep,  vec.vec[1] * lengthRep, vec.vec[2] * lengthRep, vec.vec[3] * lengthRep };
		return result;
	};

	bool operator==(const Vector4& vec1)const { return iszero(vec[0] - vec1.vec[0]) && iszero(vec[1] - vec1.vec[1]) && iszero(vec[2] - vec1.vec[2]) && iszero(vec[3] - vec1.vec[3]); }
	bool operator!=(const Vector4& vec1)const { return !(this->operator==(vec1)); }

private:
	std::array<double, 4> vec;
};

inline Vector4 operator*(double coeff, const Vector4& vec) { return vec * coeff; }
#endif // !_VECTOR_H_