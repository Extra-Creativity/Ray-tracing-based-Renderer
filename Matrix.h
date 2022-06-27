#pragma once
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "Vector.h"
class Matrix3x3
{
public:
	Matrix3x3() = default;
	Matrix3x3(Vector3 Col1, Vector3 Col2, Vector3 Col3)
	{
		for (int i = 0; i < 3; i++)
		{
			mat[i][0] = Col1[i];
			mat[i][1] = Col2[i];
			mat[i][2] = Col3[i];
		};
	}

	static Matrix3x3 zero()
	{
		Matrix3x3 result{};
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				result.mat[i][j] = 0;
			}
		return result;
	}

	static Matrix3x3 identity()
	{
		Matrix3x3 result{};
		result.mat[0][0] = result.mat[1][1] = result.mat[2][2] = 1;
		result.mat[1][0] = result.mat[2][0] = result.mat[0][1] = 0;
		result.mat[0][2] = result.mat[1][2] = result.mat[2][1] = 0;
		return result;
	}

	Matrix3x3& operator+=(Matrix3x3 mat1)
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				mat[i][j] += mat1.mat[i][j];
			}
		return *this;
	}

	Matrix3x3& operator-=(Matrix3x3 mat1)
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				mat[i][j] -= mat1.mat[i][j];
			}
		return *this;
	}

	Matrix3x3& operator*=(double coeff)
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				mat[i][j] *= coeff;
			}
		return *this;
	}

	Matrix3x3& operator/(double coeff)
	{
		double rep = 1 / coeff;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				mat[i][j] *= rep;
			}
		return *this;
	}

	Vector3 operator*(Vector3 vec) {
		return Vector3(mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2],
			mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2],
			mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2]);
	}

	Matrix3x3 operator*(Matrix3x3 mat1)
	{
		Matrix3x3 result = zero();
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
				{
					result.mat[i][j] += mat[i][k] * mat1.mat[k][j];
				}
		return result;
	}

	Matrix3x3 operator+(Matrix3x3 mat1)
	{
		Matrix3x3 result{};
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				result.mat[i][j] = mat[i][j] + mat1.mat[i][j];
		return result;
	}

	Matrix3x3 operator-(Matrix3x3 mat1)
	{
		Matrix3x3 result{};
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				result.mat[i][j] = mat[i][j] - mat1.mat[i][j];
		return result;
	}

	Vector3 Col(int index)
	{
		return Vector3(mat[0][index], mat[1][index], mat[2][index]);
	}

	double Det()
	{
		return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
			mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]) +
			mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
	}

	Matrix3x3 Inverse()
	{
		Matrix3x3 result{};
		result.mat[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
		result.mat[0][1] = -mat[0][1] * mat[2][2] + mat[2][1] * mat[0][2];
		result.mat[0][2] = mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2];
		result.mat[1][0] = -mat[1][0] * mat[2][2] + mat[1][2] * mat[2][0];
		result.mat[1][1] = mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2];
		result.mat[1][2] = -mat[0][0] * mat[1][2] + mat[1][0] * mat[0][2];
		result.mat[2][0] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
		result.mat[2][1] = -mat[0][0] * mat[2][1] + mat[0][1] * mat[2][0];
		result.mat[2][2] = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];

		double det = mat[0][0] * result.mat[0][0] + mat[0][1] * result.mat[1][0] + mat[0][2] * result.mat[2][0];
		result *= 1 / det;
		return result;
	}

	Matrix3x3 Transpose()
	{
		Matrix3x3 result{};
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				result.mat[i][j] = mat[j][i];
		return result;
	}

	double Trace()
	{
		return mat[0][0] + mat[1][1] + mat[2][2];
	}

	Vector3 Row(int index) {return mat[index];}

	Vector3& operator[](int index) { return mat[index]; }

	void Print(std::ostream& os)
	{
		os << std::format("[{:6f},{:6f},{:6f},\n {:6f},{:6f},{:6f},\n {:6f},{:6f},{:6f}]\n", mat[0][0], mat[0][1], mat[0][2], mat[1][0], mat[1][1], mat[1][2], mat[2][0], mat[2][1], mat[2][2]);
		return;
	}
private:
	std::array<Vector3, 3> mat;
};

class Matrix4x4
{
public:
	Matrix4x4() = default;
	Matrix4x4(Vector4 Col1, Vector4 Col2, Vector4 Col3, Vector4 Col4)
	{
		for (int i = 0; i < 3; i++)
		{
			mat[i][0] = Col1[i];
			mat[i][1] = Col2[i];
			mat[i][2] = Col3[i];
			mat[i][3] = Col4[i];
		};
	}

	static Matrix4x4 zero()
	{
		Matrix4x4 result{};
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
			{
				result.mat[i][j] = 0;
			}
		return result;
	}

	static Matrix4x4 identity()
	{
		Matrix4x4 result{};
		result.mat[0][0] = result.mat[1][1] = result.mat[2][2] = result.mat[3][3] = 1;
		result.mat[1][0] = result.mat[2][0] = result.mat[0][1] = result.mat[3][0] = 0;
		result.mat[0][2] = result.mat[1][2] = result.mat[2][1] = result.mat[3][2] = 0;
		result.mat[3][1] = result.mat[1][3] = result.mat[0][3] = result.mat[2][3] = 0;
		return result;
	}

	Matrix4x4& operator+=(Matrix4x4 mat1)
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
			{
				mat[i][j] += mat1.mat[i][j];
			}
		return *this;
	}

	Matrix4x4& operator-=(Matrix4x4 mat1)
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
			{
				mat[i][j] -= mat1.mat[i][j];
			}
		return *this;
	}

	Matrix4x4& operator*=(double coeff)
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
			{
				mat[i][j] *= coeff;
			}
		return *this;
	}

	Matrix4x4& operator/(double coeff)
	{
		double rep = 1 / coeff;
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
			{
				mat[i][j] *= rep;
			}
		return *this;
	}

	Vector4 operator*(Vector4 vec) {
		return Vector4(mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2] + mat[0][3] * vec[3],
			mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2] + mat[1][3] * vec[3],
			mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2] + +mat[2][3] * vec[3],
			mat[3][0] * vec[0] + mat[3][1] * vec[1] + mat[3][2] * vec[2] + +mat[3][3] * vec[3]);
	}

	Matrix4x4 operator*(Matrix4x4 mat1)
	{
		Matrix4x4 result = zero();
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int k = 0; k < 4; k++)
				{
					result.mat[i][j] += mat[i][k] * mat1.mat[k][j];
				}
		return result;
	}

	Matrix4x4 operator+(Matrix4x4 mat1)
	{
		Matrix4x4 result{};
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				result.mat[i][j] = mat[i][j] + mat1.mat[i][j];
		return result;
	}

	Matrix4x4 operator-(Matrix4x4 mat1)
	{
		Matrix4x4 result{};
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				result.mat[i][j] = mat[i][j] - mat1.mat[i][j];
		return result;
	}

	Vector4 Col(int index)
	{
		return Vector4(mat[0][index], mat[1][index], mat[2][index], mat[3][index]);
	}

	Matrix4x4 Transpose()
	{
		Matrix4x4 result{};
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				result.mat[i][j] = mat[j][i];
		return result;
	}

	double Trace()
	{
		return mat[0][0] + mat[1][1] + mat[2][2] + mat[3][3];
	}

	Vector4 Row(int index) { return mat[index]; }

	Vector4& operator[](int index) { return mat[index]; }

	void Print(std::ostream& os)
	{
		os << std::format("[{:6f},{:6f},{:6f},{:6f},\n {:6f},{:6f},{:6f},{:6f},\n {:6f},{:6f},{:6f},{:6f},\n {:6f},{:6f},{:6f},{:6f}]\n", mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0], mat[1][1], mat[1][2], \
			mat[1][3], mat[2][0], mat[2][1], mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
		return;
	}
private:
	std::array<Vector4, 4> mat;
};

#endif // !_MATRIX_H_