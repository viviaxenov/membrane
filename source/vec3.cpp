#include"../header/vec3.h"


vec3& vec3::operator=(const vec3& right)
{
	u = right.u;
	v = right.v;
	w = right.w;
	return *this;
}
vec3& vec3::operator+=(const vec3& right)
{
	u += right.u;
	v += right.v;
	w += right.w;
	return *this;
}
vec3& vec3::operator-=(const vec3& right)
{
	u -= right.u;
	v -= right.v;
	w -= right.w;
	return *this;
}
vec3& vec3::operator*=(double c)
{
	u *= c;
	v *= c;
	w *= c;
	return *this;
}
vec3 operator+(const vec3& left, const vec3& right)
{
	vec3 res;
	res += right;
	res += left;
	return res;
}

vec3 operator-(const vec3& left, const vec3& right)
{
	vec3 res;
	res -= right;
	res += left;
	return res;
}
vec3 operator*(double c, const vec3& right)
{
	vec3 res(right);
	res *= c;
	return res;
}

void vec3::Print()
{
	printf("(%lg, %lg, %lg)\n", u, v, w);
}

double vec3::DotProduct(const vec3& a, const vec3& b)
{
	return a.u*b.u + a.v*b.v + a.w*b.w; 
}

double vec3::Norm(const vec3& a)
{
	return sqrt(vec3::DotProduct(a,a));
}
double vec3::Projection(const vec3& a, const vec3& b)
{
	return vec3::DotProduct(a, b)/vec3::Norm(b);
}



