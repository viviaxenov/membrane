#ifndef VEC_INCLUDED
#define VEC_INCLUDED

#include<math.h>
#include<stdio.h>




class vec3
{
public: 
	double u;
	double v;
	double w;

	vec3(double u = 0.0, double v = 0.0, double w = 0.0) :
		u(u), v(v), w(w){}
	vec3(const vec3& orig) :
		u(orig.u), v(orig.v), w(orig.w){}

	~vec3(){};
	
	vec3& operator=(const vec3& right);
	vec3& operator+=(const vec3& right);
	vec3& operator-=(const vec3& right);
	vec3& operator*=(double c);

	void zero();
	
	void Print();

	static double DotProduct(const vec3& a, const vec3& b);
	static double Norm(const vec3& a);
	static double Projection(const vec3& a, const vec3& b);


};
vec3 operator+(const vec3& right, const vec3& left);
vec3 operator-(const vec3& right, const vec3& left);
vec3 operator*(double c, const vec3& left);



#endif /*VEC_INCLUDED */ 
