#include<stdio.h>
#include"../header/vec3.h"

int main()
{

	vec3 a(3, 4, 5);
	
	vec3 e1(1, 0, 0), e2(0, 1, 0), e3(0, 0, 1);

	printf("%lg %lg %lg\n", 
		vec3::Projection(a,e1),
		vec3::Projection(a,e2),
		vec3::Projection(a,e3)
		);

	return 0;
}


