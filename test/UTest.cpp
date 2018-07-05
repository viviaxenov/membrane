#include"../header/membrane.h"
#include"../header/uniform.h"


int main()
{
	UTask test(0.1, 0.1, 5, 1.0);
	test.Print();
	puts("");
	test.Iteration();
	test.Print();
	return 0;
}
