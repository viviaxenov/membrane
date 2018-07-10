#include<stdio.h>
#include"../header/membrane.h"

int main()
{
	unsigned X = 4, Y = 3;

	Grid2D test(X, Y);
	test.SetRect(0.2, 0.1);
	test.SetMaterialUniform(1.0, 100000.0, 0.5);
	test.OuterBorder();

	for(int i = 0; i < Y; i++)
	{
		for(int j = 0; j < X; j++)
		{
			printf("(%2.02lg, %2.02lg) ", test[i][j].x_0(), test[i][j].y_0());	
		}
		puts("");
	}

	for(int i = 0; i < Y; i++)
	{
		for(int j = 0; j < X; j++)
		{
			double A = test[i][j].A();
			double B = test[i][j].B();
			double G = test[i][j].G();
			double rho = test[i][j].rho();
			double mu = B/A;
			double E = A*(1 - mu*mu);
			printf("A = %lg, B = %lg, G = %lg, E = %lg, mu = %lg, rho = %lg\n",
				A, B, G, E, mu, rho);
		}
	}


	for(int i = 0; i < Y; i++)
	{
		for(int j = 0; j < X; j++)
		{
			char c;
			switch(test[i][j].type())
			{
				case BORDER: 
					c = '*';
					break;
				case GRID:
					c = '_';
					break;
				case INACTIVE:
					c = '!';
			}
			putc(c, stdout);
		}
		puts("");
	}

	return 0;
}
