#include<stdio.h>
#include"../header/membrane.h"

int main()
{
	unsigned X = 4, Y = 3;

	Grid2D test(X, Y);
	test.SetRect(0.2, 0.1);
	test.SetMaterialUniform(1.0, 100000.0, 0.5, 1e-3);
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
			double A = test[i][j].A;
			double B = test[i][j].B;
			double G = test[i][j].G;
			double rho = test[i][j].rho;
			double mu = B/A;
			double E = A*(1 - mu*mu);
			printf("A = %lg, B = %lg, G = %lg, E = %lg, mu = %lg, rho = %lg\n",
				A, B, G, E, mu, rho);
		}
	}


	for(int i = Y - 1; i >= 0; i--)
	{
		for(int j = 0; j < X; j++)
		{
			char c;
			switch(test[i][j].type)
			{
				case T_BORDER:
					c = '-';
					break;
				case B_BORDER: 
					c = '_';
					break;
				case GRID:
					c = 'o';
					break;
				case L_BORDER:
					c = '[';
					break;
				case R_BORDER:
					c = ']';
					break;
				case TL_CORNER:
					c = 'T';
					break;
				case TR_CORNER:
					c = '7';
					break;
				case BR_CORNER:
					c = 'J';
					break;
				case BL_CORNER:
					c = 'L';
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
