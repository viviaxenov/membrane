#include<assert.h>
#include<stdio.h>
#include"../header/uniform.h"


//---------------------------------------------------------------------------------------
//---------------@UTask------------------------------------------------------------------
//---------------------------------------------------------------------------------------

UTask::UTask(double h, double tau, unsigned nodes, double c) : grid(nodes, nodes)
{
	_h = h;
	_tau = tau;
	_nodes = nodes;
	_c = c;

	if(_c <= 0.0)
	{
		printf("ERROR: sound speed <= 0!\n");
		assert(0);
	}

	_c = c;
	grid.SetRect(_h, _h);
	grid.OuterBorder();	
}

void UTask::StartCondPoint(unsigned x, unsigned y, double amp, double v)
{
	if(x >= _nodes || y >= _nodes)
	{
		printf("ERROR: point out of range!\n");
		assert(0);
	}

	
	grid[y][x].w[0] = amp;					// previous time step
	grid[y][x].w[1] = amp + _tau*v;				// current time step ()

}

void UTask::Iteration()
{
	double coef = _c*_tau/_h;
	coef *= coef;
	
	for(int i = 1; i < _nodes - 1; i++)
	{
		for(int j = 1; j < _nodes; j++)
		{
			grid[i][j].w[2] = grid[i + 1][j].w[1] + grid[i - 1][j].w[1] 
					+ grid[i][j + 1].w[1] + grid[i][j - 1].w[1]  
					- 4.0*grid[i][j].w[1];
			grid[i][j].w[2] *= coef;					// laplace oper
			grid[i][j].w[2] += 2.0*grid[i][j].w[1] - grid[i][j].w[0]; 
		}
	}
	for(int i = 1; i < _nodes - 1; i++)
	{
		for(int j = 1; j < _nodes; j++)
		{
			grid[i][j].w[0] = grid[i][j].w[1];				// swapping buffers
			grid[i][j].w[1] = grid[i][j].w[2];				// current -> prev
			grid[i][j].w[2] = 0.0;						// next -> current
		}
	}
}

void UTask::Dump(FILE *fp)
{
	for(int i = 0; i < _nodes; i++)
	{
		for(int j = 0; j < _nodes; j++)
		{
			fprintf(fp, "(%lg, %lg, %lg) ", grid[i][j].x_0(), grid[i][j].y_0(), grid[i][j].w[1]);
		}
		fprintf(fp, "\n");
	}
}

