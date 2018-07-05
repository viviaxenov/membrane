#include<stdio.h>
#include"../header/membrane.h"

class UTask
{
private:
	Grid2D grid;			// data array

	double _h;			// node distance
	double _tau;			// time step
	unsigned _nodes;		// number of x and y nodes
					// total size _nodes^2

	double _c;			// speed of sound
public:
	UTask(double h, double tau, unsigned nodes, double c);
	~UTask(){};

	void StartCondPoint(unsigned x, unsigned y, double amp, double v);
	void Iteration();
	void Dump(FILE *fp);
	void Print() {this->Dump(stdout);}
	vtkSmartPointer<vtkStructuredGrid> vtkSGrid();
};
