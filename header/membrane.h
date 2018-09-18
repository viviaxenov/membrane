#ifndef MEMBRANE_INCL
#define MEMBRANE_INCL

#include<string>


#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>

#include"../header/vec3.h"
#include"../header/constants.h"

using std::string;

enum PointType 
{
	TL_CORNER,
	TR_CORNER,
	BR_CORNER, 
	BL_CORNER,
	B_BORDER,
	T_BORDER,
	L_BORDER,
	R_BORDER,
	GRID,
	INACTIVE
};

class Point
{
private:
	double _x_0;					// initial positions
	double _y_0;		 
public: 
	PointType type;
	double rho;					// density
	double delta;					// thickness
	double A, B, G;					// coeffs for counting \sigma(\epsilon)
		 					// A = \frac{E}{1 - \mu^2}
							// B = \mu A
							// G = \frac{E}{2(1 + \mu)}
							// where E - Young's modulus
							//	\mu - Poisson's ratio

	vec3 r;						// offsets: 	u - Ox axis

							//		v - Oy axis
							//		w - Oz (vertical)

	vec3 v;						// velocity 
							
	vec3 a;					// velocity differentials
	vec3 F_ext;					// external force
	double sigma_xx, sigma_yy, sigma_xy;
	double w_x, w_y;				// dw/dx and dw/dy needed to compute tangential basis vectors

	Point(double x_0 = 0, double y_0 = 0,
			 PointType pt = INACTIVE);	// constructor	
	Point(Point&);					// copying constructor
	~Point(){};

	double x_0() { return _x_0;}			
	double y_0() { return _y_0;}		

	void Set(double x_0, double y_0, PointType tp = GRID);
							// set point's place and type
	void SetMaterial(double rho, double E, double mu, double delta);
							// set elastic properties

	string ToString();				// dumps _x_0, _y_0, type
	
	bool IsCorrect();				
};


class Grid2D
{
private: 
	unsigned _x_nodes;				// size
	unsigned _y_nodes;

	double _h_x;					// distances between grid cells; 
							// possibly will be eliminated when arbitrary grid geometry is done
	double _h_y;

	Point *_grid;					// data array
public: 
	Grid2D(unsigned x_nodes = 0, unsigned y_nodes = 0);
//	Grid2D(Grid2D&);
	~Grid2D(){};
	
	unsigned x_nodes() {return _x_nodes;}		
	unsigned y_nodes() {return _y_nodes;}		


	Point* operator [](int y)			// data access in a[y][x] manner
	{
		return _grid + y*_x_nodes;
	}

	void SetRect(double dx , double dy);		// sets grid geometry to rectangle with steps dx and dy
	void SetMaterialUniform(double rho, double E, double mu, double delta);
	void OuterBorder();				// sets type of border cells to border correctly
							// ****	|
							// *00*	| * - border nodes
							// *00*	| 0 - inner (GRID) nodes
							// ****	|

	void DiscardOffsets();				// sets r, v, a etc to zero 

	bool IsCorrect();				
	vec3 delta_x(unsigned x, unsigned y);
	vec3 delta_y(unsigned x, unsigned y);

	void CountSigmas();
	void CountAcceleration();
};



class Task
{
private: 
	Grid2D _grid;
	
	double _tau;					// time step
	double _delta;					// membrane thickness
							// it's here temporarily, assume it to be 
							// the same across the whole membrane

	string _output_dir;
public:

	
	unsigned frames, iter_per_frame;

	Task(double tau, double h, unsigned cells);
	int SetOutputDir(string of);
	void SetFextPt(unsigned x, unsigned y, 
			double F_u, double F_v, double F_w);
	void SetOffPt(unsigned x, unsigned y, 
			double u, double v, double w);
	
	void Iteration();
	int Execute();

	vtkSmartPointer<vtkStructuredGrid> vtkSGrid();
};


#endif /* MEMBRANE_INCL*/
