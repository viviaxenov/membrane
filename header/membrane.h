#ifndef MEMBRANE_INCL
#define MEMBRANE_INCL

#include<string>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

#include"../header/vec3.h"
#include"../header/constants.h"

using std::string;

enum PointType 
{
	BORDER,
	GRID,
	INACTIVE
};

class Point
{
private:
	double _x_0;					// initial positions
	double _y_0;		 
	PointType _type;
public: 
	double rho;					// density
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
							
	vec3 dv;					// velocity differentials
	vec3 F_ext;					// external force
	double sigma_xx, sigma_yy, sigma_xy;

	Point(double x_0 = 0, double y_0 = 0,
			 PointType pt = INACTIVE);	// constructor	
	Point(Point&);					// copying constructor
	~Point(){};

	PointType type() {return _type;}		// getter functions
	double x_0() { return _x_0;}			
	double y_0() { return _y_0;}		

	void Set(double x_0, double y_0, PointType tp = GRID);
							// set point's place and type
	void SetMaterial(double rho, double E, double mu);
							// set elastic properties

	void ToBorder();				// changes _type to BORDER

	string ToString();				// dumps _x_0, _y_0, type
	
	bool IsCorrect();				
};


class Grid2D
{
private: 
	unsigned _x_nodes;				// size
	unsigned _y_nodes;

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
	void SetMaterialUniform(double rho, double E, double mu);
	void OuterBorder();				// sets type of border cells to border correctly
							// ****	|
							// *00*	| * - border nodes
							// *00*	| 0 - inner (GRID) nodes
							// ****	|

	void DiscardOffsets();				// sets r, v, dv etc to zero 

	bool IsCorrect();				

};



class Task
{
private: 
	Grid2D _grid;
	
	double _tau;					// time step
	double _h;					// distance btw greed points
	// assume membrane to be rectangular; will make other geometry
	// when file input is done
	double _delta;					// membrane thickness
							// it's here temporarily, assume it to be 
							// the same across the whole membrane

	void CountSigmas();
	void DvInner();
	void DvTopBot();
	void DvRightLeft();
	void DvCorners();
public:
	Task(double tau, double h, double delta, unsigned cells);
	void SetFextPt(unsigned x, unsigned y, 
			double F_u, double F_v, double F_w);
	void SetOffPt(unsigned x, unsigned y, 
			double u, double v, double w);
	
	void Iteration();

	vtkSmartPointer<vtkStructuredGrid> vtkSGrid();
};


#endif /* MEMBRANE_INCL*/
