#ifndef MEMBRANE_INCL
#define MEMBRANE_INCL

#include<string>
//#include <vtkCellArray.h>
//#include <vtkPoints.h>
//#include <vtkStructuredGrid.h>
//#include <vtkSmartPointer.h>

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
	double _rho;					// density
	double _A, _B, _G;				// coeffs for counting \sigma(\epsilon)
		 					// A = \frac{E}{1 - \mu^2}
							// B = \mu A
							// G = \frac{E}{2(1 + \mu)}
							// where E - Young's modulus
							//	\mu - Poisson's ratio
	PointType _type;
public: 

	double u;					// offsets: 	u - Ox axis
	double v;					//		v - Oy axis
	double w;					//		w - Oz (vertical)
	
	double v_u;					// speeds: 	u - Ox axis
	double v_v;					//		v - Oy axis
	double v_w;					//		w - Oz (vertical)
							
	double dv_u;					// speed differentials
	double dv_v;					//		
	double dv_w;					//		


	

	Point(double x_0 = 0, double y_0 = 0,
			 PointType pt = INACTIVE);	// constructor	
	Point(Point&);					// copying constructor
	~Point(){};

	PointType type() {return _type;}		// getter functions
	double x_0() { return _x_0;}			
	double y_0() { return _y_0;}		
	double rho() { return _rho;}		
	double A() { return _A;}		
	double B() { return _B;}		
	double G() { return _G;}		

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

	void DiscardOffsets();				// sets u, v_u, dv_u etc to zero 

	bool IsCorrect();				

//	vtkSmartPointer<vtkStructuredGrid> vtkSGrid();
};

class Task
{
private: 
	Grid2D _grid;
	
	double tau;
	double h;
}


#endif /* MEMBRANE_INCL*/
