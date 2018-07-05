#ifndef MEMBRANE_INCL
#define MEMBRANE_INCL

#include<string>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

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

	double u[3];					// offsets: 	u - Ox axis
	double v[3];					//		v - Oy axis
	double w[3];					//		w - Oz (vertical)
							// 3 values for previous, current and next time steps

//	double _sigma_xx, _sigma_yy, _sigma_zz;		// stress tensor
//	double _tau_xy, _tau_xz, _tau_yz;
	

	Point(double x_0 = 0, double y_0 = 0,
			 PointType pt = INACTIVE);	// constructor	
	Point(Point&);					// copying constructor
	~Point(){};

	PointType type() {return _type;}
	double x_0() { return _x_0;}			
	double y_0() { return _y_0;}		

	void Set(double x_0, double y_0, PointType tp = GRID);
	void ToBorder();				// changes _type to BORDER

	string ToString();				// dumps _x_0, _y_0, type
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

	Point* operator [](int y)			// data access in a[y]x[] manner
	{
		return _grid + y*_x_nodes;
	}

	void SetRect(double dx , double dy);		// sets grid geometry to rectangle with steps dx and dy
	void OuterBorder();				// sets type of border cells to border correctly

	vtkSmartPointer<vtkStructuredGrid> vtkSGrid();
							// ****	|
							// *00*	| * - border nodes
							// *00*	| 0 - inner (GRID) nodes
							// ****	|

	void DiscardOffsets();				// sets u[i], v[i], w[i] to zero
};

#endif /* MEMBRANE_INCL*/
