#include<string>

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

//	double _sigma_xx, _sigma_yy, _sigma_zz;		// stress tensor
//	double _tau_xy, _tau_xz, _tau_yz;
	

	Point(double x_0 = 0, double y_0 = 0, PointType pt = INACTIVE);
	Point(Point&);
	~Point(){};

	PointType type() {return _type;}
	double x_0() { return _x_0;}			
	double y_0() { return _y_0;}		

	void Set(double x_0, double y_0, PointType tp = GRID);
	void ToBorder();

	string ToString();
};


class Grid2D
{
private: 
	unsigned _x_nodes;
	unsigned _y_nodes;

	Point *_grid;
public: 
	Grid2D(unsigned x_nodes = 0, unsigned y_nodes = 0);
//	Grid2D(Grid2D&);
	~Grid2D(){};

	Point* operator [](int y)
	{
		return _grid + y*_x_nodes;
	}

	void SetRect(double dx , double dy);
	void OuterBorder();
	void DiscardOffsets();
};


