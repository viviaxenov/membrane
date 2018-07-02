
enum Point_Type 
{
	BORDER,
	GRID,
	INACTIVE;
}



class Point
{
private:
	double _x_0;					// initial positions
	double _y_0;		 

	double _u;					// offsets: 	u - Ox axis
	double _v;					//		v - Oy axis
	double _w;					//		w - Oz (vertical)

//	double _sigma_xx, _sigma_yy, _sigma_zz;		// stress tensor
//	double _tau_xy, _tau_xz, _tau_yz;
	
	PointType _type;

public: 

	Point(double x_0, double y_0, PointType pt);
	~Point(){};

	double x_0 {return _x_0;}
	double y_0 {return _y_0;}
	double z_0 {return _z_0;}

	double u {return _u;}
	double v {return _v;}
	double w {return _w;}
	
	PointType _type {return _type};
}


class Grid2D
{
private: 
	unsigned _x_nodes;
	unsigned _y_nodes;

	Point *_grid;
public: 
	Grid2D(unsigned x_nodes, unsigned y_nodes, double x_step, double y_step);
	~Grid2D();
	
	
}
Grid2D::Point* operator [](int y)
{
	return _grid + y*_x_nodes;
}


