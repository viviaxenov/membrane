#include<stdio.h>
#include<strings.h>
#include<assert.h>


#include"../header/membrane.h"


//-------------------------------------------------------------------------------
//---------------Point-----------------------------------------------------------
//-------------------------------------------------------------------------------

Point::Point(double x_0, double y_0, PointType pt)
{
	_x_0 = x_0; 	
	_y_0 = y_0; 	
	_type = pt;
	bzero(u, sizeof(u));
	bzero(v, sizeof(v));
	bzero(w, sizeof(w));
}

Point::Point(Point& orig)
{
	_x_0 = orig._x_0; 	
	_y_0 = orig._y_0; 	
	_type = orig._type;
	bzero(u, sizeof(u));
	bzero(v, sizeof(v));
	bzero(w, sizeof(w));
}

string Point::ToString()
{
	char *s = new char[50];
	sprintf(s,	"x = %lg\n"
			"y = %lg\n"
			"Type %d\n",
			_x_0, _y_0, _type);
	string res(s);
	delete s;
	return res;
}

void Point::Set(double x_0, double y_0, PointType pt)
{
	_x_0 = x_0; 	
	_y_0 = y_0; 	
	_type = pt;
	bzero(u, sizeof(u));
	bzero(v, sizeof(v));
	bzero(w, sizeof(w));
}
void Point::ToBorder()
{
	_type = BORDER;
}

//-------------------------------------------------------------------------------
//---------------Grid2D----------------------------------------------------------
//-------------------------------------------------------------------------------

Grid2D::Grid2D(unsigned x_nodes, unsigned y_nodes)
{

	if(x_nodes == 0 || y_nodes == 0)
	{
		printf("ERROR: zero dimension!\n");
		assert(0);
	}
/*	if(x_nodes == 0 || y_nodes == 0)
	{
		_x_nodes = 0;
		_y_nodes = 0;
		_grid = NULL;
		return;
	}
*/
	_x_nodes = x_nodes;
	_y_nodes = y_nodes;
	
	_grid = new Point[_x_nodes*_y_nodes];
}



void Grid2D::SetRect(double dx , double dy)
{
	double x = 0, y = 0;	
	for(unsigned i = 0; i < _y_nodes; i++)
	{
		x = 0;
		for(unsigned j = 0; j < _x_nodes; j++)
		{			
			(*this)[i][j].Set(x, y);
			x += dx;	
		}
		y += dy;
	}
}

void Grid2D::OuterBorder()
{
	for(int i = 0; i < _y_nodes; i++)
	{
		(*this)[i][0].ToBorder();		
		(*this)[i][_x_nodes - 1].ToBorder();		
	}
	for(int i = 0; i < _x_nodes; i++)
	{
		(*this)[0][i].ToBorder();		
		(*this)[_y_nodes - 1][i].ToBorder();		
	}
}

void Grid2D::DiscardOffsets()
{
	for(unsigned i = 0; i < _y_nodes; i++)
	{
		for(unsigned j = 0; j < _x_nodes; j++)
		{			
			bzero((*this)[i][j].u, sizeof((*this)[i][j].u));
			bzero((*this)[i][j].v, sizeof((*this)[i][j].v));
			bzero((*this)[i][j].w, sizeof((*this)[i][j].w));
		}
	}
}
