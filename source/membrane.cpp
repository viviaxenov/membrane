#include<stdio.h>
#include<strings.h>
#include<assert.h>
#include<math.h>

#include"../header/membrane.h"

const double ZERO_TOL = 1e-9;




//-------------------------------------------------------------------------------
//---------------Point-----------------------------------------------------------
//-------------------------------------------------------------------------------

Point::Point(double x_0, double y_0, PointType pt)
{
	_x_0 = x_0; 	
	_y_0 = y_0; 	
	_type = pt;
	_rho = 0; _A = 0; _B = 0; _G = 0;
	u = 0; v_u = 0; dv_u = 0;
	v = 0; v_v = 0; dv_v = 0;
	w = 0; v_w = 0; dv_w = 0;
	
}

Point::Point(Point& orig)
{
	_x_0 = orig._x_0; 	
	_y_0 = orig._y_0; 	
	_A = orig._A; 	
	_B = orig._B; 	
	_y_0 = orig._y_0; 	
	_rho = orig._rho; 	
	_type = orig._type;
	
	u = 0; v_u = 0; dv_u = 0;
	v = 0; v_v = 0; dv_v = 0;
	w = 0; v_w = 0; dv_w = 0;
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
	u = 0; v_u = 0; dv_u = 0;
	v = 0; v_v = 0; dv_v = 0;
	w = 0; v_w = 0; dv_w = 0;
}

void Point::SetMaterial(double rho, double E, double mu)
{
	if(mu > 0.5 || mu <= - 1.0)
	{
		puts("ERROR: wrong Poisson's coef\n");
		assert(0);
	}	
	if(E <= 0)
	{
		puts("ERROR: negative Young's modulus\n");
		assert(0);
	}	
	if(rho <= 0)
	{
		puts("ERROR: negative density\n");
		assert(0);
	}	

	_rho = rho;
	_A = E/(1.0 - mu*mu);
	_B = mu*_A;
	_G = (_A - _B)/2.0;
}

void Point::ToBorder()
{
	_type = BORDER;
}

bool Point::IsCorrect()
{
	double mu = _B/_A;
	if((_rho > 0)&&(mu > -1.0)&&(mu <= 0.5)&&(_A > 0)&&(fabs((_A - _B)/2.0 - _G) <= ZERO_TOL))
		return true;
	return false;
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


void Grid2D::SetMaterialUniform(double rho, double E, double mu)
{
	for(unsigned i = 0; i < _y_nodes; i++)
	{
		for(unsigned j = 0; j < _x_nodes; j++)
		{			
			(*this)[i][j].SetMaterial(rho, E, mu);
		}
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
			(*this)[i][j].u = 0; (*this)[i][j].v_u = 0; (*this)[i][j].dv_u = 0;
			(*this)[i][j].v = 0; (*this)[i][j].v_v = 0; (*this)[i][j].dv_v = 0;
			(*this)[i][j].w = 0; (*this)[i][j].v_w = 0; (*this)[i][j].dv_w = 0;
		}
	}
}

bool Grid2D::IsCorrect()
{
	if(_x_nodes == 0 || _y_nodes == 0 || _grid == NULL)
		return false;
	bool status = true;
	for(unsigned i = 0; status&&(i < _y_nodes); i++)	
	{
		for(unsigned j = 0; status&&(j < _x_nodes); j++)	
		{
			status = (*this)[i][j].IsCorrect(); 
		}
	}
	return status;
}	

//vtkSmartPointer<vtkStructuredGrid> Grid2D::vtkSGrid()
//{
//	vtkSmartPointer<vtkStructuredGrid> structuredGrid =
//		    vtkSmartPointer<vtkStructuredGrid>::New();
// 
//	vtkSmartPointer<vtkPoints> points =
//	    vtkSmartPointer<vtkPoints>::New();
//	for(int i = 0; i < _y_nodes; i++)
//	{
//		for(int j = 0; j < _x_nodes; j++)
//		{
//			double x = (*this)[i][j].x_0()
//					+ (*this)[i][j].u[1]; 
//			double y = (*this)[i][j].y_0()
//					+ (*this)[i][j].v[1]; 
//			double z = (*this)[i][j].w[1]; 
//			points->InsertNextPoint(x, y, z);
//		}
//	}
//	
//	structuredGrid->SetDimensions(_x_nodes, _y_nodes, 1);
//	structuredGrid->SetPoints(points);
//	return structuredGrid;
//}





