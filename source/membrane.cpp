#include<stdio.h>
#include<strings.h>
#include<assert.h>
#include<math.h>

#include"../header/membrane.h"

const double ZERO_TOL = 1e-9;




//-------------------------------------------------------------------------------
//--------------@Point-----------------------------------------------------------
//-------------------------------------------------------------------------------

Point::Point(double x_0, double y_0, PointType pt)
{
	_x_0 = x_0; 	
	_y_0 = y_0; 	
	_type = pt;
	rho = 0; A = 0; B = 0; G = 0;

	r.zero();
	v.zero();
	dv.zero();
	F_ext.zero();
	sigma_xx = 0;
	sigma_xy = 0;
	sigma_yy = 0;
	
}

Point::Point(Point& orig)
{
	_x_0 = orig._x_0; 	
	_y_0 = orig._y_0; 	
	A = orig.A; 	
	B = orig.B; 	
	_y_0 = orig._y_0; 	
	rho = orig.rho; 	
	_type = orig._type;
	
	r.zero();
	v.zero();
	dv.zero();
	F_ext.zero();
	sigma_xx = 0;
	sigma_xy = 0;
	sigma_yy = 0;
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

	r.zero();
	v.zero();
	dv.zero();
	F_ext.zero();
	sigma_xx = 0;
	sigma_xy = 0;
	sigma_yy = 0;
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

	this->rho = rho;
	A = E/(1.0 - mu*mu);
	B = mu*A;
	G = (A - B)/2.0;
}

void Point::ToBorder()
{
	_type = BORDER;
}

bool Point::IsCorrect()
{
	double mu = B/A;
	if((rho > 0)&&(mu > -1.0)&&(mu <= 0.5)&&(A > 0)&&(fabs((A - B)/2.0 - G) <= ZERO_TOL))
		return true;
	return false;
}


//-------------------------------------------------------------------------------
//---------------@Grid2D----------------------------------------------------------
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
			(*this)[i][j].r.zero();
			(*this)[i][j].v.zero();
			(*this)[i][j].dv.zero();
//			(*this)[i][j].F_ext.zero();
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

vtkSmartPointer<vtkStructuredGrid> Task::vtkSGrid()
{
	vtkSmartPointer<vtkStructuredGrid> structuredGrid =
		    vtkSmartPointer<vtkStructuredGrid>::New();
 
	vtkSmartPointer<vtkPoints> points =
	    vtkSmartPointer<vtkPoints>::New();

	unsigned _x_nodes = _grid.x_nodes();
	unsigned _y_nodes = _grid.y_nodes();

	for(int i = 0; i < _y_nodes; i++)
	{
		for(int j = 0; j < _x_nodes; j++)
		{
			double x = _grid[i][j].x_0()
					+ _grid[i][j].r.u; 
			double y = _grid[i][j].y_0()
					+ _grid[i][j].r.v; 
			double z = _grid[i][j].r.w; 
			points->InsertNextPoint(x, y, z);
		}
	}
	
	structuredGrid->SetDimensions(_x_nodes, _y_nodes, 1);
	structuredGrid->SetPoints(points);
	return structuredGrid;
}


//-------------------------------------------------------------------------------
//---------------@Task-----------------------------------------------------------
//-------------------------------------------------------------------------------

Task::Task(double tau, double h, double delta,  unsigned cells) :
	_grid(cells, cells), _tau(tau), _h(h), _delta(delta)
{
	_grid.SetRect(_h, _h);
	_grid.OuterBorder();

	_grid.SetMaterialUniform(DEF_RHO, DEF_E, DEF_MU);
	// will be removed when task input from file is done
}


void Task::CountSigmas()						// assuming free border (sigma = 0 on outer border)
{									// counting only for inner parts
	unsigned X = _grid.x_nodes() - 1;
	unsigned Y = _grid.y_nodes() - 1;
	for(int i = 1; i < Y; i++)
	{
		for(int j = 1; j < X; j++)
		{
			vec3 delta_x = (_grid[i][j + 1].r - _grid[i][j - 1].r);
			vec3 delta_y = (_grid[i + 1][j].r - _grid[i - 1][j].r);

			double w_x = delta_x.w/2.0/_h;
			double w_y = delta_y.w/2.0/_h;

			vec3 r_x(1, 0, w_x);
			vec3 r_y(0, 1, w_y);
			
			double Pi_xx = vec3::Projection(delta_x, r_x);
			double Pi_xy = vec3::Projection(delta_x, r_y);
			double Pi_yx = vec3::Projection(delta_y, r_x);
			double Pi_yy = vec3::Projection(delta_y, r_y);
			
			double eps_xx = Pi_xx/2.0/_h;
			double eps_yy = Pi_yy/2.0/_h;
			double eps_xy = (Pi_xy + Pi_yx)/4.0/_h;

			_grid[i][j].sigma_xx = _grid[i][j].A*eps_xx + _grid[i][j].B*eps_yy;	
			_grid[i][j].sigma_yy = _grid[i][j].A*eps_yy + _grid[i][j].B*eps_xx;	
			_grid[i][j].sigma_xy = _grid[i][j].G*eps_xy; 
		}
	}
}	

void Task::DvInner()
{
	unsigned X = _grid.x_nodes() - 1;
	unsigned Y = _grid.y_nodes() - 1;
	for(int i = 1; i < Y; i++)					// inner part 
	{
		for(int j = 1; j < X; j++)
		{
			double vx_t = 0.0;
			vx_t += (_grid[i][j + 1].sigma_xx - _grid[i][j - 1].sigma_xx)/_h; 	
			vx_t += (_grid[i + 1][j].sigma_xy - _grid[i - 1][j].sigma_xy)/_h; 	

			double vy_t = 0.0;
			vy_t += (_grid[i + 1][j].sigma_yy - _grid[i - 1][j].sigma_yy)/_h; 	
			vy_t += (_grid[i][j + 1].sigma_xy - _grid[i][j - 1].sigma_xy)/_h; 	
			
			
			double w_x = (_grid[i][j + 1].r.w - _grid[i][j - 1].r.w)/2.0/_h;
			double w_y = (_grid[i + 1][j].r.w - _grid[i - 1][j].r.w)/2.0/_h;

			double norm_x = sqrt(1 + w_x*w_x);
			double norm_y = sqrt(1 + w_y*w_y);
			
			_grid[i][j].dv.u = vx_t/norm_x;
			_grid[i][j].dv.v = vy_t/norm_y;
			_grid[i][j].dv.w = _grid[i][j].dv.u*w_x + _grid[i][j].dv.v*w_y;
			
			_grid[i][j].dv += (1.0/_h/_h/_delta)*_grid[i][j].F_ext;

			_grid[i][j].dv *= _tau/_grid[i][j].rho;
		}
	}

}

void Task::DvTopBot()
{
	unsigned X = _grid.x_nodes() - 1;
	unsigned Y = _grid.y_nodes() - 1;
	for(int i = 1; i < X; i++)					// top and bottom borders
	{
		// bottom border (y = 0, x = i)
		double vx_t = 0.0;
		vx_t += (_grid[0][i + 1].sigma_xx - _grid[0][i - 1].sigma_xx)/_h; 	
		vx_t += (_grid[1][i].sigma_xy)/_h; 	

		double vy_t = 0.0;
		vy_t += (_grid[1][i].sigma_yy)/_h; 	
		vy_t += (_grid[0][i + 1].sigma_xy - _grid[0][i - 1].sigma_xy)/_h; 	

		double w_x = (_grid[0][i + 1].r.w - _grid[0][i - 1].r.w)/2.0/_h;
		double w_y = (_grid[1][i].r.w - _grid[0][i].r.w)/_h;
		double norm_x = sqrt(1 + w_x*w_x);
		double norm_y = sqrt(1 + w_y*w_y);

		_grid[0][i].dv.u = vx_t/norm_x;
		_grid[0][i].dv.v = vy_t/norm_y;
		_grid[0][i].dv.w = _grid[0][i].dv.u*w_x + _grid[0][i].dv.v*w_y;
		
		_grid[0][i].dv += (1.0/_h/_h/_delta)*_grid[0][i].F_ext;

		_grid[0][i].dv *= _tau/_grid[0][i].rho;
		//top border ( y = Y, x = i)
		vx_t = 0.0;
		vx_t += (_grid[Y][i + 1].sigma_xx - _grid[Y][i - 1].sigma_xx)/_h; 	
		vx_t -= (_grid[Y - 1][i].sigma_xy)/_h; 	

		vy_t = 0.0;
		vy_t -= (_grid[Y - 1][i].sigma_yy)/_h; 	
		vy_t += (_grid[Y][i + 1].sigma_xy - _grid[0][i - 1].sigma_xy)/_h; 	

		w_x = (_grid[Y][i + 1].r.w - _grid[Y][i - 1].r.w)/2.0/_h;
		w_y = (_grid[Y][i].r.w - _grid[Y - 1][i].r.w)/_h;
		norm_x = sqrt(1 + w_x*w_x);
		norm_y = sqrt(1 + w_y*w_y);

		_grid[Y][i].dv.u = vx_t/norm_x;
		_grid[Y][i].dv.v = vy_t/norm_y;
		_grid[Y][i].dv.w = _grid[Y][i].dv.u*w_x + _grid[Y][i].dv.v*w_y;
		
		_grid[Y][i].dv += (1.0/_h/_h/_delta)*_grid[Y][i].F_ext;

		_grid[Y][i].dv *= _tau/_grid[Y][i].rho;
	} 
}

void Task::DvRightLeft()
{
	unsigned X = _grid.x_nodes() - 1;
	unsigned Y = _grid.y_nodes() - 1;

	for(int i = 1; i < Y; i++)					// right and left borders
	{
		//left border (y = i, x = 0)
		double vx_t = 0.0;
		vx_t += (_grid[i][1].sigma_xx)/_h; 	
		vx_t += (_grid[i + 1][0].sigma_xy - _grid[i - 1][0].sigma_xy)/_h; 	

		double vy_t = 0.0;
		vy_t += (_grid[i + 1][0].sigma_yy - _grid[i - 1][0].sigma_yy)/_h; 	
		vy_t += (_grid[i][1].sigma_xy)/_h; 	
		
		double w_x = (_grid[i][1].r.w - _grid[i][0].r.w)/_h;
		double w_y = (_grid[i + 1][0].r.w - _grid[i - 1][0].r.w)/2.0/_h;

		double norm_x = sqrt(1 + w_x*w_x);
		double norm_y = sqrt(1 + w_y*w_y);
		
		_grid[i][0].dv.u = vx_t/norm_x;
		_grid[i][0].dv.v = vy_t/norm_y;
		_grid[i][0].dv.w = _grid[i][0].dv.u*w_x + _grid[i][0].dv.v*w_y;
		
		_grid[i][0].dv += (1.0/_h/_h/_delta)*_grid[i][0].F_ext;
		_grid[i][0].dv *= _tau/_grid[i][0].rho;

		// right border (y = i, x = X)
		vx_t = 0.0;
		vx_t -=  (_grid[i][X - 1].sigma_xx)/_h; 	
		vx_t += (_grid[i + 1][X].sigma_xy - _grid[i - 1][X].sigma_xy)/_h; 	

		vy_t = 0.0;
		vy_t += (_grid[i + 1][X].sigma_yy - _grid[i - 1][X].sigma_yy)/_h; 	
		vy_t -= (_grid[i][X - 1].sigma_xy)/_h; 	
		
		w_x = (_grid[i][X].r.w - _grid[i][X - 1].r.w)/_h;
		w_y = (_grid[i + 1][X].r.w - _grid[i - 1][X].r.w)/2.0/_h;

		norm_x = sqrt(1 + w_x*w_x);
		norm_y = sqrt(1 + w_y*w_y);
		
		_grid[i][X].dv.u = vx_t/norm_x;
		_grid[i][X].dv.v = vy_t/norm_y;
		_grid[i][X].dv.w = _grid[i][X].dv.u*w_x + _grid[i][X].dv.v*w_y;
		
		_grid[i][X].dv += (1.0/_h/_h/_delta)*_grid[i][X].F_ext;
		_grid[i][X].dv *= _tau/_grid[i][X].rho;

	} 
}

void Task::DvCorners()
{
	unsigned X = _grid.x_nodes() - 1;
	unsigned Y = _grid.y_nodes() - 1;
	// left bottom corner (y = 0, x = 0)
	double vx_t = (_grid[1][0].sigma_xy + _grid[0][1].sigma_xx)/_h;
	double vy_t = (_grid[1][0].sigma_yy + _grid[0][1].sigma_xy)/_h;

	double w_x = (_grid[0][1].r.w - _grid[0][0].r.w)/_h;
	double w_y = (_grid[1][0].r.w - _grid[0][0].r.w)/_h;
	
	double norm_x = sqrt(1 + w_x*w_x);
	double norm_y = sqrt(1 + w_y*w_y);
	_grid[0][0].dv.u = vx_t/norm_x;
	_grid[0][0].dv.v = vy_t/norm_y;
	_grid[0][0].dv.w = _grid[0][0].dv.u*w_x + _grid[0][0].dv.v*w_y;
	_grid[0][0].dv += (1.0/_h/_h/_delta)*_grid[0][0].F_ext;
	_grid[0][0].dv *= _tau/_grid[0][0].rho;
	// right bottom corner (y = 0, x = X)
	vx_t = (_grid[1][X].sigma_xy - _grid[0][X - 1].sigma_xx)/_h;
	vy_t = (_grid[1][X].sigma_yy - _grid[0][X - 1].sigma_xy)/_h;

	w_x = (- _grid[0][X - 1].r.w + _grid[0][X].r.w)/_h;
	w_y = (_grid[1][X].r.w - _grid[0][X].r.w)/_h;
	
	norm_x = sqrt(1 + w_x*w_x);
	norm_y = sqrt(1 + w_y*w_y);
	_grid[0][X].dv.u = vx_t/norm_x;
	_grid[0][X].dv.v = vy_t/norm_y;
	_grid[0][X].dv.w = _grid[0][X].dv.u*w_x + _grid[0][X].dv.v*w_y;
	_grid[0][X].dv += (1.0/_h/_h/_delta)*_grid[0][X].F_ext;
	_grid[0][X].dv *= _tau/_grid[0][X].rho;
	// right top corner (y = Y, x = X)
	vx_t = -(_grid[Y - 1][X].sigma_xy + _grid[Y][X - 1].sigma_xx)/_h;
	vy_t = -(_grid[Y - 1][X].sigma_yy + _grid[Y][X - 1].sigma_xy)/_h;

	w_x = ( _grid[Y][X].r.w - _grid[Y][X - 1].r.w)/_h;
	w_y = (_grid[Y][X].r.w - _grid[Y - 1][X].r.w)/_h;
	
	norm_x = sqrt(1 + w_x*w_x);
	norm_y = sqrt(1 + w_y*w_y);
	_grid[Y][X].dv.u = vx_t/norm_x;
	_grid[Y][X].dv.v = vy_t/norm_y;
	_grid[Y][X].dv.w = _grid[Y][X].dv.u*w_x + _grid[Y][X].dv.v*w_y;
	_grid[Y][X].dv += (1.0/_h/_h/_delta)*_grid[Y][X].F_ext;
	_grid[Y][X].dv *= _tau/_grid[Y][X].rho;
	// left top corner (y = Y, x = 0)
	vx_t = ( -_grid[Y - 1][0].sigma_xy + _grid[Y][1].sigma_xx)/_h;
	vy_t = ( -_grid[Y - 1][0].sigma_yy + _grid[Y][1].sigma_xy)/_h;

	w_x = ( -_grid[Y][0].r.w + _grid[Y][1].r.w)/_h;
	w_y = (_grid[Y][0].r.w - _grid[Y - 1][0].r.w)/_h;
	
	norm_x = sqrt(1 + w_x*w_x);
	norm_y = sqrt(1 + w_y*w_y);
	_grid[Y][0].dv.u = vx_t/norm_x;
	_grid[Y][0].dv.v = vy_t/norm_y;
	_grid[Y][0].dv.w = _grid[Y][0].dv.u*w_x + _grid[Y][0].dv.v*w_y;
	_grid[Y][0].dv += (1.0/_h/_h/_delta)*_grid[Y][0].F_ext;
	_grid[Y][0].dv *= _tau/_grid[Y][0].rho;
}

void Task::Iteration()
{
	this->CountSigmas();
	this->DvInner();
	this->DvTopBot();
	this->DvRightLeft();
	this->DvCorners();

	unsigned X = _grid.x_nodes();
	unsigned Y = _grid.y_nodes();
	for(int i = 0; i < Y; i++)
	{
		for(int j = 0; j < X; j++)
		{
			_grid[i][j].v += _grid[i][j].dv;
			_grid[i][j].r += _tau*_grid[i][j].v;
		}
	}
}


void Task::SetFextPt(unsigned x, unsigned y, double F_u, double F_v, double F_w)
{
	if(x >= _grid.x_nodes() || y >= _grid.y_nodes())
	{
		printf("Error: point out of range!\n");
		assert(0);
	}
	_grid[y][x].F_ext.u = F_u;
	_grid[y][x].F_ext.v = F_v;
	_grid[y][x].F_ext.w = F_w;
}
			
















