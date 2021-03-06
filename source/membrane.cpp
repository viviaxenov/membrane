#include<stdio.h>
#include<strings.h>
#include<sstream>
#include<assert.h>
#include<math.h>

#include<sys/types.h>
#include<sys/stat.h>

#include"../header/membrane.h"
#include"../header/constants.h"

const double ZERO_TOL = 1e-9;




//-------------------------------------------------------------------------------
//--------------@Point-----------------------------------------------------------
//-------------------------------------------------------------------------------

Point::Point(double x_0, double y_0, PointType pt)
{
	_x_0 = x_0; 	
	_y_0 = y_0; 	
	type = pt;
	rho = 0; A = 0; B = 0; G = 0;

	r.zero();
	v.zero();
	a.zero();
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
	type = orig.type;
	
	r.zero();
	v.zero();
	a.zero();
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
			_x_0, _y_0, type);
	string res(s);
	delete s;
	return res;
}

void Point::Set(double x_0, double y_0, PointType pt)
{
	_x_0 = x_0; 	
	_y_0 = y_0; 	
	type = pt;

	r.zero();
	v.zero();
	a.zero();
	F_ext.zero();
	sigma_xx = 0;
	sigma_xy = 0;
	sigma_yy = 0;
}

void Point::SetMaterial(double rho, double E, double mu, double delta)
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
	this->delta = delta;
	A = E/(1.0 - mu*mu);
	B = mu*A;
	G = (A - B)/2.0;
}


bool Point::IsCorrect()
{
	if(type == INACTIVE)
		return true;
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



void Grid2D::SetRect(double h_x , double h_y)
{
	double x = 0, y = 0;	
	for(unsigned i = 0; i < _y_nodes; i++)
	{
		x = 0;
		for(unsigned j = 0; j < _x_nodes; j++)
		{			
			(*this)[i][j].Set(x, y);
			x += h_x;	
		}
		y += h_y;
	}
	_h_x = h_x;
	_h_y = h_y;
}


void Grid2D::SetMaterialUniform(double rho, double E, double mu, double delta)
{
	for(unsigned i = 0; i < _y_nodes; i++)
	{
		for(unsigned j = 0; j < _x_nodes; j++)
		{			
			(*this)[i][j].SetMaterial(rho, E, mu, delta);
		}
	}
}


void Grid2D::OuterBorder()
{
	for(int i = 1; i < _y_nodes - 1; i++)
	{
		(*this)[i][0].type = L_BORDER;		
		(*this)[i][_x_nodes - 1].type = R_BORDER;
	}
	for(int i = 1; i < _x_nodes - 1; i++)
	{
		(*this)[0][i].type = B_BORDER;		
		(*this)[_y_nodes - 1][i].type = T_BORDER;		
	}
	(*this)[0][0].type = BL_CORNER;		
	(*this)[0][_x_nodes - 1].type = BR_CORNER;		
	(*this)[_y_nodes - 1][0].type = TL_CORNER;		
	(*this)[_y_nodes - 1][_x_nodes - 1].type = TR_CORNER;		
}

void Grid2D::DiscardOffsets()
{
	for(unsigned i = 0; i < _y_nodes; i++)
	{
		for(unsigned j = 0; j < _x_nodes; j++)
		{			
			(*this)[i][j].r.zero();
			(*this)[i][j].v.zero();
			(*this)[i][j].a.zero();
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
	// dumping coordinates

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

	// dumping velocity vectors

	vtkSmartPointer<vtkDoubleArray> vel = 
		vtkSmartPointer<vtkDoubleArray>::New();
	vel->SetName("Velocity");
	vel->SetNumberOfComponents(3);
	vel->SetNumberOfTuples(_x_nodes*_y_nodes);
	for(int i = 0; i < _y_nodes; i++)
	{
		for(int j = 0; j < _x_nodes; j++)
		{
			vel->SetTuple3(_x_nodes*i + j, 
				_grid[i][j].v.u,
				_grid[i][j].v.v,
				_grid[i][j].v.w);
		}
	}
	structuredGrid->GetPointData()->AddArray(vel);
	
	// dumping elastic force data
	vtkSmartPointer<vtkDoubleArray> sigma_xx = 
		vtkSmartPointer<vtkDoubleArray>::New(),
		sigma_yy = vtkSmartPointer<vtkDoubleArray>::New(),
		sigma_xy = vtkSmartPointer<vtkDoubleArray>::New();
	sigma_xx->SetName("sigma_xx");
	sigma_xy->SetName("sigma_xy");
	sigma_yy->SetName("sigma_yy");
	sigma_xx->SetNumberOfComponents(1);
	sigma_xy->SetNumberOfComponents(1);
	sigma_yy->SetNumberOfComponents(1);
	sigma_xx->SetNumberOfTuples(_x_nodes*_y_nodes);
	sigma_yy->SetNumberOfTuples(_x_nodes*_y_nodes);
	sigma_xy->SetNumberOfTuples(_x_nodes*_y_nodes);
	for(int i = 0; i < _y_nodes; i++)
	{
		for(int j = 0; j < _x_nodes; j++)
		{
			sigma_xx->SetValue(_x_nodes*i + j, 
						_grid[i][j].sigma_xx);	
			sigma_xy->SetValue(_x_nodes*i + j, 
						_grid[i][j].sigma_xy);	
			sigma_yy->SetValue(_x_nodes*i + j, 
						_grid[i][j].sigma_yy);	
		}
	}
	 
	structuredGrid->GetPointData()->AddArray(sigma_xx);
	structuredGrid->GetPointData()->AddArray(sigma_xy);
	structuredGrid->GetPointData()->AddArray(sigma_yy);

	return structuredGrid;
}

vec3 Grid2D::delta_x(unsigned x, unsigned y)				// NOTE: count distance between cells regerdless of 
									// _h_x, _h_y so fit for arbitrary grid
{
	if(x >= _x_nodes || y >= _y_nodes)
	{
		assert(!"ERROR: Index out of range!");
	}
	
	double dx;
	vec3 delta;
	switch((*this)[y][x].type)
	{
		case GRID:
		case T_BORDER:
		case B_BORDER:
			dx = (*this)[y][x + 1].x_0() - (*this)[y][x - 1].x_0();
			delta = (*this)[y][x + 1].r - (*this)[y][x - 1].r; 
			break;
		case L_BORDER:
		case TL_CORNER:
		case BL_CORNER:
			dx = (*this)[y][x + 1].x_0() - (*this)[y][x].x_0(); 
			delta = (*this)[y][x + 1].r - (*this)[y][x].r; 
			break;
		case R_BORDER:
		case TR_CORNER:
		case BR_CORNER:
			dx = (*this)[y][x].x_0() - (*this)[y][x - 1].x_0(); 
			delta = (*this)[y][x].r - (*this)[y][x - 1].r; 
			break;
		case INACTIVE: 
			assert(!"ERROR: Trying to calculate derivative for inactive cell");
	}
	return (1.0/dx)*delta;
}

vec3 Grid2D::delta_y(unsigned x, unsigned y)
{
	if(x >= _x_nodes || y >= _y_nodes)
	{
		assert(!"ERROR: Index out of range!");
	}
	
	double dy;
	vec3 delta;
	switch((*this)[y][x].type)
	{
		case GRID:
		case L_BORDER:
		case R_BORDER:
			dy = (*this)[y + 1][x].y_0() - (*this)[y - 1][x].y_0();
			delta = (*this)[y + 1][x].r - (*this)[y - 1][x].r; 
			break;
		case B_BORDER:
		case BL_CORNER:
		case BR_CORNER:
			dy = (*this)[y + 1][x].y_0() - (*this)[y][x].y_0(); 
			delta = (*this)[y + 1][x].r - (*this)[y][x].r; 
			break;
		case T_BORDER:
		case TR_CORNER:
		case TL_CORNER:
			dy = (*this)[y][x].y_0() - (*this)[y - 1][x].y_0(); 
			delta = (*this)[y][x].r - (*this)[y - 1][x].r; 
			break;
		case INACTIVE: 
			assert(!"ERROR: Trying to calculate derivative for inactive cell");
	}
	return (1.0/dy)*delta;
}

void Grid2D::CountSigmas()							// assuming free border 
{									
	unsigned X = (*this).x_nodes();
	unsigned Y = (*this).y_nodes();
	for(int i = 0; i < Y; i++)
	{
		for(int j = 0; j < X; j++)
		{
			vec3 delta_x = (*this).delta_x(j, i);
			vec3 delta_y = (*this).delta_y(j, i);

			double w_x = delta_x.w;
			double w_y = delta_y.w;
			(*this)[i][j].w_x = w_x;				// saving for further use
			(*this)[i][j].w_y = w_y;

			vec3 r_x(1, 0, w_x);
			vec3 r_y(0, 1, w_y);
			
			double Pi_xx = vec3::Projection(delta_x, r_x);
			double Pi_xy = vec3::Projection(delta_x, r_y);
			double Pi_yx = vec3::Projection(delta_y, r_x);
			double Pi_yy = vec3::Projection(delta_y, r_y);
			
			double eps_xx = Pi_xx;
			double eps_yy = Pi_yy;
			double eps_xy = (Pi_xy + Pi_yx)/2.0;

			(*this)[i][j].sigma_xx = (*this)[i][j].A*eps_xx + (*this)[i][j].B*eps_yy;	
			(*this)[i][j].sigma_yy = (*this)[i][j].A*eps_yy + (*this)[i][j].B*eps_xx;	
			(*this)[i][j].sigma_xy = (*this)[i][j].G*eps_xy; 
		}
	}
}	


void Grid2D::CountAcceleration()
{
	unsigned X = (*this).x_nodes() - 1;
	unsigned Y = (*this).y_nodes() - 1;
	for(int i = 0; i < Y; i++)					
	{
		for(int j = 0; j < X; j++)
		{

			double vx_t = 0.0;
			double vy_t = 0.0;
			PointType pt = (*this)[i][j].type;
			if(pt == INACTIVE)
				continue;
			// TODO: count dx and dy somehow for case of irregular grid (like i did it for sigmas)
			if(pt != TL_CORNER && pt != L_BORDER && pt != BL_CORNER)// if has left neighbor
			{
				vx_t -= (*this)[i][j - 1].sigma_xx/_h_x;
				vy_t -=	(*this)[i][j - 1].sigma_xy/_h_x;
			}
			if(pt != TR_CORNER && pt != R_BORDER && pt != BR_CORNER)// if has right neighbor
			{
				vx_t += (*this)[i][j + 1].sigma_xx/_h_x;
				vy_t +=	(*this)[i][j + 1].sigma_xy/_h_x;
			}
			if(pt != TR_CORNER && pt != T_BORDER && pt != TL_CORNER)// if has top neighbor
			{
				vx_t += (*this)[i + 1][j].sigma_xy/_h_y;
				vy_t +=	(*this)[i + 1][j].sigma_yy/_h_y;
			}
			if(pt != BR_CORNER && pt != B_BORDER && pt != BL_CORNER)// if has bottom neighbor
			{
				vx_t -= (*this)[i - 1][j].sigma_xy/_h_y;
				vy_t -=	(*this)[i - 1][j].sigma_yy/_h_y;
			}


	
			double w_x = (*this)[i][j].w_x;
			double w_y = (*this)[i][j].w_y;
			
			double norm_x = sqrt(1 + w_x*w_x);
			double norm_y = sqrt(1 + w_y*w_y);
			
			(*this)[i][j].a.u = vx_t/norm_x;
			(*this)[i][j].a.v = vy_t/norm_y;
			(*this)[i][j].a.w =// vx_t*w_x/norm_x + vy_t*w_y/norm_y;
					 (*this)[i][j].a.u*w_x + (*this)[i][j].a.v*w_y;
			
			double delta = (*this)[i][j].delta;
			(*this)[i][j].a += (1.0/_h_x/_h_y/delta)*(*this)[i][j].F_ext;

			(*this)[i][j].a *= 1.0/(*this)[i][j].rho;
		}
	}
}



//-------------------------------------------------------------------------------
//---------------@Task-----------------------------------------------------------
//-------------------------------------------------------------------------------

Task::Task(double tau, double h, unsigned cells) :
	_grid(cells, cells), _tau(tau)
{
	_grid.SetRect(h, h);
	// assume membrane to be square; will make other geometry
	// when file input is done
	_grid.OuterBorder();

	_grid.SetMaterialUniform(DEF_RHO, DEF_E, DEF_MU, DEF_DELTA);
	// will be removed when task input from file is done
	_output_dir = "";						

	frames = 0;
	iter_per_frame = 0;
}

int Task::SetOutputDir(string of)
{
	struct stat sb;
	if(stat(of.c_str(), &sb) == -1 || ((sb.st_mode & S_IFMT) != S_IFDIR))
	{
		if(mkdir(of.c_str(), S_IRUSR | S_IWUSR | S_IXUSR) != 0)		// not sure of those flags
			return -1;
	}
	
	_output_dir = of;
	return 0;
}

void Task::Iteration()
{
	_grid.CountSigmas();
	_grid.CountAcceleration();

	unsigned X = _grid.x_nodes();
	unsigned Y = _grid.y_nodes();
	for(int i = 0; i < Y; i++)
	{
		for(int j = 0; j < X; j++)
		{
			_grid[i][j].v += _tau*_grid[i][j].a;
			_grid[i][j].r += _tau*_grid[i][j].v;
		}
	}
}

int Task::Execute()
{
	if(frames == 0 || iter_per_frame == 0 || _output_dir == "")
		return -1;

	vtkSmartPointer<vtkStructuredGrid> sg = (*this).vtkSGrid();		
	vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

	for(int i = 0; i <= frames; i++)
	{
		for(int j = 0; (j < iter_per_frame) && (i != 0); j++)
			(*this).Iteration();

		sg = (*this).vtkSGrid();		
		
		std::ostringstream path;
		path << _output_dir << "/f" << i << ".vts";

		writer->SetFileName(path.str().c_str());
		writer->SetInputData(sg);
		writer->Write();
	}

	return 0;
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
			


void Task::SetOffPt(unsigned x, unsigned y, double u, double v, double w)
{
	
	if(x >= _grid.x_nodes() || y >= _grid.y_nodes())
	{
		printf("Error: point out of range!\n");
		assert(0);
	}
	_grid[y][x].r.u = u;
	_grid[y][x].r.v = v;
	_grid[y][x].r.w = w;
}
			














