#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

#include"../header/membrane.h"




int main()
{
	Task test(TAU, H, DEF_DELTA, N_CELLS);

	test.SetFextPt(N_CELLS/2, N_CELLS/2, 0, 0, F_EXT);
//	test.SetFextPt(N_CELLS/2 - 1, N_CELLS/2, 0, 0, F_EXT);
//	test.SetFextPt(N_CELLS/2, N_CELLS/2 - 1, 0, 0, F_EXT);
//	test.SetFextPt(N_CELLS/2 - 1, N_CELLS/2 - 1, 0, 0, F_EXT);

	vtkSmartPointer<vtkStructuredGrid> sg = test.vtkSGrid();		
	vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

	char path[50];
	int i = 0;
	for(i = 0; i < FRAMES; i++)
	{
		for(int j = 0; j < ITER_PER_FRAME; j++)
			test.Iteration();
		sg = test.vtkSGrid();		

		sprintf(path, "output/f%d.vts", i);
		writer->SetFileName(path);
		writer->SetInputData(sg);
		writer->Write();
	}

	test.SetFextPt(N_CELLS/2, N_CELLS/2, 0, 0, 0);
	for(; i < 4*FRAMES; i++)
	{
		for(int j = 0; j < ITER_PER_FRAME; j++)
			test.Iteration();
		sg = test.vtkSGrid();		

		sprintf(path, "output/f%d.vts", i);
		writer->SetFileName(path);
		writer->SetInputData(sg);
		writer->Write();
	}
	
	return 0;
}
