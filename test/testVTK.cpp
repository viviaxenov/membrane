#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

#include"../header/membrane.h"




int main()
{
	Task test(TAU, H, N_CELLS);

	test.SetOffPt(N_CELLS/2, N_CELLS/2, 0, 0, 0.5*H);
	test.SetOffPt(N_CELLS/2 - 1, N_CELLS/2, 0, 0, 0.5*H);
	test.SetOffPt(N_CELLS/2, N_CELLS/2 - 1, 0, 0, 0.5*H);
	test.SetOffPt(N_CELLS/2 - 1, N_CELLS/2 - 1, 0, 0, 0.5*H);

	vtkSmartPointer<vtkStructuredGrid> sg = test.vtkSGrid();		
	vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

	char path[50];
	int i = 0;
	for(i = 1; i <= FRAMES; i++)
	{
		for(int j = 0; j < ITER_PER_FRAME; j++)
			test.Iteration();
		sg = test.vtkSGrid();		

		sprintf(path, "output/f%d.vts", i);
		writer->SetFileName(path);
		writer->SetInputData(sg);
		writer->Write();
	}

//	test.SetFextPt(N_CELLS/2, N_CELLS/2, 0, 0, 0);
//	for(; i < 4*FRAMES; i++)
//	{
//		for(int j = 0; j < ITER_PER_FRAME; j++)
//			test.Iteration();
//		sg = test.vtkSGrid();		
//
//		sprintf(path, "output/f%d.vts", i);
//		writer->SetFileName(path);
//		writer->SetInputData(sg);
//		writer->Write();
//	}
	
	return 0;
}
