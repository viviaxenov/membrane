#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

#include"../header/membrane.h"


int main()
{
	Task test(0.00001, 0.01, 0.005, 60);
	test.SetFextPt(30, 30, 0, 0, 1);
	vtkSmartPointer<vtkStructuredGrid> sg = test.vtkSGrid();		
	vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

	writer->SetFileName("output/output.vts");
	writer->SetInputData(sg);
	writer->Write();

	char path[50];
	for(int i = 0; i < 400; i++)
	{
		test.Iteration();
		sg = test.vtkSGrid();		

		sprintf(path, "output/p%d.vts", i);
		writer->SetFileName(path);
		writer->SetInputData(sg);
		writer->Write();
	}
	test.SetFextPt(30, 30, 0, 0, 0.0);
	for(int i = 400; i < 800; i++)
	{
		test.Iteration();
		sg = test.vtkSGrid();		

		sprintf(path, "output/p%d.vts", i);
		writer->SetFileName(path);
		writer->SetInputData(sg);
		writer->Write();
	}

	return 0;
}
