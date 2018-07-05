#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

#include"./header/membrane.h"
#include"./header/uniform.h"


int main()
{
	UTask test(0.1, 0.05, 50, 1.0);
	puts("");

	vtkSmartPointer<vtkStructuredGrid> sg = test.vtkSGrid();		
	vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

	writer->SetFileName("output.vts");
	writer->SetInputData(sg);
	writer->Write();
	
	char path[20];	
	test.StartCondPoint(25, 25, 0.3, 0);

	for(int i = 0; i < 200; i++)
	{
		sg = test.vtkSGrid();		

		sprintf(path, "output/p%d.vts", i);
		writer->SetFileName(path);
		writer->SetInputData(sg);
		writer->Write();
	
		test.Iteration();
	}	

	return 0;
}
