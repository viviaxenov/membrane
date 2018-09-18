#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

#include"../header/membrane.h"




int main(int argc, char *argv[])
{

	if(argc != 4)
	{
		printf("USAGE: TestVTK frames iter_per_frame output_dir\n");
		return -1;
	}

	Task test(TAU, H, N_CELLS);

	test.frames = atoi(argv[1]);
	test.iter_per_frame = atoi(argv[2]);
	test.SetOutputDir(string(argv[3]));


	double amp = 0.1;

	test.SetOffPt(N_CELLS/2, N_CELLS/2, 0, 0, amp*H);
	test.SetOffPt(N_CELLS/2 - 1, N_CELLS/2, 0, 0, amp*H);
	test.SetOffPt(N_CELLS/2, N_CELLS/2 - 1, 0, 0, amp*H);
	test.SetOffPt(N_CELLS/2 - 1, N_CELLS/2 - 1, 0, 0, amp*H);

	test.Execute();

	return 0;
}
