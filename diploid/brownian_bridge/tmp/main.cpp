// main() function: reads parameter values in input file,
// and runs the simulation.

#include <iostream>
using namespace std;

// input and output files:

FILE * pFile;

// random number generator (Mersenne Twister):


int main(int argc, char * argv[])
{
	// definitions of variables
	int lv = 8;
	int nv = 2;
	int mut1 = 0;
	int mut2 = 0;
	double mut3 = 0.0;
	//double dom = 0.0;
	double * muteff = new double [lv];
	pFile = fopen(argv[1], "r");
	for (int i = 0; i < (lv * nv); i++)
	{
		fscanf(pFile, "%d", &mut1);
		fscanf(pFile, "%d", &mut2);
		fscanf(pFile, "%lf", &mut3);
		//fscanf(pFile, "%lf", &dom);
		cout << "mut1: " << mut1 << ", mut2: " << mut2 << ", mut3: " << mut3 << "\n";	
	}
	fclose(pFile);
	return(0);
}
