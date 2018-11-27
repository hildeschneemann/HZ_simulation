#include <iostream>
#include <cmath>
using namespace std;


int main()
{
	int geno;
	int lv = 4;

	for (int loc = 0; loc < lv; loc++)
	{
		geno = 0;
		for (int t = 0; t < pow(3,loc); t++)
		{
			for (int i = 0; i < 3; i++)//3 states: 0 homo P1, 1 hetero, 2 homo P2
			{
				for (int j = 0; j < pow(3, (lv - loc - 1)); j++)
				{
					cout << i;
					geno += 1;
				}
			}
		}
		cout << "\n";
	}

	return(0);

}
