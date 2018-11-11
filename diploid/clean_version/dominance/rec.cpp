#include "fisher.h"
#include <vector>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <bitset>
#include <string>
#include <random>
using namespace std;

extern MTRand rnd;

// rec function: generates recombinant chromosome "res"
// from parental chromosomes "c1" and "c2"
// nbCo is the number of cross-overs,
// nS the number of selected loci

void rec(chr &res, chr &c1, chr &c2, double R, int nS)
{
	vector<int> pos;
	int j;
	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;

	res.sel.clear();

	// vector "pos" holds the positions of cross-overs:

    int nbCo = int(poisdev(R));
    for (j = 0; j < nbCo; j++)
		pos.push_back(rnd.randInt(nS));
	sort(pos.begin(), pos.end());

	// creates recombination mask:

	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 1 : 0);

	rec.resize(nS, (nbCo % 2) == 0 ? 1 : 0);

	// creates recombinant chromosome:

	off1 = (c1.sel & rec);
	rec.flip();
	off2 = (c2.sel & rec);
	res.sel = (off1 | off2);
}


// generate a random mask for recombination
// used only for free recombining option with option L = -1
boost::dynamic_bitset<> RandomMask(int N)
{
	int num = floor(N/64)+1;
	string s_mask = "";
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<unsigned long int> dis;
	for (int i = 0 ; i < num ; i++)
	{
		bitset<64> tt(dis(gen));
		s_mask += tt.to_string();
	}
	boost::dynamic_bitset<> maskfromstring(s_mask);
	maskfromstring.resize(N);
	return maskfromstring;
}


void freerec(chr &res, chr &c1, chr &c2, int nS)
{
	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;

	res.sel.clear();

	// recombination mask:
	rec = RandomMask(nS);

	off1 = (c1.sel & rec);
	rec.flip();
	off2 = (c2.sel & rec);

	res.sel = (off1 | off2);
	//cout << "freerec at work!";
}
