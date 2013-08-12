#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cstdlib>
#include <iterator>
#include <ctime>
#include <unordered_set>
#include <sys/time.h>
#include <chrono>
#include <algorithm>
#include <sstream>
#include "ograph.h"


using namespace std;

int main()
{
	ifstream in("fa");
	ofstream out("lm");
	string data;
	string * genome=new string();
	int k=20;
	while(!in.eof())
	{
		getline(in,data);
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		*genome+=data;
	}
	for(int64_t i(0);i+k<genome->size();i++)
	{
		out<<genome->substr(i,k)<<";"<<endl;
	}
	return 0;
}

