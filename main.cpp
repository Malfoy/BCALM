#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include "lm.h"
#include "ograph.h"
#include "debug.h"

using namespace std;


int main(int argc, char ** argv)
{
	int sys(0);
	if(argc==1)
	{
		srand(time(NULL));
		bool b(true);
		int n(1000);
		while(b)
		{

			int sys(0);
			int k(20);
			int m(2);

			cout<<"Test "<<"n:"<<n<<" k:"<<k<<" m:"<<2*m<<endl;
			createinputlm(n,k,"randomgenome");
			sys+=system("cat randomgenome | sort | uniq > input.dot");
			cout<<"GO!"<<endl;
			graph G(k),G2(k);
			G.importg("input.dot");
			G.debruijn();
			G.compress();
			G.print("output.dot");
			cout<<"GO!"<<endl;
			createoutfile("input.dot","outputlowmemory.dot",k,m);
			if(!checkfile("output.dot","outputlowmemory.dot",k)){
				cout<<"Errors occurred !"<<endl;
				b=false;
			}
			else{
				cout<<"Success !"<<endl;
				n*=10;
				cout<<n<<endl;
			}
			b=false;
		}

	}
	if(argc==2)
	{
		string input(argv[1]);
		string output("compacted.dot");
		int m(5);
		//Should put 4 in case of not enought ulimit -n
		int k(detectk(input));
		createoutfile(input.c_str(),output.c_str(),k,m);
	}
	if(argc==3)
	{
		string input(argv[1]);
		string output(argv[2]);
		int m(5);
		//Should put 4 in case of not enought ulimit -n
		int k(detectk(input));
		createoutfile(input.c_str(),output.c_str(),k,m);
	}
	if(argc==4)
	{
		string input(argv[1]);
		string output(argv[2]);
		int m(atoi(argv[3])/2);
		int k(detectk(input));
		createoutfile(input.c_str(),output.c_str(),k,m);
	}
	if(argc==5)
	{
		string input(argv[1]);
		string output(argv[2]);
		int m(atoi(argv[3])/2);
		int k(atoi(argv[4]));

		createoutfile(input.c_str(),output.c_str(),k,m);
	}

	return sys;
}
