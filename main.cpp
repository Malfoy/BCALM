#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <chrono>
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
		uint64_t n(100000);
		int k(50);
		int m(4);
		
		cout<<"Test "<<"n:"<<n<<" k:"<<k<<" m:"<<2*m<<endl;
		createinputlm(n,k,"randomgenome");
		sys+=system("cat randomgenome | sort | uniq > input.dot");
		
		graph G(k);
		G.importg("input.dot");
		G.debruijn();
		G.compress();
		G.print("output.dot");
		createoutfile("input.dot","outputlowmemory.dot",k,m);
		if(!checkfile("output.dot","outputlowmemory.dot"))
		{
			cout<<"Errors occurred !"<<endl;
		}
		else
		{
			cout<<"Success !"<<endl;
		}
	}
	if(argc==2)
	{
		srand(time(NULL)); 
		uint64_t n(atoi(argv[1]));
		int k(50);
		int m(4);
		
		cout<<"Test "<<"n:"<<n<<" k:"<<k<<" m:"<<2*m<<endl;
		createinputlm(n,k,"randomgenome");
		sys+=system("cat randomgenome | sort | uniq > input.dot");
		
		graph G(k);
		G.importg("input.dot");
		G.debruijn();
		G.compress();
		G.print("output.dot");
		createoutfile("input.dot","outputlowmemory.dot",k,m);
		if(!checkfile("output.dot","outputlowmemory.dot"))
		{
			cout<<"Errors occurred !"<<endl;
		}
		else
		{
			cout<<"Success !"<<endl;
		}
	}
	if(argc==3)
	{
		srand(time(NULL)); 
		uint64_t n(atoi(argv[1]));
		int k(atoi(argv[2]));
		int m(4);
		
		cout<<"Test "<<"n:"<<n<<" k:"<<k<<" m:"<<2*m<<endl;
		createinputlm(n,k,"randomgenome");
		sys+=system("cat randomgenome | sort | uniq > input.dot");
		
		graph G(k);
		G.importg("input.dot");
		G.debruijn();
		G.compress();
		G.print("output.dot");
		createoutfile("input.dot","outputlowmemory.dot",k,m);
		if(!checkfile("output.dot","outputlowmemory.dot"))
		{
			cout<<"Errors occurred !"<<endl;
		}
		else
		{
			cout<<"Success !"<<endl;
		}
	}
	if(argc==4)
	{
		srand(time(NULL)); 
		uint64_t n(atoi(argv[1]));
		int k(atoi(argv[2]));
		int m(atoi(argv[3])/2);
		
		cout<<"Test "<<"n:"<<n<<" k:"<<k<<" m:"<<2*m<<endl;
		createinputlm(n,k,"randomgenome");
		sys+=system("cat randomgenome | sort | uniq > input.dot");
		
		graph G(k);
		G.importg("input.dot");
		G.debruijn();
		G.compress();
		G.print("output.dot");
		createoutfile("input.dot","outputlowmemory.dot",k,m);
		if(!checkfile("output.dot","outputlowmemory.dot"))
		{
			cout<<"Errors occurred !"<<endl;
		}
		else
		{
			cout<<"Success !"<<endl;
		}
	}
	if(argc==5)
	{
		string input(argv[1]);
		string output(argv[2]);
		int k(atoi(argv[3]));
		int m(atoi(argv[4])/2);
		createoutfile(input.c_str(),output.c_str(),k,m);
	}

	return sys;
}
