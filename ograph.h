#ifndef OGRAPH
#define OGRAPH

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cstdlib>
#include <unordered_map>

using namespace std;

bool adjacent (const string& node1,const  string& node2,int k);

string readn(ifstream *file,uint64_t n);

int chartoint(char c);

string minimalsub(const string &w, const int &p,const int &k);

string minimalsub2(const string &w, const int &p,const int &k);

class neighbour
{
	public:
		array<uint64_t,4> son;
		array<uint64_t,4> father;
		
		uint64_t nbson();
		uint64_t nbfather();
		uint64_t gson();
		uint64_t gfather();
		
		void add(uint64_t p, array<uint64_t,4> *array);
		void addfather(uint64_t p);
		void addson(uint64_t p);
		void remove(uint64_t v);
		
		
};

class graph
{
	public:
		uint64_t n;
		int k;
		vector<string> nodes;
		unordered_multimap<uint64_t,uint64_t> map;
		vector<neighbour> neighbor;
		
		graph(const int ni)
		{
			k=ni;
			n=1;
			nodes.push_back("");
		}
		
		uint64_t getkey(string str);
		int weight();
		void addvertex(const string str);
		void debruijn();
		void compress(const string& str="");
		void importg(const char *name);
		void print(const char *name);
		void printedges(const char *name);
		
};

#endif
