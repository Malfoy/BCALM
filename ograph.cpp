#include "ograph.h"

/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */

using namespace std;



//read n character
string readn(ifstream *file,uint64_t n)
{
	string contents(n,'0');
	file->read(&contents[0],n);
	return(contents);
}



bool adjacent (const string& node1,const  string& node2,int k)
{
	return(node1.substr(node1.size()-k+1,k-1)==node2.substr(0,k-1));
}



int chartoint(char c)
{
	switch(c)
	{
		case 'a':
		return 0;
		case 'c':
		return 1;
		case 'g':
		return 2;
		case 't':
		return 3;
		default:
		cout<<"pb:"<<c<<endl;
		return 0;
		
	}
}

uint64_t stringtotint(const string& str)
{
	uint64_t res(0);
	for(uint64_t i(0);i<str.size();i++)
	{
		res<<=2;
		res+=chartoint(str[i]);
	}
	return res;
} 



//return left minimiser
string minimalsub2(const string &w, const int &p,const int &k)
{
	int s(w.size());
	string min2(w.substr(s-k+1,p));
	for(int i(s-k+2);i<s-p+1;i++)
		if(min2.compare(0,p,w,i,p)>0)
			min2=w.substr(i,p);
	return(min2);
}



//return right minimiser
string minimalsub(const string &w, const int &p,const int &k)
{
	string min(w.substr(0,p));
	for(int i(1);i<k-p;i++)
		if(min.compare(0,p,w,i,p)>0)
			min=w.substr(i,p);
	return(min);
}



void graph::addvertex(string str)
{
	n++;
	nodes.push_back(str);
	uint i(nodes.size()-1);
	uint64_t key(getkey(str));
	map.insert({key,i});
}



void graph::debruijn()
{
	vector<uint> empty;
	neighbor=vector<neighbour> (n);
	string node,kmmer;
	uint64_t key;
	for(uint64_t i(1);i<n;i++)
	{
		node=nodes[i];
		kmmer=node.substr(node.size()-k+1,k-1);
		key=getkey(kmmer);
		auto it(map.equal_range(key));
		for(auto j(it.first);j!=it.second;j++)
			//if k>32 collision can occur
			if(adjacent(node,nodes[j->second],k))
			{
				neighbor[i].addson(j->second);
				neighbor[j->second].addfather(i);
			}
	}
	map.clear();
}



//Compact the graph but not the nodes that should be compacted in an other bucket 
void graph::compress(const string& min)
{
	uint64_t p(min.size()),sonindice;
	string newnode,son,node;
	for(uint64_t nodeindice(1);nodeindice<n;nodeindice++)
		if(nodes[nodeindice]!="")
			if(neighbor[nodeindice].nbson()==1)
			{
				sonindice=(neighbor[nodeindice].gson());
				if(neighbor[sonindice].nbfather()==1)
				{
					son=nodes[sonindice];
					if(son!="")
					{
						node=(nodes[nodeindice]);
						if(minimalsub2(node,p,k)==min || p==0)
						{
							newnode=node+son.substr(k-1);
							nodes[sonindice]=newnode;
							nodes[nodeindice]="";
							neighbor[sonindice].father=neighbor[nodeindice].father;
							for(auto k(neighbor[nodeindice].father.begin());k!=neighbor[nodeindice].father.end();k++)
							{
								neighbor[*k].remove(nodeindice);
								neighbor[*k].addson(sonindice);
								
							}
						}
				}
				}
			}
}



//import graph from file
void graph::importg(const char *name)
{
	ifstream fichier(name,ios::in);
	string line,vertex;
	if(fichier)
		while(!fichier.eof())
		{
			getline(fichier,line); 
			vertex=line.substr(0,line.size()-1);
			if(vertex!="")
				addvertex(vertex);
		}
	else
		cerr<<"no file"<<endl;
}



int graph::weight()
{
	int res(0);
	for(auto k(nodes.begin());k!=nodes.end();k++)
		res+=k->size();
	return res;
}



void graph::print(const char *name)
{
	ofstream fichier(name, ios::out | ios::trunc);
	if(fichier)
		for(uint64_t i(1);i<n;i++)
		{
			string s(nodes[i]);
			if(s!="")
				fichier<<s<<";"<<endl;
		}
}

void graph::printedges(const char *name)
{
	ofstream fichier(name, ios::out | ios::trunc);
	if(fichier)
	{      
		fichier << "digraph test {" <<endl;
		for(uint64_t i(1);i<n;i++)
		{
			string s(nodes[i]);
			if(s!="")
			{
				fichier<<s<<";"<<endl;
				auto v=neighbor[i].son;
				for(auto j=v.begin();j!=v.end();j++)
					if(*j!=0)
						if(nodes[*j]!="")
							fichier<<s<<"->"<<nodes[*j]<<";"<<endl;
			}
		}
		fichier << "}"<<endl;
		fichier.close();
	}
}


uint64_t graph::getkey(string str)
{
	return stringtotint(str.substr(0,k-1));
}


void neighbour::add(uint64_t p, array<uint64_t,4> *array)
{
	if((*array)[0]==0)
		(*array)[0]=p;
	else
		if((*array)[1]==0)
			(*array)[1]=p;
		else
			if((*array)[2]==0)
				(*array)[2]=p;
			else
				(*array)[3]=p;
}

void neighbour::addfather(uint64_t p)
{
	add(p,&father);
}


void neighbour::addson(uint64_t p)
{
	add(p,&son);
}



uint64_t neighbour::nbson()
{
	uint64_t ret(0);
	for(int i(0);i<4;i++)
		if(son[i]!=0)
			ret++;
	return ret;
}



uint64_t neighbour::nbfather()
{
	uint64_t ret(0);
	for(int i(0);i<4;i++)
		if(father[i]!=0)
			ret++;
	return ret;
}



uint64_t neighbour::gson()
{
	for(int i(0);i<4;i++)
		if(son[i]!=0)
			return son[i];
	cout<<"bug neighbour"<<endl;		
	return 0;
}


//orphan
uint64_t neighbour::gfather()
{
	for(int i(0);i<4;i++)
		if(father[i]!=0)
			return father[i];
	return 0;
}



void neighbour::remove(uint64_t v)
{
	for(int i(0);i<4;i++)
		if(son[i]==v)
			son[i]=0;
}
