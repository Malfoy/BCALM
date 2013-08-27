#include "ograph.h"
#include <algorithm>

/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */

using namespace std;


string reversecompletment(const string& str ){
	string res;
	for(int i(str.size()-1); i > -1; i--){
		switch(str[i]){
			case 'g':
				res.push_back('c');
				break;
			case 'a':
				res.push_back('t');
				break;
			case 't':
				res.push_back('a');
				break;
			case 'c':
				res.push_back('g');
				break;
			case '+':
				res.push_back('-');
				break;
			case '-':
				res.push_back('+');
				break;
			default:
				res.push_back(str[i]);
				break;
		}
	}
	return res;
}

//read n character
string readn(ifstream *file,uint64_t n){
	string contents(n,'0');
	file->read(&contents[0],n);
	return(contents);
}



bool adjacent (const string& node1,const  string& node2,int k){
	return(node1.substr(node1.size()-k+1,k-1)==node2.substr(0,k-1));
}



int chartoint(char c){
	switch(c){
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

uint64_t stringtotint(const string& str){
	uint64_t res(0);
	for(uint64_t i(0);i<str.size();i++){
		res<<=2;
		res+=chartoint(str[i]);
	}
	return res;
} 

uint64_t stringtotintc(const string& str){
	uint64_t res(0);
	for(int64_t i(str.size()-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[i]);
	}
	return res;
} 



//return left minimiser
string minimalsub2(const string &w, const int &p,const int &k){
	int s(w.size());
	string min2(w.substr(s-k+1,p));
	for(int i(s-k+2);i<s-p+1;i++)
		if(min2.compare(0,p,w,i,p)>0)
			min2=w.substr(i,p);
	return(min2);
}



//return right minimiser
string minimalsub(const string &w, const int &p,const int &k){
	string min(w.substr(0,p));
	for(int i(1);i<k-p;i++)
		if(min.compare(0,p,w,i,p)>0)
			min=w.substr(i,p);
	return(min);
}



void graph::addvertex(string str){
	n++;
	nodes.push_back(str);
	uint64_t i(nodes.size()-1);
	uint64_t key(getkey(str));
	uint64_t keyrc(getkeyrevc(str));
	map.insert({key,i});
	maprev.insert({keyrc,i});
}



void graph::debruijn(){
	neighbor=vector<neighbour> (n);
	string node,kmmer,kmmerr;
	uint64_t key,keyrc;
	for(uint64_t i(1);i<n;i++){
		node=nodes[i];
		kmmer=node.substr(node.size()-k+1,k-1);
		kmmerr=node.substr(0,k-1);
		key=getkey(kmmer);
		keyrc=getkeyrevc(kmmerr);
		auto it(map.equal_range(key));
		for(auto j(it.first);j!=it.second;j++)
			//if k>32 collision can occur
			if(adjacent(node,nodes[j->second],k)){
				neighbor[i].add(j->second,1);
				neighbor[j->second].add(i,4);
			}
		it=(maprev.equal_range(key));
		for(auto j(it.first);j!=it.second;j++)
			if(adjacent(node,reversecompletment(nodes[j->second]),k)){
				neighbor[i].add(j->second,2);
				neighbor[j->second].add(i,2);
			}
		it=(map.equal_range(keyrc));
		for(auto j(it.first);j!=it.second;j++)
			if(adjacent(reversecompletment(node),nodes[j->second],k)){
				neighbor[i].add(j->second,3);
				neighbor[j->second].add(i,3);
			}
	}
	map.clear();
}


uint64_t graph::becompacted(uint64_t nodeindice,const string& min, unsigned char *type)
{
	*type=0;
	int m=min.size();
	string node=nodes[nodeindice];
	
	if(node.empty())
		return 0;
	
	auto neigh(neighbor[nodeindice]); 
	int one(neigh.nbtype(1)),two(neigh.nbtype(2)),three(neigh.nbtype(3)),four(neigh.nbtype(4));
	int in(three+four);
	int out(one+two);
	if(out==1 && (minimalsub2(node,m,k)==min || min.empty())){
		if(one==1){
			uint64_t sonindice(neigh.gtype(1));
			*type=1;
			if(neighbor[sonindice].nbtype(3)+neighbor[sonindice].nbtype(4)==1 && sonindice!=nodeindice)
				return sonindice;
		}
		else{
			uint64_t sonindice(neigh.gtype(2));
			*type=2;
			if(neighbor[sonindice].nbtype(2)+neighbor[sonindice].nbtype(1)==1 && sonindice!=nodeindice)
				return sonindice;
		}
	}
	if(in==1){
		if(three==1 && (minimalsub2(reversecompletment(node),m,k)==min || min.empty())){
			uint64_t sonindice(neigh.gtype(3));
			*type=3;
			if(neighbor[sonindice].nbtype(3)+neighbor[sonindice].nbtype(4)==1 && sonindice!=nodeindice)
				return sonindice;
		}
		if(four==1 && (minimalsub(node,m,k)==min || min.empty())){
			uint64_t sonindice(neigh.gtype(4));
			*type=4;
			if(neighbor[sonindice].nbtype(1)+neighbor[sonindice].nbtype(2)==1 && sonindice!=nodeindice)
				return sonindice;
		}
	}
	
	return 0;
}

void graph::reverse(int64_t with){
	string newnode(nodes[with]);
	int64_t indice,type;
	if(newnode>(reversecompletment(newnode)))	{
		nodes[with]=reversecompletment(newnode);
		for(auto i(0);i<8;i++)		{
			indice=neighbor[with].list[i].first;
			type=neighbor[with].list[i].second;
			
			if(type==1)			{
				neighbor[indice].removep(with,4);
				neighbor[with].removep(indice,1);
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
				continue;
			}
			if(type==2 )			{
				neighbor[indice].removep(with,2);
				neighbor[with].removep(indice,2);
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
				continue;
			}
			if(type==3 )			{
				neighbor[indice].removep(with,3);
				neighbor[with].removep(indice,3);
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
				continue;
			}
			if(type==4 )			{
				neighbor[indice].removep(with,1);
				neighbor[with].removep(indice,4);
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
				continue;
			}
		}
	}
}


void graph::compact(uint64_t nodeindice,uint64_t with, unsigned char c){
	
	string newnode,node(nodes[nodeindice]),son(nodes[with]);
	uint64_t indice;
	if(nodeindice==with || node.empty() || son.empty())
		return;
	unsigned char type;
	switch(c){
		case 1:
		newnode=node+son.substr(k-1);
		nodes[nodeindice]="";
		nodes[with]=newnode;
		for(int  i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].remove(nodeindice);
			if(indice==nodeindice){
				continue;
			}
			if(type==3 ){
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==4 ){
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
			}
		}
		if(newnode>reversecompletment(newnode)){
			reverse(with);
		}
		break;
		
		case 2:
		newnode=node+reversecompletment(son).substr(k-1);
		nodes[nodeindice]="";
		nodes[with]=newnode;
		for(auto i(0);i<8;i++){
			indice=neighbor[with].list[i].first;
			type=neighbor[with].list[i].second;
			neighbor[indice].remove(with);
			neighbor[with].remove(indice);
			if(type==3){
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
			}
			if(type==4){
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
			}
		}
		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].remove(nodeindice);
			if(type==3){
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==4 ){
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
			}
		}
		if(newnode>reversecompletment(newnode))
			reverse(with);
		break;
		
		case 3:
		newnode=reversecompletment(node)+son.substr(k-1);
		nodes[nodeindice]="";
		nodes[with]=newnode;
		neighbor[with].removep(nodeindice,3);
		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			if(type==1){
				neighbor[indice].removep(nodeindice,4);
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==2 ){
				neighbor[indice].removep(nodeindice,2);
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
				//warning echange =bug
			}
		}
		if(newnode>reversecompletment(newnode))
			reverse(with);
		break;
		
		case 4:
		newnode=son+node.substr(k-1);
		nodes[nodeindice]="";
		nodes[with]=newnode;
		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].removep(nodeindice,1);
			if(type==1){
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
			}
			if(type==2 ){
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
			}
		}
	}

}

//Compact the graph but not the nodes that should be compacted in an other bucket 
void graph::compress(const string& min)
{
	for(uint64_t nodeindice(1);nodeindice<n;nodeindice++){
		unsigned char type(0);
		uint64_t with=becompacted(nodeindice,min,&type);
		if(with!=0){
			compact(nodeindice,with,type);
		}
	}
}



//import graph from file
void graph::importg(const char *name){
	ifstream fichier(name,ios::in);
	string line,vertex;
	if(fichier)
		while(!fichier.eof()){
			getline(fichier,line); 
			vertex=line.substr(0,line.size()-1);
			if(vertex!="")
				addvertex(vertex);
		}
	else
		cerr<<"no file"<<endl;
}



int graph::weight(){
	int res(0);
	for(auto k(nodes.begin());k!=nodes.end();k++)
		res+=k->size();
	return res;
}



void graph::print(const char *name){
	ofstream fichier(name, ios::out | ios::trunc);
	if(fichier)
		for(uint64_t i(1);i<n;i++){
			string s(nodes[i]);
			if(s!="")
				fichier<<s<<";"<<endl;
		}
}

void graph::printedges(const char *name){
	ofstream fichier(name, ios::out | ios::trunc);
	if(fichier){      
		fichier << "digraph test {" <<endl;
		for(uint64_t i(1);i<n;i++){
			string s(nodes[i]);
			if(s!=""){
				fichier<<s<<";"<<endl;
				//~ auto v=neighbor[i].son;
				//~ for(auto j=v.begin();j!=v.end();j++)
					//~ if(j->first!=0)
						//~ if(nodes[j->first]!="")
							//~ fichier<<s<<"->"<<nodes[j->first]<<";"<<endl;
			}
		}
		fichier << "}"<<endl;
		fichier.close();
	}
}


uint64_t graph::getkey(string str){
	return stringtotint(str.substr(0,k-1));
}


uint64_t graph::getkeyrevc(string str){
	return stringtotintc(str.substr(str.size()-k+1,k-1));
}


void neighbour::add(uint64_t p,unsigned char c)
{
	for(int i(0);i<8;i++){
		if(list[i].first==0  ){
			list[i]=make_pair(p,c);
			return;
		}
		if(list[i].first==p && list[i].second==c ){
			return;
		}
	}
}



uint64_t neighbour::nbtype(unsigned char c)
{
	uint64_t ret(0);
	for(int i(0);i<8;i++)
		if(list[i].second==c){
			ret++;
		}
	return ret;
}




uint64_t neighbour::gtype(unsigned char c)
{
	for(int i(0);i<8;i++)
		if(list[i].second==c)
			return list[i].first;
	cout<<"bug neighbour"<<endl;		
	return 0;
}




unsigned char neighbour::remove(uint64_t v)
{
	for(int i(0);i<8;i++){
		if(list[i].first==v){
			list[i].first=0;
			list[i].second=0;
			return 0;
		}
	}
	return 0;
}

unsigned char neighbour::removep(uint64_t v,unsigned char c)
{
	for(int i(0);i<8;i++){
		if(list[i].first==v && list[i].second==c){
			list[i].first=0;
			list[i].second=0;
			return 0;
		}
	}
	return 0;
}


unsigned char neighbour::removetype(unsigned char c)
{
	for(int i(0);i<8;i++)
		if(list[i].second==c){
			list[i].first=0;
		}
	return 0;
}

