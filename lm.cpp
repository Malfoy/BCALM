#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <chrono>
#include <algorithm>
#include "ograph.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>

/*
 * Constuct the compacted de bruijn graph from list of distinct kmers
 */

using namespace std;



int stringtointb(const string& str)
{
	int res(0);
	for(uint64_t i(0);i<str.size();i++)
	{
		res<<=2;
		res+=chartoint(str[i]);
	}
	return min(res,1019);
} 



string minbutbiggerthan(const string& minimiser1,const string& minimiser2,string namebucket)
{
	return (minimiser1<minimiser2) ? ((minimiser1<namebucket) ? minimiser2 : minimiser1) : ((minimiser2<namebucket) ? minimiser1 : minimiser2);
}



bool nextchar(char * c)
{
	switch(*c)
	{
		case 'a':
		*c='c';
		return true;
		case 'c':
		*c='g';
		return true;
		case 'g':
		*c='t';
		return true;
		case 't':
		*c='a';
		return false;
		default :
		cout<<"pb";
		return false;
	}
}



string nextstring(string s)
{
	int i(0),size(s.size());
	bool again(true);
	while(again && i!=size)
	{
		if(nextchar(&s[size-i-1]))
			again=false;
		i++;
	}
	return s;
}



//return minimiser
string minimiser(const string &node,const int &minimisersize)
{
	string minimiser(node.substr(0,minimisersize));
	for(uint64_t i(1);i<node.size()-minimisersize+1;i++)
		if(minimiser.compare(0,minimisersize,node,i,minimisersize)>0)
			minimiser=node.substr(i,minimisersize);
	return(minimiser);
}



//Put kmers in superbuckets
void sortentry(string namefile,const int k,const int m)
{
	string superbucketname(m,'a');
	int numbersuperbucket(min ((int)pow(4,m),1020));
	ofstream out[1020];
	for(int i(0);i<numbersuperbucket;i++)
	{
		out[i].open(".bcalmtmp/z"+superbucketname,ofstream::app);
		superbucketname=nextstring(superbucketname);
	}
	
	ifstream in(namefile); 
	in.seekg(0,ios_base::end);
	int64_t size(in.tellg()), buffsize(10000*(k+2)), n(size/buffsize);
	in.seekg(0,ios::beg);
	string* buffer=new string();
	string kmer,min,prefix; 
	
	for(int j(0);j<=n;j++)
	{
		if(j==n)
			*buffer=readn(&in,size-n*buffsize);
		else
			*buffer=readn(&in,buffsize);
		for(uint64_t i(0);i<buffer->size();i+=k+2)
		{
			kmer=buffer->substr(i,k);
			min=minimiser(kmer,2*m);
			prefix=min.substr(0,m);
			out[stringtointb(prefix)]<<kmer<<";";
		}
	}
	
	delete buffer;
}



//copy n characters from in to out
void copylm(ifstream* in,uint64_t n,ofstream* out)
{
	int buffsize(1000000),nbbuffer(n/buffsize);
	string* str=new string();
	for(int j(0);j<nbbuffer;j++)
	{
		*str=readn(in,buffsize);
		*out<<*str;
	}
	*str=readn(in,n-nbbuffer*buffsize);
	*out<<*str;
	delete(str);
}


//Put nodes from superbuckets to buckets
void createbucket(const string superbucketname,const int k,const int m)
{
	ifstream in(".bcalmtmp/z"+superbucketname); 
	if(!in.is_open())
	{
		cerr<<"createbucket fail"<<endl;
		return;
	}
	in.seekg(0, ios_base::end);
	int64_t size(in.tellg()), buffsize(1000000), numberbuffer(size/buffsize),nb(min((int)pow(4,m),1020));;
	if(size==0)
	{
		remove((".bcalmtmp/z"+superbucketname).c_str());
		return;
	}
	in.seekg(0,ios::beg);
	string suffix(m,'a');
	ofstream out[1020];
	for(int i(0);i<nb;i++)
	{
		out[i].open(".bcalmtmp/"+superbucketname+suffix,ofstream::app);
		suffix=nextstring(suffix);
	}
	
	int64_t lastposition(-1),position(0),point(0);
	string * buffer= new string();
	string kmerbeg(readn(&in,k-1)),kmerend,minleft,minright,min,prefix;
	in.seekg(0);
	
	for(int j(0);j<=numberbuffer;j++)
	{
		if(j==numberbuffer )
			if(size-numberbuffer*buffsize-1!=-1)
			{
				*buffer=readn(&in,size-numberbuffer*buffsize-1);
				*buffer+=";";
			}
			else
				*buffer="";
		else
		{
			*buffer=readn(&in,buffsize);
			point+=buffsize;
		}
		for (uint64_t i(0); i<buffer->size(); i++,position++) 
			if((*buffer)[i]==';')
			{
				if(i>=(uint64_t)k-1)
					kmerend=buffer->substr(i-k+1,k-1);
				else
				{
					in.seekg(position-k+1);
					kmerend=readn(&in,k-1);
				}
				minleft=minimalsub(kmerbeg,2*m,k);
				minright=minimalsub2(kmerend,2*m,k);
				min=minbutbiggerthan(minleft,minright,superbucketname);
				prefix=min.substr(m,m);
				in.seekg(lastposition+1,ios_base::beg);
				copylm(&in,position-lastposition,&out[stringtointb(prefix)]);
				lastposition=position;
				if(buffer->size()>i+k)
					kmerbeg=buffer->substr(i+1,k-1);
				else
				{
					in.seekg(position+1);
					kmerbeg=readn(&in,k-1);
				}
				
			}
		in.seekg(point);
	}
	
	delete(buffer);
	remove((".bcalmtmp/z"+superbucketname).c_str());
}



//count the length of each node
vector<int64_t> countbucket(const string& name)
{
	vector<int64_t> count;
	ifstream in(".bcalmtmp/"+name); 
	if(in)
	{
		in.seekg( 0 , ios_base::end );
		int64_t size(in.tellg()),buffsize(1000000),numberbuffer(size/buffsize),lastposition(-1),position(0);
		if(size<2)
			return count;
		in.seekg(0,ios::beg);
		string* buffer=new string ();
		
		for(int j(0);j<numberbuffer;j++)
		{
			*buffer=readn(&in,buffsize);
			for (uint64_t i(0); i<buffer->size(); i++,position++) 
				if((*buffer)[i]==';')
				{
					count.push_back(position-lastposition-1);
					lastposition=position;
				}
		}
		
		*buffer=readn(&in,size-numberbuffer*buffsize);
		for (uint64_t i(0); i<buffer->size(); i++,position++) 
			if((*buffer)[i]==';')
			{
				count.push_back(position-lastposition-1);
				lastposition=position;
			}
			
		delete buffer;
	}
	return count;
}



//true iff node does not contain tag
bool notag(const string& node,const int64_t start,int64_t* n)
{
	for(uint64_t i(start);i<node.size();i++)
	{		
		if (node[i] >= '0' && node[i] <= '9')
		{
			*n=i;
			return false;
		}
	}

	return true;
}



//return length of tag
int taglength(const string& node, int64_t j)
{
	int n=1;
	for(uint64_t i(j+1);i<node.size();i++)
		if (node[i] >= '0' && node[i] <= '9')
			n++;
		else
			return n;
	return n;
}



//Write a node remplacing tags by their sequences
void writeit(const string& outfile,const string& node,vector<pair<int64_t,int64_t>>* tagsposition,ifstream* tagfile,int64_t j,const string& fout)
{
	ofstream out(".bcalmtmp/"+outfile,ios::app);
	if(out)
	{
		int64_t lastposition(0),tag,tagl,position,length;
		pair<int64_t,int64_t> pair;
		do
		{
			out<<node.substr(lastposition,j-lastposition);
			tagl=taglength(node,j);
			tag=stoi(node.substr(j,tagl));
			lastposition=j+tagl;
			pair=(*tagsposition)[tag];
			position=pair.first;
			length= pair.second;
			tagfile->seekg(position,ios_base::beg);
			copylm(tagfile,length,&out);
		}
		while(!notag(node,lastposition,&j));
		if(outfile!=fout)
			out<<node.substr(lastposition)<<";";
		else
			out<<node.substr(lastposition)<<";"<<endl;
	}
	else
		cerr<<"writeitbug"<<endl;
}

void put(const string& outfile,const string& node,const string& fout)
{
	ofstream out(".bcalmtmp/"+outfile,ios::app);
	out<<node<<";";
	if(outfile==fout)
	{
		out<<endl;
	}
}

void putorwrite(const string& outfile, const string& node, vector<pair<int64_t,int64_t>>* tagsposition , ifstream* tagfile,const string& fout)
{
	int64_t i;
	if(notag(node,0,&i))
	{
		put(outfile,node,fout);
	}
	else
	{
		writeit(outfile,node,tagsposition,tagfile,i,fout);
	}
}

//Decide where to put a node 
void goodplace(const string& node,int k, const string& bucketname,vector<pair<int64_t,int64_t>>* tagsposition,ifstream* tagfile,const int m,const string& nameout)
{
	string suffix(bucketname.substr(0,m)),next(nextstring(suffix)),begin(m,'a');
	if(next==begin)
		next=suffix+"ttttttt";
	string leftminimiser(minimalsub(node,2*m,k)),rightminimiser(minimalsub2(node,2*m,k));
	string minimum(min(leftminimiser,rightminimiser)),maximum(max(leftminimiser,rightminimiser));
	string minprefix('z'+minimum.substr(0,m)),maxprefix('z'+maximum.substr(0,m));
	
	if(minimum>bucketname)
		putorwrite( ((minimum>=next) ? minprefix : minimum), node,tagsposition,tagfile,nameout);
	else
		if(maximum>bucketname)
			putorwrite( ((maximum>=next) ? maxprefix : maximum), node,tagsposition,tagfile,nameout);
		else
			putorwrite( nameout, node,tagsposition,tagfile,nameout);
}



//Compact a bucket and put the nodes on the right place
void compactbucket(const string prefix,const string suffix,const int k,const char *nameout,const int m)
{
	int64_t buffsize(2*k),postags(0),length;
	long long tagnumber(0);
	string fullname(prefix+suffix),node,tag,end;
	auto count(countbucket(fullname));
	if(count.size()==0)
	{
		remove((".bcalmtmp/"+fullname).c_str());
		return;
	}
	ifstream in(".bcalmtmp/"+fullname);
	ofstream tagfile(".bcalmtmp/tags"),out(".bcalmtmp/"+(string)nameout,ios_base::app);
	graph* g=new graph(k);
	vector<pair<int64_t,int64_t>> tagsposition;
	
	if(in && tagfile && out)
	{
		for(auto it=count.begin();it!=count.end();it++)
		{
			length=*it;
			if(length<=2*buffsize)
			{
				node=readn(&in,length+1);
				g->addvertex(node.substr(0,length));
			}
			else
			{
				node=readn(&in,buffsize);
				tag=to_string(tagnumber);
				node+=tag;
				tagnumber++;
				copylm(&in,length-2*buffsize,&tagfile);
				tagsposition.push_back(make_pair(postags,length-2*buffsize));
				postags+=length-2*buffsize;
				end=readn(&in,buffsize+1);
				node+=end.substr(0,buffsize);
				g->addvertex(node);
			}
		}
		tagfile.close();
	}
	
	remove((".bcalmtmp/"+fullname).c_str());
	g->debruijn();
	g->compress(fullname);
	ifstream fichiertagin(".bcalmtmp/tags");
	for(auto it(g->nodes.begin());it!=g->nodes.end();it++)
		if(it->size()!=0)
			goodplace(*it,k,fullname,&tagsposition,&fichiertagin,m,nameout);
	
	delete g;
	remove(".bcalmtmp/tags");
	return;
}

//Create a file with the nodes of the compacted graph
void createoutfile(const char *namein,const char *nameout,const int k,const int m)
{
	int sys(0);
	mkdir(".bcalmtmp",0777);
	sys+=system("rm -rf .bcalmtmp/*");
	int64_t nbsuperbucket(min((int)pow(4,m),1020));
	remove(nameout);
	sortentry(namein,k,m);
	string seed(m,'a');
	string next,prefix,suffix;
	prefix=seed;
	next=nextstring(prefix);
	for(int i(0);i<nbsuperbucket;i++)
	{
		createbucket(prefix,k,m);
		suffix=string(m,'a');
		for(int j(0);j<pow(4,m);j++)
		{
			compactbucket(prefix,suffix,k,nameout,m);
			suffix=nextstring(suffix);
		}
		prefix=nextstring(prefix);
		if((i+1)%10==0)
			cout<<'-'<<flush;
	}	
	cout<<endl;
	sys+=system(("mv .bcalmtmp/"+(string)nameout+" "+(string)nameout).c_str());
	if(sys!=0)
		cerr<<"system call failed"<<endl;
}

