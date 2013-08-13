#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <map>
#include <string>
#include <cstdlib>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <sys/time.h>
#include <chrono>
#include <algorithm>
#include <memory>
#include "ograph.h"

using namespace std;


void fastatodot(const char * name,const char * namout)
{
	ifstream in(name);
	ofstream out(namout);
	string data;
	while(!in.eof())
	{
		getline(in,data);
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		if(data!="")
		out<<data<<";"<<endl;
	}
}

string trunc(const string str, const int t)
{
	return(str.substr(0,t));
}

string minbutbiggerthan(const string& minimiser1,const string& minimiser2,string namebucket)
{
	if(minimiser1<minimiser2)
	{
		if(minimiser1<namebucket)
		{
			return(minimiser2);
		}
		else
		{
			return minimiser1;
		}
	}
	else
	{
		if(minimiser2<namebucket)
		{
			return(minimiser1);
		}
		else
		{
			return minimiser2;
		}
	}
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
	bool again=true;
	while(again && i!=size)
	{
		if(nextchar(&s[size-i-1]))
		{
			again=false;
		}
		i++;
	}
	return s;
}

//return minimiser
string minimiser(const string &s,const int &minimisersize)
{
	int size=s.size();
	string minimiser=s.substr(0,minimisersize);
	for(int i=1;i<size-minimisersize+1;i++)
	{
		if(minimiser.compare(0,minimisersize,s,i,minimisersize)>0)
		{
			minimiser=s.substr(i,minimisersize);
		}
	}
	return(minimiser);
}

//testing function
void createinputlm(int64_t lr,int k,const char *name)
{
	ofstream out(name,ios::trunc);
	int r;
	string c;
	string kmer(k,'a');
	for(int b=0;b<k;b++)
	{
		r=rand()%4;
		switch(r)
		{
			case 1:
			{
				kmer[b]='a';
				break;
			}
			case 2:
			{
				kmer[b]='c';
				break;
			}
			case 3:
			{
				kmer[b]='g';
				break;
			}
			case 0:
			{
				kmer[b]='t';
				break;
			}
		}
	}
	for(int64_t b=0;b<lr;b++)
	{
		kmer=kmer.substr(1,k-1);
		r=rand()%4;
		switch(r)
		{
			case 1:
			{
				c='a';
				break;
			}
			case 2:
			{
				c='c';
				break;
			}
			case 3:
			{
				c='g';
				break;
			}
			case 0:
			{
				c='t';
				break;
			}
		}
		kmer+=c;
		out<<kmer<<";"<<endl;
	}
}

//Put kmers in superbuckets
void sortentry(string namefile,const int k,const int m)
{
	string superbucketname(m,'a');
	int numbersuperbucket=pow(4,m);
	ofstream out[256];
	//~ ofstream out[1024];
	for(int i(0);i<numbersuperbucket;i++)
	{
		out[i].open("z"+superbucketname,ofstream::app);
		superbucketname=nextstring(superbucketname);
	}
	
	ifstream in(namefile); 
	in.seekg(0,ios_base::end);
	int64_t size=in.tellg();
	in.seekg(0,ios::beg);
	int buffsize=10000*(k+2);
	int n=size/buffsize;
	string* buffer=new string();
	string kmer,min,prefix; 
	
	for(int j(0);j<n;j++)
	{
		*buffer=readn(&in,buffsize);
		for(uint64_t i=0;i<buffer->size();i+=k+2)
		{
			kmer=buffer->substr(i,k);
			min=minimiser(kmer,2*m);
			prefix=min.substr(0,m);
			out[stringtotint(prefix)]<<kmer<<";";
		}
	}
	*buffer=readn(&in,size-n*buffsize);
	
	for(uint64_t i=0;i<buffer->size();i+=k+2)
	{
		kmer=buffer->substr(i,k);
		min=minimiser(kmer,2*m);
		prefix=min.substr(0,m);
		out[stringtotint(prefix)]<<kmer<<";";
	}
	
	delete buffer;
}



//copy n character from in to out
void copylm(ifstream* in,uint64_t l,ofstream* out)
{
	int buffsize=1000000;
	int n=l/buffsize;
	string* str=new string();
	for(int j(0);j<n;j++)
	{
		*str=readn(in,buffsize);
		*out<<*str;
	}
	*str=readn(in,l-n*buffsize);
	*out<<*str;
	delete(str);
}

//Put nodes from superbuckets to buckets
void createbucket(const string bname,const int k,const int m)
{
	ifstream in("z"+bname); 
	if(!in.is_open())
	{
		cout<<"createbucket fail"<<endl;
		return;
	}
	in.seekg(0,ios_base::end );
	int64_t size=in.tellg();
	if(size==0)
	{
		remove(("z"+bname).c_str());
		return;
	}
	in.seekg(0,ios::beg);
	int buffsize=1000000;
	int numberbuffer=size/buffsize;
	
	string fname(m,'a');
	int nb=pow(4,m);
	ofstream out[256];
	//~ ofstream out[1024];
	for(int i(0);i<nb;i++)
	{
		out[i].open(bname+fname,ofstream::app);
		fname=nextstring(fname);
	}
	
	int64_t lastposition(-1),position(0);
	int64_t point(0);
	string kmerend;
	string * buffer= new string();
	string kmerbeg=readn(&in,k-1);
	in.seekg(0);
	string minleft,minright,min,prefix;
	for(int j(0);j<numberbuffer;j++)
	{
		*buffer=readn(&in,buffsize);
		point+=buffsize;
		for (uint64_t i(0); i<buffer->size(); i++,position++) 
		{
			if((*buffer)[i]==';')
			{
				if(i>=(uint64_t)k-1)
				{
					kmerend=buffer->substr(i-k+1,k-1);
				}
				else
				{
					in.seekg(position-k+1);
					kmerend=readn(&in,k-1);
				}
				minleft=minimalsub(kmerbeg,2*m,k);
				minright=minimalsub2(kmerend,2*m,k);
				min=minbutbiggerthan(minleft,minright,bname);
				prefix=min.substr(m,m);
				in.seekg(lastposition+1,ios_base::beg);
				copylm(&in,position-lastposition,&out[stringtotint(prefix)]);
				lastposition=position;
				if(buffer->size()>i+k)
				{
					kmerbeg=buffer->substr(i+1,k-1);
				}
				else
				{
					in.seekg(position+1);
					kmerbeg=readn(&in,k-1);
				}
				
			}
		}
		in.seekg(point);
	}
	if(size-numberbuffer*buffsize-1!=-1)
	{
		*buffer=readn(&in,size-numberbuffer*buffsize-1);
		*buffer+=";";
		for (uint64_t i(0); i<buffer->size(); i++,position++) 
		{
			if((*buffer)[i]==';')
			{
				if(i>=(uint64_t)k-1)
				{
					kmerend=buffer->substr(i-k+1,k-1);
				}
				else
				{
					in.seekg(position-k+1);
					kmerend=readn(&in,k-1);
				}
				minleft=minimalsub(kmerbeg,2*m,k);
				minright=minimalsub2(kmerend,2*m,k);
				min=minbutbiggerthan(minleft,minright,bname);
				prefix=min.substr(m,m);
				in.seekg(lastposition+1,ios_base::beg);
				copylm(&in,position-lastposition,&out[stringtotint(prefix)]);
				lastposition=position;
				if(buffer->size()>i+k)
				{
					kmerbeg=buffer->substr(i+1,k-1);
				}
				else
				{
					in.seekg(position+1);
					kmerbeg=readn(&in,k-1);
				}
			}
		}
	}
	delete(buffer);
	remove(("z"+bname).c_str());
}

//count the length of each node
vector<int64_t> countbucket(const string name)
{
	vector<int64_t> count;
	ifstream in(name); 
	if(in)
	{
		in.seekg( 0 , ios_base::end );
		int64_t size=in.tellg();
		if(size<2)
		{
			return count;
		}
		in.seekg(0,ios::beg);
		int buffsize=1000000;
		int numberbuffer=size/buffsize;
		
		int64_t lastposition(-1),position(0);
		string* buffer=new string ();
		for(int j(0);j<numberbuffer;j++)
		{
			*buffer=readn(&in,buffsize);
			for (uint64_t i(0); i<buffer->size(); i++,position++) 
			{
				if((*buffer)[i]==';')
				{
					count.push_back(position-lastposition-1);
					lastposition=position;
				}
			}
		}
		*buffer=readn(&in,size-numberbuffer*buffsize);
		for (uint64_t i(0); i<buffer->size(); i++,position++) 
		{
			if((*buffer)[i]==';')
			{
				count.push_back(position-lastposition-1);
				lastposition=position;
			}
		}
		delete buffer;
	}
	else
	{
	}
	return count;
}

//true iff node does not contain tag
bool notag(const string& node,const int64_t start,int64_t* n)
{
	for(uint64_t i(start);i<node.size();i++)
	{
		switch (node[i]) 
		{
			case '0' :
			*n=i;
			return false;
			case '1' :
			*n=i;
			return false;
			case '2' :
			*n=i;
			return false;
			case '3' :
			*n=i;
			return false;
			case '4' :
			*n=i;
			return false;
			case '5' :
			*n=i;
			return false;
			case '6' :
			*n=i;
			return false;
			case '7' :
			*n=i;
			return false;
			case '8' :
			*n=i;
			return false;
			case '9' :
			*n=i;
			return false;
			default :
			break;
			
		}
	}
	return true;
}

//return length of tag
int taglength(const string& node, int64_t j)
{
	int n=1;
	for(uint64_t i(j+1);i<node.size();i++)
	{
		switch (node[i]) 
		{
			case '0' :
			n++;
			break;
			case '1' :
			n++;
			break;
			case '2' :
			n++;
			break;
			case '3' :
			n++;
			break;
			case '4' :
			n++;
			break;
			case '5' :
			n++;
			break;
			case '6' :
			n++;
			break;
			case '7' :
			n++;
			break;
			case '8' :
			n++;
			break;
			case '9' :
			n++;
			break;
			default:
			return n;
		}
	}
	return n;
}

//Write a node remplacing tags by their sequences
void writeit(const string& outfile,const string& node,vector<pair<int64_t,int64_t>>* tagsposition,ifstream* tagfile,int64_t j)
{
	ofstream out(outfile,ios::app);
	if(out)
	{
		int64_t lastposition(0),tag,tagl,position,length;
		do
		{
			out<<node.substr(lastposition,j-lastposition);
			tagl=taglength(node,j);
			tag=stoi(node.substr(j,tagl));
			lastposition=j+tagl;
			auto pair=(*tagsposition)[tag];
			position=pair.first;
			length= pair.second;
			tagfile->seekg(position,ios_base::beg);
			copylm(tagfile,length,&out);
		}
		while(!notag(node,lastposition,&j));
		if(outfile!="inter.dot")
		{
			out<<node.substr(lastposition)<<";";
		}
		else
		{
			out<<node.substr(lastposition)<<";"<<endl;
		}
	}
	else
	{
		cerr<<"writeitbug"<<endl;
	}
}

//Decide where to put a node 
void goodplace(const string& node,int k, const string bucketname,ofstream* out,vector<pair<int64_t,int64_t>>* tagsposition,ifstream* tagfile,const int m)
{
	string suffix=bucketname.substr(0,m);
	string next=nextstring(suffix);
	string begin(m,'a');
	if(next==begin)
	{
		next=suffix+"ttttttt";
	}
	
	string leftminimiser=minimalsub(node,2*m,k);
	string rightminimiser=minimalsub2(node,2*m,k);
	string leftprefix='z'+trunc(leftminimiser,m);
	string rightprefix='z'+trunc(rightminimiser,m);
	int64_t i;
	
	if(notag(node,0,&i))
	{
		if(leftminimiser<rightminimiser)
		{
			if(leftminimiser>bucketname)
			{
				if(leftminimiser>=next)
				{
					ofstream fichierb(leftprefix,ios::app);
					fichierb<<node<<";";
					fichierb.close();
				}
				else
				{
					ofstream fichierb(leftminimiser,ios::app);
					fichierb<<node<<";";
					fichierb.close();
				}
			}
			else
			{
				if(rightminimiser>bucketname)
				{
					if(rightminimiser>=next)
					{
						ofstream fichierb(rightprefix,ios::app);
						fichierb<<node<<";";
						fichierb.close();
					}
					else
					{
						ofstream fichierb(rightminimiser,ios::app);
						fichierb<<node<<";";
						fichierb.close();
					}
				}
				else
				{
					*out<<node<<";"<<endl;
				}
			}
		}
		else
		{
			if(rightminimiser>bucketname)
			{
				if(rightminimiser>=next)
				{
					ofstream fichierb(rightprefix,ios::app);
					fichierb<<node<<";";
					fichierb.close();
				}
				else
				{
					ofstream fichierb(rightminimiser,ios::app);
					fichierb<<node<<";";
					fichierb.close();
				}
			}
			else
			{
				if(leftminimiser>bucketname)
				{
					if(leftminimiser>=next)
					{
						ofstream fichierb(leftprefix,ios::app);
						fichierb<<node<<";";
						fichierb.close();
					}
					else
					{
						ofstream fichierb(leftminimiser,ios::app);
						fichierb<<node<<";";
						fichierb.close();
					}
				}
				else
				{
					*out<<node<<";"<<endl;
				}
			}
		}
	}
	else
	{
		if(leftminimiser<rightminimiser)
		{
			if(leftminimiser>bucketname)
			{
				if(leftminimiser>=next)
				{
					writeit(leftprefix,node,tagsposition,tagfile,i);
				}
				else
				{
					writeit(leftminimiser,node,tagsposition,tagfile,i);
				}
			}
			else
			{
				if(rightminimiser>bucketname)
				{
					if(rightminimiser>=next)
					{
						writeit(rightprefix,node,tagsposition,tagfile,i);
					}
					else
					{
						writeit(rightminimiser,node,tagsposition,tagfile,i);
					}
				}
				else
				{
					writeit("inter.dot",node,tagsposition,tagfile,i);
				}
			}
		}
		else
		{
			if(rightminimiser>bucketname)
			{
				if(rightminimiser>=next)
				{
					writeit(rightprefix,node,tagsposition,tagfile,i);
				}
				else
				{
					writeit(rightminimiser,node,tagsposition,tagfile,i);
				}
			}
			else
			{
				if(leftminimiser>bucketname)
				{
					if(leftminimiser>=next)
					{
						writeit(leftprefix,node,tagsposition,tagfile,i);
					}
					else
					{
						writeit(leftminimiser,node,tagsposition,tagfile,i);
					}
				}
				else
				{
					writeit("inter.dot",node,tagsposition,tagfile,i);
				}
			}
		}
	}
}

//Compact a bucket and put the nodes on the right place
void compactbucket(const string prefix,const string suffix,const int k,const char *nameout,const int m)
{
	int buffsize=k;
	string fullname=prefix+suffix;
	auto count=countbucket(fullname);
	if(count.size()==0)
	{
		return make_pair(sizebucket,nbnode);
	}
	ifstream in(fullname);
	ofstream tagfile("tags");
	ofstream out(nameout,ios_base::app);
	graph* g=new graph(k);
	int tagnumber=0;
	vector<pair<int64_t,int64_t>> tagsposition;
	int64_t postags=0;
	string node,tag;
	
	if(in && tagfile && out)
	{
		for(auto it=count.begin();it!=count.end();it++)
		{
			int64_t length=*it;
			if(length<=2*buffsize)
			{
				node=readn(&in,length+1);
				g->addvertex(node.substr(0,length));
				sizebucket+=length;
				nbnode++;
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
				string t2=readn(&in,buffsize+1);
				node+=t2.substr(0,buffsize);
				g->addvertex(node);
				sizebucket+=node.size();
				nbnode++;
			}
		}
		tagfile.close();
	}
	else
	{
	}
	
	remove(fullname.c_str());
	g->debruijn();
	g->compress(fullname);
	ifstream fichiertagin("tags");
	for(auto it=g->nodes.begin();it!=g->nodes.end();it++)
	{
		if(it->size()!=0)
		{
			goodplace(*it,k,fullname,&out,&tagsposition,&fichiertagin,m);
		}
	}
	
	delete g;
	remove("tags");
	return make_pair(sizebucket,nbnode);
}

//debug function
bool checkfile(string name1, string name2,int k)
{
	int fail(0);
	ifstream t1(name1), t2(name2);
	string line;
	unordered_map<string,bool> s1,s2,e1,e2;
	while(!t1.eof())
	{
		getline(t1,line);
		if(line.size()>2)
		s1.insert(make_pair(line,false));
	}
	while(!t2.eof())
	{
		getline(t2,line);
		if(line.size()>2)
		s2.insert(make_pair(line,false));
	}
	for(auto it=s1.begin(); it!=s1.end(); it++)
	{
			string tf=it->first;
			auto smt=s2.find(tf);
			if(smt==s2.end())
			{
				string str=it->first;
				for(uint64_t i(0);i+k<str.size();i++)
				{
					e1.insert(make_pair(str.substr(i,k),false));
				}
			}
			else
			{
				smt->second=true;
			}
	}
	for(auto it=s2.begin(); it!=s2.end(); it++)
	{
		if(!it->second)
		{
			string str=it->first;
			for(uint64_t i(0);i+k<str.size();i++)
			{
				e2.insert(make_pair(str.substr(i,k),false));
			}
		}
	}
	for(auto it=e1.begin(); it!=e1.end(); it++)
	{
			string tf=it->first;
			auto smt=e2.find(tf);
			if(smt==s2.end())
			{
				fail++;
			}
			else
			{
				smt->second=true;
			}
	}
	for(auto it=e2.begin(); it!=e2.end(); it++)
	{
		if(!it->second)
		{
			fail++;
		}
	}
	cout<<"fail:"<<fail<<endl;
	return(fail==0);
}

//Create a file with the nodes of the compacted graph
void createoutfile(const char *namein,const char *nameout,const int k,const int m)
{
	int64_t max=0;
	int64_t nb=0;
	remove(nameout);
	auto temps1 = chrono::system_clock::now();
	sortentry(namein,k,m);
	auto temps2 = chrono::system_clock::now();
	chrono::nanoseconds timesort = (temps2 - temps1);
	chrono::nanoseconds timecreatebucket,timecompactbucket;
	auto nano=timesort.count();
	cout<<"sort "<<nano/1000000000<<" secondes."<<endl;
	nano=0;
	string seed(m,'a');
	string next,prefix,suffix;
	prefix=seed;
	next=nextstring(prefix);
	for(int i(0);i<pow(4,m);i++)
	{
		temps1 = chrono::system_clock::now();
		createbucket(prefix,k,m);
		temps2 = chrono::system_clock::now();
		auto nbrSecondes = (temps2 - temps1);
		timecreatebucket+=nbrSecondes;
		suffix=string(m,'a');
		temps1 = chrono::system_clock::now();
		for(int j(0);j<pow(4,m);j++)
		{
			auto p=compactbucket(prefix,suffix,k,nameout,m);
			if(max<p.first)
			{
				max=p.first;
			}
			if(nb<p.second)
			{
				nb=p.second;
			}
			suffix=nextstring(suffix);
		}
		temps2 = chrono::system_clock::now();
		nbrSecondes = (temps2 - temps1);
		timecompactbucket+=nbrSecondes;
		prefix=nextstring(prefix);
		if((i+1)%10==0)
			cout<<'-'<<flush;
	}	
	cout<<endl;
	nano=timecreatebucket.count();
	cout<<"createbucket "<<nano/1000000000<<"second"<<endl;
	nano=timecompactbucket.count();
	cout<<"compactbucket "<<nano/1000000000<<"second"<<endl;
	cout<<"maxbucket "<<max<<endl;
	cout<<"sizebucket "<<nb<<endl;
}

