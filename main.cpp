#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include "lm.h"
#include "ograph.h"
#include "debug.h"

using namespace std;

int main(int argc, char ** argv)
{
	int sys(0);
	if(argc==1){
        printf("usage: <input> [output.dot] [minimizer length]\n");
        printf("Note: maximum minimizer length (minimizer length = 10) requires that you type 'ulimit -n 1100' in your shell prior to running bcalm, else the software will crash\n");
        exit(1);
	}
	if(argc==2){
		string input(argv[1]);
		string output("compacted.dot");
		int m(4);
		if(testulimit(300)){
			int k(detectk(input));
			if(k<=2*m){
				cout<<"k too low"<<endl;
			}
			else{
			createoutfile(input.c_str(),output.c_str(),k,m);
			}
		}else{
		cout<<"ulimit too low"<<endl;
		}
	}
	if(argc==3){
		string input(argv[1]);
		string output(argv[2]);
		int m(4);
		if(testulimit(300)){
			int k(detectk(input));
			if(k<=2*m){
				cout<<"k too low"<<endl;
			}
			else{
			createoutfile(input.c_str(),output.c_str(),k,m);
			}
		}else{
		cout<<"ulimit too low"<<endl;
		}
	}
	if(argc==4){
		string input(argv[1]);
		string output(argv[2]);
		int m(atoi(argv[3])/2);
		if(testulimit(pow(4,m)+50))
		{
			int k(detectk(input));
			if(k<=2*m){
				cout<<"k too low"<<endl;
			}
			else{
			createoutfile(input.c_str(),output.c_str(),k,m);
			}
		}else{
		cout<<"ulimit too low"<<endl;
		}
	}


	return sys;
}
