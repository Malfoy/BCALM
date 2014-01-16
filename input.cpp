#include <string>
#include <iostream>
#include <fstream>

#include "input.h"
#include "ograph.h"



using namespace std;




void InputDot::init_input(string namefile, const int k){
    ifstream in(namefile);
	in.seekg(0,ios_base::end);
	size = in.tellg();
    buffer_size = 10000*(k+2);
    nb_buffers = size/buffer_size;
	in.seekg(0,ios::beg);
    index = 0;
    if(nb_buffers==0)
    {
		buffer=readn(&in,size);
	}else{
		buffer=readn(&in,buffer_size);
	}
}

string InputDot::next_input(const int k){
    if (index >= buffer.size())
    {
        index ++;
        if(index >= nb_buffers*buffer_size)
			buffer=readn(&in,size-nb_buffers*buffer_size);
		else
			buffer=readn(&in,buffer_size);
    }
    if(buffer.size()>=index+k)
    {
		kmer=buffer.substr(index,k);
	}else{
		kmer="";
	}
	index += k+2;
    return kmer;
}

// l'idee sera d'avoir des classes ayant la meme intreface que InputDot pour gerer d'autres format d'entree
