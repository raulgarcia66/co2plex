
// Algorithm 1 for max weight Co2Plex problem
//#include "Graph.cpp"
#include <iostream>
#include <fstream>
//#include <sstream>
#include <string>
#include <cstring>
//#include <ctime>
#include <vector>
#include <cmath>
//#include <cstdlib>
//#include <cstdio>

using namespace std;

int main(int argc, char* argv[]){
    //cout << argv[1] << argv[2];
    if(argc != 3) {
        cerr<<" USAGE: ./executable graphAdjacency graphWeights"<<endl;
        exit(0);
    }
    //Graph g(argv[1],argv[2]);
    //clock_t t1,t2;
    //t1=clock();
    
    ifstream f_in, f_win;
    f_in.open(argv[1]);
    f_win.open(argv[2]);
    
    if(!f_in.is_open() && !f_win.is_open()){
        cerr<< "Can't open files\n";
    }
    else{
        Graph g(argv[1], argv[2]);
        f_in.close();
        f_win.close();
        //g.initMinty();
    }
    

    // ostringstream fname;
	// fname << "verticeCo2PlexResults.txt";
	// //fname << ".txt";
	// ofstream f_out;
	// f_out.open( fname.str().c_str(),ofstream::app);
    
    // t2=clock();
    // float diff (((float)t2-(float)t1)/CLOCKS_PER_SEC);
    
    
	// f_out<<"time for: "<< argv[1]  << " was: "<< diff <<" s"<<endl;
	
    
    // //cout<<"time for: "<< argv[1]  << " was: "<< diff <<" s"<< endl;

    return 0;
}