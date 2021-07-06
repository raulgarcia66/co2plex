/* 
Minty's Generalized Algorithm for the Max Weight Co2Plex on a {claw,bull}-free Graph
*/

//#ifndef KGRAPH_H
//#define KGRAPH_H
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <string>
//#include <map>
//#include <ctime>

#include "modcppstd.hh"

//extern double copy_time;
#ifdef __cplusplus
extern "C" {
#endif

	struct CCdatagroup;
	int perfect_match(int ncount, CCdatagroup *dat, int ecount,
					int **elist, int **elen, int **ourMatch, int **ourWeights, char *blo_filename,
					char *mat_filename, int just_fractional, int no_fractional,
					int use_all_trees, int partialprice,
					double *totalzeit) ;

	
#ifdef __cplusplus
}        /* end extern "C" */
#endif

#ifdef CC_PROTOTYPE_ANSI
static void
    usage (char *name);
static int
    parseargs (int ac, char **av);

#else

static void
    usage ();
static int
    parseargs ();

#endif

static int seed = 0;
static char *blossom_file = (char *) NULL;
static char *match_file = (char *) NULL;
static char *datfilename = (char *) NULL;
static char *edgefilename = (char *) NULL;
static char *edgegenfname = (char *) NULL;
static int no_frac = 0;
static int just_frac = 0;
static int use_all_trees = 0;
static int binary_in = 0;
static int tsplib_in = 0;
static int nnodes_want = 0;
static int partialprice = 0;
//static int norm = CC_EUCLIDEAN;
static int no_price = 0;

void PrintVector(const vector<int> arr) {
	for (int i = 0; i < arr.size(); i++) {
		cout<<arr.at(i)<<" ";
	}
	cout<<endl;
}

void PrintVectorVector(const vector< vector<int> > arr) {
	
    for (int i = 0; i <arr.size(); i++) {
        for (int j = 0; j < arr[i].size(); j++) {
            cerr<<arr[i].at(j)<<" ";
        }
        cerr<<endl;
    }
    
}

void PrintMatrix(int** M,int NumberOfVertices) {
	
	for(int i = 0; i < NumberOfVertices; i++){
		cerr<<endl;//to get in the form of an nxn matrix
		for(int j = 0; j < NumberOfVertices; j++)
			cerr<<M[i][j]<<" ";
	}
	cerr<<endl;
	
}

int* copyArray(const int* data, int dsize) {
	int* copy = new int[dsize];
	memcpy(copy, data, dsize * sizeof(int));
	return copy;
}


void testing( int& playing, vector<int>& play) {
	playing = 100;
	
	play.push_back(12345);
}

bool equalFunction (int i, int j) {
    return (i==j);
}

class Graph{
	int n,m;
	int **C;
    vector<int> weights;
    vector<int> nodes;
    
	
public:
	Graph(char* InputFile, char* weightFile);
	~Graph();
	void Minty(vector<int> S);
    void initMinty();
    int numMatchingProblem;
    int numberOfBlackNeighbors(int vertex, vector<int> blackVertices);
    int weightOfaPath(vector<int> blackPath, vector<int> whitePath);
    bool violateCo2PlexProp( int vertex, vector<int> blackVertices);
    bool co2PlexAfterSymmDiff( int vertex, vector<int> blackVertices);
	bool adjacent(vector<int> A, vector<int> B);
    
	vector<int> edmondsCorrection1(vector<int> regNeib, vector<int> whiteWings, vector<int> irreg, vector< vector<int> > irregPairs); 
	int edmondsCorrection3(int y11, int yl1, int y12, int yl2, vector<int> whiteWings, vector<int> irreg, vector< vector<int> > irregPairs, vector<int> &k, int &wP11i, int &wP22i, int &wP11j, int &wP22j);
	vector<bool> iwapBF(vector<int> &weight, vector<int> &whiteReachIrreg, vector<int> whiteSource, vector<int> whiteSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices, vector< vector<int> > irregularPairs);
	vector<int> bellmanFordVariant(vector<int> regularSource, vector< vector<int> > regularSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices, vector< vector<int> > irregularPairs);
    vector<int> maxWeightWhiteAugPath(vector<int> SFT1, vector<int> SFT2, vector<int> FT1, vector<int> FT2, vector<int> blackVertices, vector<int> BT1, vector<int> BT2, vector<int> BT3);
    vector< vector<int> > Edmonds(int freeVertexA, int freeVertexB, vector<int> xa, vector<int> xb, vector<int> blackSetSingle, vector<int> blackSetPairs, vector<int> boundedVerticesT1, vector<int> boundedVerticesT2, vector<int>boundedVerticesT3, int &N);
    
};

//CODE

int main(int argc, char* argv[]){
    //cout << argv[1] << argv[2];
    if(argc != 3) {
        cerr<<" USAGE: ./executable edges.txt weights.txt"<<endl;
        exit(0);
    }
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
		g.initMinty();
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