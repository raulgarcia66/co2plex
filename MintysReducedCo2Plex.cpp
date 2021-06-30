/*
Minty's Algorithm on a Reduced Graph for the Max Weight Co2Plex on a {claw, bull}-free Graph
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

	struct CCdatagroup;  //might need file that this struct
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

void PrintMatrix(int** M,int NumberOfVertices) {	
	for(int i = 0; i < NumberOfVertices; i++){
		cerr<<endl;//to get in the form of an nxn matrix
		for(int j = 0; j < NumberOfVertices; j++)
			cerr<<M[i][j]<<" ";
	}
	cerr<<endl;
}

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

class Graph
{
    int n; // number of nodes
    int m; // number of edges
    int **C; //adjacency matrix
    vector<int> nodes; // list of nodes; indexing of nodes starts at 0
    vector<int> weights; // list of weights associated with each node
    //vector<int> *adj; //store edges in a vector where components are lists of integers

    public:
    Graph(char* inputFile, char* weightFile); //constructor
    ~Graph();
    void initMinty();
    void Minty(vector<int> S);
    int numMatchingProblem;
    bool adjacent(vector<int> A, vector<int> B);
    int numberOfBlackNeighbors(int vertex, vector<int> blackVertices);
    int weightOfaPath(vector<int> blackPath, vector<int> whitePath);

    vector<int> maxWeightWhiteAugPath(vector<int> SF, vector<int> F2, vector<int> B, vector<int> blackVertices);
	vector< vector<int> > Edmonds(int freeVertexA, int freeVertexB, int xa, int xb, vector<int> blackVertices, vector<int> boundedVertices, int &N);
    vector<int> edmondsCorrection1(vector<int> regNeib, vector<int> whiteWings, vector<int> irreg); 
	int edmondsCorrection3(int y11, int yl1, int y12, int yl2, vector<int> whiteWings, vector<int> irreg, vector< vector<int> > irregPairs, vector<int> &k, int &wP11i, int &wP22i, int &wP11j, int &wP22j);
	vector<bool> iwapBF(vector<int> &weight, vector<int> &whiteReachIrreg, vector<int> whiteSource, vector<int> whiteSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices);
	//vector<int> bellmanFordVariant(vector<int> regularSource, vector< vector<int> > regularSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices, vector< vector<int> > irregularPairs);
    
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
		//G = createExtGraph(g); // creates the extended graph we apply Minty's to
        f_in.close();
        f_win.close();
		g.initMinty();
        //G.initMinty();
		// Then need to extract vertices of g from G
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