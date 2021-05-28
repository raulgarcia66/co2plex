/*
    Algorithm 1 for finding a Max Weight Co2Plex
    Code borrowed from Cynthia Wood's 
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
#include "time.h"

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
    //vector< vector<int> > Edmonds(int freeVertexA, int freeVertexB, vector<int> xa, vector<int> xb, vector<int> blackSetSingle, vector<int> blackSetPairs, vector<int> boundedVerticesT1, vector<int> boundedVerticesT2, vector<int>boundedVerticesT3, int &N);
    vector<int> edmondsCorrection1(vector<int> regNeib, vector<int> whiteWings, vector<int> irreg, vector< vector<int> > irregPairs); 
	int edmondsCorrection3(int y11, int yl1, int y12, int yl2, vector<int> whiteWings, vector<int> irreg, vector< vector<int> > irregPairs, vector<int> &k, int &wP11i, int &wP22i, int &wP11j, int &wP22j);
	vector<bool> iwapBF(vector<int> &weight, vector<int> &whiteReachIrreg, vector<int> whiteSource, vector<int> whiteSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices, vector< vector<int> > irregularPairs);
	vector<int> bellmanFordVariant(vector<int> regularSource, vector< vector<int> > regularSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices, vector< vector<int> > irregularPairs);
    
    
};

Graph::Graph(char* inputFile, char* weightFile){
    ifstream f_in, f_win;
	f_in.open(inputFile);
	
	f_in>>n;//number of vertices
	f_in>>m;//number of edges
    //cout << "Number of Nodes: " << n << "\nNumber of Edges: " << m << "\n";

    int *x = new int[m];
	int *y = new int[m];

    for(int i = 0; i < m; i++){
		f_in>>x[i]>>y[i];  //Cynthias code stores a third element in a temp variable z
	}
	f_in.close();
    // for(int i = 0; i < m; i++){
	// 	cout << x[i] << " " << y[i] << "\n";
	// }

    int z,w;
    f_win.open(weightFile);
    for (int i = 0; i < n; i++) {
        f_win>>z>>w;
        //nodes.push_back(z); // nodes is meant to store the current white nodes
        weights.push_back(w);
    }
    // for(int i = 0; i < n; i++){
	// 	cout << nodes[i] << " " << weights[i] << "\n";
	// }
    
    f_win.close();
	
	//Create the attribute C matrix
	C = new int*[n];
	for( int i = 0; i < n; i++) {
		C[i] = new int[n];
	}
	
    // initialize values of C
	for(int i = 0; i < n; i++){
		for(int j = i; j < n; j++) {
			C[i][j] = 0;
			C[j][i] = 0;
		}
	}

    //create adjacency matrix
	for(int i = 0; i < m; i++){
		C[x[i]][y[i]] = 1;
		C[y[i]][x[i]] = 1;
        //cout << C[x[i]-1][y[i]-1] << " ";
	}
    PrintMatrix(C,n);

	delete[] x;
	delete[] y;

}

Graph::~Graph(){
	for(int i = 0; i < n; i++){
		delete[] C[i];
	}
	delete[] C;
}

void Graph::initMinty(){
    //1. S is originally empty
    vector<int> S;
   // Store all the nodes (all vertices are white on the first step)
    for (int i = 0; i < n; i++) {
        nodes.push_back(i);
        //S.push_back(i);
    }
    //PrintVector(S); // if you want to see the nodes pushed
    Minty(S);
}

bool Graph::adjacent(vector<int> A, vector<int> B){
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < B.size(); j++) {
			if (C[A.at(i)][B.at(j)]) {
				return true;
			}
		}
	}
	return false;
}	

int Graph::numberOfBlackNeighbors(int vertex, vector<int> blackVertices){
	// Can't have more than 2 black neighbors
    int numOfN = 0;

    for (int i = 0; i < blackVertices.size(); i++) {
        if ( C[vertex][blackVertices.at(i)]) {
            numOfN++;
        }
        if (numOfN == 3) {
            break;
        }
    }
    return numOfN;
}

int Graph::weightOfaPath(vector<int> weightOfBlackPath, vector<int> weightofWhitePath){ 
    
	int weightWhite = accumulate(weightofWhitePath.begin(),weightofWhitePath.end(),0); 
    int weightBlack = accumulate(weightOfBlackPath.begin(), weightOfBlackPath.end(), 0);   
    int weightPath = weightWhite - weightBlack;
    
    //cout<<"weightWhite: "<<weightWhite<<endl;
    //cout<<"weightBlack: "<<weightBlack<<endl;    
    return weightPath;
}

vector<int> Graph::maxWeightWhiteAugPath(vector<int> SF, vector<int> F, vector<int> B, vector<int> blackVertices){
	
	sort (blackVertices.begin(),blackVertices.end());
    // vector<int> mergeFT12(FT1.size()+FT2.size()+SFT2.size());//Contains FT1, FT2 & SFT2
    // vector<int> mergeF(FT1.size()+FT2.size());
    // merge (FT1.begin(),FT1.end(),FT2.begin(),FT2.end(),mergeF.begin());
    // merge (mergeF.begin(),mergeF.end(),SFT2.begin(),SFT2.end(),mergeFT12.begin());
    //cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
    //PrintVector(mergeFT12);

    //2.1 Generate all white alternating paths of length 0 or 2
    
    //Save SF vertex with maximum weight only
    int maxWeightl0 = -1;
	int maxWeightln = -1;
    vector<int> maxWWAPl0;
	vector<int> maxWWAPln;
    for (int i = 0; i < SF.size(); i++) {
        int tempMax = weights.at(SF.at(i));
        if (tempMax > maxWeightl0) {
            maxWWAPl0.clear();
            maxWWAPl0.push_back(SF.at(i));
            maxWeightl0 = tempMax;
        }
    }

	// Ignore this block as black pairs can't exist in a stable set --------------------------------------
    // //Now we will generate WAP of length 2 and save it only if its weight is greater than the current max weight
    // vector<int> blackPairs;//Will always be even. An odd entry will always be paired with the even entry that preceeds it
    // vector<int> sortBlackPairs;
    // for (int i = 0; i < blackVertices.size(); i++) {
    //     for (int j = i+1; j < blackVertices.size(); j++) {
    //         if (C[blackVertices.at(i)][blackVertices.at(j)]) {
    //             blackPairs.push_back(blackVertices.at(i));
    //             blackPairs.push_back(blackVertices.at(j));
    //             sortBlackPairs.push_back(blackVertices.at(i));
    //             sortBlackPairs.push_back(blackVertices.at(j));
    //             break;// Since it is only adjacent to one other black vertex
    //         }
    //     }
    // }
    
    // //Black pairs are not contained in white alternating paths of length 2
    // vector<int> blackSingle(n);                      // 0  0  0 ... 0 (n of them)
    // vector<int>::iterator it;
    
	// //sort (blackVertices.begin(), blackVertices.end());you sort them at the beginning of the function
    // sort (sortBlackPairs.begin(),sortBlackPairs.end());
    
    // it=std::set_symmetric_difference (blackVertices.begin(),blackVertices.end(), sortBlackPairs.begin(),sortBlackPairs.end(), blackSingle.begin());
    // //  blackvertices ... 0  0  0  0
    // blackSingle.resize(it-blackSingle.begin());//only black vertices
	//---------------------------------------------------------------------------------------------------


	//All black vertices are single, hence replaced blackSingle with blackVertices
	//Now we will generate WAP of length 2 and save it only if its weight is greater than the current max weight
    vector<int> weightWhite;
    vector<int> weightBlack;
    for (int i = 0; i < blackVertices.size(); i++) {
        int adjacent1 = -1;
        int adjacent2 = -1;
        weightBlack.push_back(weights.at(blackVertices.at(i)));
        //for (int j = 0; j < mergeFT12.size(); j++) {
		for (int j = 0; j < F.size(); j++) { 
            if (adjacent1 == -1 && C[blackVertices.at(i)][F.at(j)]) {
                adjacent1 = F.at(j);
                for (int k = j+1; k < F.size(); k++) {
                    if (C[blackVertices.at(i)][F.at(k)]  && !C[F.at(j)][F.at(k)] ) {
                        adjacent2 = F.at(k);
                        weightWhite.push_back(weights.at(F.at(j)));
                        weightWhite.push_back(weights.at(F.at(k)));
                        
                        int tempw = weightOfaPath(weightBlack, weightWhite);
                        if (tempw > maxWeightl0) {
                            maxWWAPl0.clear();
                            maxWWAPl0.push_back(F.at(j));
                            maxWWAPl0.push_back(blackVertices.at(i));
                            maxWWAPl0.push_back(F.at(k));
                            maxWeightl0 = tempw;
                        }
                        weightWhite.clear();
                    }
                    
                }
            }
        }
        weightBlack.clear();
    }
    
    //2.2 For each pair of nonadjacent free vertices a and b do
    
    for (int a = 0; a < F.size(); a++) {
        for (int b = a+1; b < F.size(); b++) {
            if (!C[F.at(a)][F.at(b)]) {
                //vector<int> xa;// black vertex or vertices adjacent to a // old code
                //vector<int> xb;// black vertex or vertices adjacent to b
				//Let xa and xb be the black vertices andjacent to a and b respectively (only 1 because free)
				int xa = -1;
				int xb = -1;
                for (int i = 0; i < blackVertices.size(); i++) {
                    if (C[a][blackVertices.at(i)] && xa == -1) {
                        //xa.push_back(blackVertices.at(i));
						xa = blackVertices.at(i);
                    }
                    
                    if (C[b][blackVertices.at(i) && xb == -1]) {
                        //xb.push_back(blackVertices.at(i));
						xb = blackVertices.at(i);
                    }
					if (xa != -1 && xb != -1) {
						break;
					}
                }
                
                bool xEqual = false;
                
				// obsolete
                // for (int i = 0; i < xa.size(); i++) {
                //     for (int j = 0; j < xb.size(); j++) {
                //         //Want to make sure the black vertices to which they are joined are distinct
                //         if (xa.at(i) == xb.at(j)) {
                //             xEqual = true;
                //             break;
                //         }
                //     }
                    
                //     if (xEqual) { //needed?
                //         break;
                //     }
                // }

				if (xa == xb) {
					xEqual = true;
				}
                
                //if xa != xb
                if (!xEqual) {
				
					vector< vector<int> > EdmondsG;
                    //find MWWAP
					//N = 2*rbsN+2
					int N = 0;
                    //EdmondsG = Edmonds(F.at(a), F.at(b), xa, xb, blackSingle, blackPairs, BT1, BT2, BT3,N);
					EdmondsG = Edmonds(F.at(a), F.at(b), xa, xb, blackVertices, B, N);
					
					sort(EdmondsG.begin(), EdmondsG.end());
                    
                   // cout<<"LINE 460\n";
					
					if ( N > EdmondsG.size()) {
						double matzeit = 0.0;
						double genzeit = 0.0;
						int ncount, ecount;
						//long l; // Cynthia's code
						time_t l; // my fix
						
						//include the right header for this one otherwise it gives you a warning( "time.h")
						seed = time (&l);
						
						
						ncount = N;
						//mydat = NULL;
						ecount = EdmondsG.size();
						
						
						//cout << "Edmonds number of vertices: "<< N<<endl;
						//cout << "Edmonds number of edges: "<< EdmondsG.size()<<endl;
						//PrintVectorVector(EdmondsG);
						
						
						int *elen  = (int*)malloc(sizeof(int) * ecount);
						int *ourMatch = (int*)malloc(sizeof(int) * ncount);
						int *ourWeights = (int*)malloc(sizeof(int) * ncount/2);
						
						for (int i = 0; i < ecount; i++) {
                                //Negative since we are solving the Max weighted co-2-plex and blossom solves
                                //for min weight perfect matching
								elen[i] = -EdmondsG[i].at(2);

						}
		
						//elist = CC_SAFE_MALLOC (2*ecount, int);
						//int *elist = new int [2*ecount];
						int *elist = (int*)malloc(sizeof(int) * 2 * ecount);
						
						int idx = 0;
						for (int i = 0; i < ecount; i++) {
							elist[idx++] = EdmondsG[i].at(0);
							elist[idx++] = EdmondsG[i].at(1);
						}
						
						/*cout << "elist:"<<endl;  
						for (int i = 0; i < 2*ecount; i++) {
							cout << elist[i]<<" ";
						}
						cout << "\n";
						
						cout << "elen:"<<endl;
						for (int i = 0; i < ecount; i++) {
							cout << elen[i]<<" ";
						}
						cout << "\n";*/
						
                        //cout<<"LINE 515\n";
                        
						/*Returns 0 if it worked and 1 otherwise (for example, when one
						of the mallocs failed). The nodes in the graph should be named 
						0 through #nodes - 1 */
						if (perfect_match (ncount, NULL, ecount, &elist, &elen, &ourMatch, &ourWeights,
						   blossom_file, match_file, just_frac, no_frac, 
						   use_all_trees, partialprice, &matzeit)) {
							fprintf (stderr, "perfect_match failed\n");
			
						}
						
						//The mapping is a pair entry(i)-> entries 2i and 2i+1 weight
						//Entries 2i and 2i+1 correspond to node(i) in the original graph
						for (int i = 0; i < ncount; i++) {
							if (ourMatch[i]%2 == 0) {
								maxWWAPln. push_back(ourMatch[i]/2);
							
							}
							else {
								maxWWAPln. push_back((ourMatch[i]-1)/2);
							}

						}
						
						int tempMaxln = 0;
						for (int i = 0; i < maxWWAPln. size(); i++) {
							tempMaxln += weights. at(maxWWAPln. at(i));
						}
						
						if (tempMaxln > 0 && tempMaxln > maxWeightln) {
							maxWeightln = tempMaxln;
						}
						
						
						
						fflush (stdout);
						
						free(ourMatch);
						free(ourWeights);
					}
					
                    
					
                }
                
            }
        }
    }
    
	if (maxWWAPln > maxWWAPl0) {  //checks element by element as needed (both could have same weight)
		return maxWWAPln;
	}
	
    return maxWWAPl0;
}

vector < vector<int> > Graph::Edmonds(int freeVertexA, int freeVertexB, int xa, int xb, vector<int> blackVertices, vector<int> boundedVertices, int &N){
    
    vector< vector<int> > edmondsG;
   
    // sort(xa.begin(),xa.end());
    // sort(xb.begin(),xb.end());
    //sort(blackSetSingle.begin(),blackSetSingle.end());//Don't really need to do this anymore, but leave it for now
	//int weightxA = 0;
	//int weightxB = 0;
	
	//for (int i = 0 ; i < xa.size(); i++) {
	int weightxA = weights.at(xa);
	//}
	
	//for (int i = 0 ; i < xb.size(); i++) {
	int weightxB = weights.at(xb);
	//}
	
    /*We want to construct the reduced basic structure
     1.1 Ignore all the SF vertices (just don't pass them to the function)
     1.2 Ignore free vertices except a and b (not passed either)
     1.3 Ignore all white vertices adjacent to a or b, since they will never appear in an alternating path; 
	 	 Note: These are necessarily bounded vertices
     */
    
    // vector<int> whiteverticesNonAdjacentAorB( boundedVerticesT1.size() + boundedVerticesT2.size() + boundedVerticesT3.size());
    // vector<int> whitev2(boundedVerticesT1.size()+ boundedVerticesT2.size());
    // merge (boundedVerticesT1.begin(),boundedVerticesT1.end(),boundedVerticesT2.begin(),boundedVerticesT2.end(),whitev2.begin());
    // merge (whitev2.begin(),whitev2.end(),boundedVerticesT3.begin(),boundedVerticesT3.end(),whiteverticesNonAdjacentAorB.begin());
	vector <int> whiteverticesNonAdjacentAorB;

	whiteverticesNonAdjacentAorB.resize(boundedVertices.size(),'\0'); // does not contain A or B
	copy(boundedVertices.begin(),boundedVertices.end(),whiteverticesNonAdjacentAorB.begin());
    sort(whiteverticesNonAdjacentAorB.begin(), whiteverticesNonAdjacentAorB.end());
    
    //cout<<"#########################################################"<<endl;
    //PrintVector(whiteverticesNonAdjacentAorB);
   
   // remove white vertices adjacent to A or B; A and B are not included either above
    for (int i = 0; i < whiteverticesNonAdjacentAorB.size(); i++) {
		
        if (C[whiteverticesNonAdjacentAorB.at(i)][freeVertexA] || C[whiteverticesNonAdjacentAorB.at(i)][freeVertexB]) {
            whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + i);
            //When you remove you decrease the size of the vector
            i--;
        }
    }
    //cout<<"int A:"<< freeVertexA<<endl;
    //cout<<"int B:"<< freeVertexB<<endl;
    //cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
    //PrintVector(whiteverticesNonAdjacentAorB);
    
    //Now we form the reduced basic structure that is also the weight function. However, we dont have to get a copy of the weights of the nodes since we have global access to them (only for black vertices)
	vector<int> whiteRBS; //Contains the white vertices in the RBS
    
    whiteRBS.resize( whiteverticesNonAdjacentAorB.size(), '\0');
    copy(whiteverticesNonAdjacentAorB.begin(), whiteverticesNonAdjacentAorB.end(), whiteRBS.begin());
    whiteRBS.push_back(freeVertexA);
    whiteRBS.push_back(freeVertexB);
    
	vector<int> blackRBS; //Contains the black vertices in the RBS

	blackRBS.resize(blackVertices.size(),'\0');
	copy(blackVertices.begin(),blackVertices.end(),blackRBS.begin());

    // vector<int> blackSingleRBS;
    
    // blackSingleRBS.resize( blackSetSingle.size(), '\0');
    // copy(blackSetSingle.begin(), blackSetSingle.end(), blackSingleRBS.begin());
    
    // vector< vector<int> > blackPairsRBS;
    // vector<int> weightsBlackPairs;
    
    // for (int i = 0; i < blackPairsRBS.size(); i++) {
    //     //The original black set comes sorted, which implies entry[i] < entry[i+1]
    //     if (i%2 == 0) {
    //         vector<int> temp;
    //         temp.push_back(blackSetPairs.at(i));
    //         temp.push_back(blackSetPairs.at(i+1));
    //         sort(temp.begin(),temp.end());
            
    //         if (temp != xa && temp != xb) {
    //             blackPairsRBS.push_back(temp);
    //             int weight = 0;
    //             weight = weights.at(blackSetPairs.at(i)) + weights.at(blackSetPairs.at(i+1)) ;
    //             weightsBlackPairs.push_back(weight);
    //         }
    //         temp.clear();
    //     }
    // }
    
    
    /***********************CLASIFICATION OF BLACK VERTICES*************************/
    
    // A nonempty set of all bounded vertices, which are adjacent to the same two black objects x and y is called a wing
    // First figure out which black vertices and objects are adjacent to bounded vertices
    
    //int sinPairSize = blackSingleRBS.size() + blackPairsRBS.size();
    vector< int > numWings(blackRBS.size()); //keeps a counter to the number of wings each black vertex is adjacent to. They are 0 originally.
    
    vector< vector<int> > enumerateWings;//Keeps track of tip of the wings as pairs (ie, the two black vertices). Its size also gives how many different wings we have
	vector<int> whiteInWings; //Will store only one vertex that is in each distinct wing
	// I think enumerateWings and whiteInWings are the same size, and the i-th component of whiteInWings
	// is one vertex that is in the wing that is the i-th component of enumerateWings
	
    for (int i = 0; i < blackRBS.size(); i++) {
        int numi = numWings.at(i);
        for (int j = i+1; j < blackRBS.size(); j++) {
            
            int numj = numWings.at(j);
            for (int k = 0; k < whiteverticesNonAdjacentAorB.size(); k++) {
                if (C[whiteverticesNonAdjacentAorB.at(k)][blackRBS.at(j)] && C[whiteverticesNonAdjacentAorB.at(k)][blackRBS.at(i)]) {
                    numi++;
                    numj++;
					
                    whiteInWings.push_back(whiteverticesNonAdjacentAorB.at(k));
                    vector<int> tempWing;
                    tempWing.push_back(blackRBS.at(i));
                    tempWing.push_back(blackRBS.at(j));
                    sort(tempWing.begin(), tempWing.end());
                    enumerateWings.push_back(tempWing);
                    tempWing.clear();
                    break;//as long as the set = wing is nonempty then we know a particular vertex is part of a wing
    
                }
                
              
            }
            numWings.at(j) = numj;
        }
        
        numWings.at(i) = numi;
    }
    
    // for (int i = 0; i < blackSingleRBS.size(); i++) {
    //     int numi = numWings.at(i);
    //     for (int j = 0; j < blackPairsRBS.size(); j++) {
            
    //         int numj = numWings.at(blackSingleRBS.size() + j);
    //         for (int k = 0; k < whiteverticesNonAdjacentAorB.size(); k++) {
    //             if (C[whiteverticesNonAdjacentAorB.at(k)][blackSingleRBS.at(j)] && (C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[j].at(0)] || C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[j].at(1)])) {
    //                 numi++;
    //                 numj++;
	// 				whiteInWings.push_back(whiteverticesNonAdjacentAorB.at(k));
    //                 vector<int> tempWing1;
    //                 tempWing1.push_back(blackSingleRBS.at(i));
    //                 sort(blackPairsRBS[j].begin(), blackPairsRBS[j].end());
    //                 tempWing1.push_back(blackPairsRBS[j].at(0));
    //                 tempWing1.push_back(blackPairsRBS[j].at(1));
    //                 enumerateWings.push_back(tempWing1);
    //                 tempWing1.clear();
                    
                    
                    
    //                 break;//as long as the set = wing is nonempty then we know a particular vertex is part of a wing
    //             }
    //         }
    //         numWings.at(blackSingleRBS.size() + j) = numj;
    //     }
        
    //     numWings.at(i) = numi;
    // }
    
    
    // for (int i = 0; i < blackPairsRBS.size(); i++) {
    //     int numi = numWings.at(blackSingleRBS.size() + i);
	// 	sort(blackPairsRBS[i].begin(), blackPairsRBS[i].end());
    //     for (int j = i+1; j < blackPairsRBS.size(); j++) {
            
    //         int numj = numWings.at(blackSingleRBS.size() + j);
    //         for (int k = 0; k < whiteverticesNonAdjacentAorB.size(); k++) {
    //             if ((C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[i].at(0)] || C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[i].at(1)]) && (C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[j].at(0)] || C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[j].at(1)])) {
    //                 numi++;
    //                 numj++;
					
	// 				whiteInWings.push_back(whiteverticesNonAdjacentAorB.at(k));
    //                 vector<int> tempWing2;
                    
	// 				sort(blackPairsRBS[j].begin(), blackPairsRBS[j].end());
    //                 tempWing2.push_back(blackPairsRBS[i].at(0));
    //                 tempWing2.push_back(blackPairsRBS[i].at(1));
    //                 tempWing2.push_back(blackPairsRBS[j].at(0));
    //                 tempWing2.push_back(blackPairsRBS[j].at(1));
    //                 enumerateWings.push_back(tempWing2);
    //                 tempWing2.clear();
    //                 break;//as long as the set = wing is nonempty then we know a particular vertex is part of a wing
    //             }
    //         }
    //         numWings.at(blackSingleRBS.size() + j) = numj;
    //     }
        
    //     numWings.at(blackSingleRBS.size() + i)= numi;
    // }
    
	sort(enumerateWings.begin(), enumerateWings.end());
	sort(whiteInWings.begin(), whiteInWings.end());
	vector<int>::iterator itWiW;
	itWiW = unique( whiteInWings.begin(), whiteInWings.end());
	whiteInWings.resize( std::distance(whiteInWings.begin(),itWiW)); // gets rid of repeats
	
	
    //cout<<"***************************************************************"<<endl;
    //PrintVector(numWings);
    
    //Regular I: vectors xa and xb. Make sure you don't count this in black pairs
			// I'm guessing they may also be classified as regularII or irregular
    //Regular II: A black vertex adjacent to 3 or more wings
    vector<int> regularII; //vector< vector<int> > regularIIPairs;
    //Irregular: Adjacent to exactly 2 wings
    //useless otherwise
    vector<int> irregular; //vector< vector<int> > irregularPairs; int irregularPairsSize =0;
	
    for (int i = 0; i < blackRBS.size(); i++) {
        if (numWings.at(i) > 2) { // was 3 for some reason
			//cout << "regular"<<endl;
            regularII.push_back(blackRBS.at(i));
        }
        if (numWings.at(i) == 2) {
			//cout << "irregular"<<endl;
            irregular.push_back(blackRBS.at(i));
        }
        
        //A useless vertex along with all adjoining white vertices can be deleted
		//Not sure why. Keep in code for now, might manifest itself later
        if (numWings.at(i) <2) {
            
            for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
                if (C[blackRBS.at(i)][whiteverticesNonAdjacentAorB.at(j)]) {
				
					vector<int>::iterator itw;
					itw = find (whiteInWings.begin(), whiteInWings.end(), whiteverticesNonAdjacentAorB.at(j));
            
					if (itw != whiteInWings.end()) {
						whiteInWings.erase(itw);
					}
					
					whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + j);
                    j--;
	
					
                }
            }
            blackRBS.erase(blackRBS.begin() + i);
            numWings.erase(numWings.begin() + i);
            i--;
        }
    }
    
    // int bsp = blackSingleRBS.size() + blackPairsRBS.size();
    // for (int i = blackSingleRBS.size(); i < bsp; i++) {
        
    //     if (numWings.at(i) > 3) {
	// 		cout << "regular pair"<<endl;
    //         regularIIPairs.push_back(blackPairsRBS[i - blackSingleRBS.size()]);
            
    //     }
    //     if (numWings.at(i) == 2) {
	// 		cout << "irregular pair"<<endl;
    //         irregularPairs.push_back(blackPairsRBS[i - blackSingleRBS.size()]);
	// 		irregularPairsSize++;
    //     }
        
    //     //A useless object along with all adjoining white vertices can be deleted
    //     // don't include them
    //     if (numWings.at(i) <2) {
            
    //         for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
    //             if (C[blackPairsRBS[i - blackSingleRBS.size()].at(0)][whiteverticesNonAdjacentAorB.at(j)] || C[blackPairsRBS[i - blackSingleRBS.size()].at(1)][whiteverticesNonAdjacentAorB.at(j)]) {
    //                 vector<int>::iterator itw2;
	// 				itw2 = find (whiteInWings.begin(), whiteInWings.end(), whiteverticesNonAdjacentAorB.at(j));
            
	// 				if (itw2 != whiteInWings.end()) {
	// 					whiteInWings.erase(itw2);
	// 				}
	
	// 				whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + j);
    //                 j--;
    //             }
    //         }
            
    //         //just keep them empty here, then you don't have to decrease i, change this to remove it for real 10.15.14
    //         blackPairsRBS[i - blackSingleRBS.size()].clear();
    //         blackPairsRBS.erase(blackPairsRBS.begin() + i - blackSingleRBS.size());
    //         i--;
    //     }
       
    // }
    
    //NEW  white RBS whiteverticesNonAdjacentAorB U freevertexA U freeVertexB instead of white RBS
    //N1(xa) = freevertexA and N1(xb) = freevertexB
    // for a regular vertex of the first kind put the free vertex into one node and the rest of its neighbors into the other class
    vector<int> n2xa;
    vector<int> n2xb;
    for(int i = 0 ; i < whiteverticesNonAdjacentAorB.size(); i++) {
		
		//vector<int> tempw;
		//tempw.push_back(whiteverticesNonAdjacentAorB.at(i));
        // if (xa.size() == 2) {
			
		// 	if (adjacent(xa, tempw)){
		// 		//cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
		// 		n2xa.push_back(whiteverticesNonAdjacentAorB.at(i));
		// 		//PrintVector(n2xa);
		// 	}
        // }
		
		// if (xb.size() == 2) {
		// 	//cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
		// 	if (adjacent(xb, tempw))n2xb.push_back(whiteverticesNonAdjacentAorB.at(i));
    
        // }
		
		//tempw.clear();
        //if (xa.size() == 1) {
        if (C[xa][whiteverticesNonAdjacentAorB.at(i)]) {
			n2xa.push_back(whiteverticesNonAdjacentAorB.at(i));    
        }
		//if (xb.size() == 1) {
        if (C[xb][whiteverticesNonAdjacentAorB.at(i)]) {
			n2xb.push_back(whiteverticesNonAdjacentAorB.at(i));  
        }
		
    }
    
    
    //For any regularII black vertex or pair v  or {v,w} get the neighborhood of v or {v,w}
    //int regIIsP = regularII.size() + regularIIPairs.size();
    //vector< vector<int> > neighborsV;//These are the neighbors of regularII vertices and pairs
	vector< vector<int> > neighborsV;//These are the neighbors of regularII vertices, all have to be white obv.
    vector<int> tempRegII;

    for (int i = 0;  i < regularII.size(); i++) {
        
        //if (i < regularII.size()) {
            
		for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
			if (C[regularII.at(i)][whiteverticesNonAdjacentAorB.at(j)]) {
				tempRegII.push_back(whiteverticesNonAdjacentAorB.at(j));
			}
		}
		
		neighborsV.push_back(tempRegII); //neigborsV.at(i) then contains a vector of neighbors of regularII node i, which may empty
		tempRegII.clear();
           
        //}
        // else{
            
        //     for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
        //         if (C[regularIIPairs[i - regularII.size()].at(0)][whiteverticesNonAdjacentAorB.at(j)] || C[regularIIPairs[i - regularII.size()].at(1)][whiteverticesNonAdjacentAorB.at(j)]) {
        //             tempRegII.push_back(whiteverticesNonAdjacentAorB.at(j));
        //         }
        //     }
            
        //     neighborsV.push_back(tempRegII);
        //     tempRegII.clear();
            
        // }
        
    }
	
    
    //Fact 2.2: For any regularII vertex v or pair, N(v) is uniquely partitioned into N1(v) and N2(v) so that for any x and y in distinct wings
    //Now we want to find all the pairs x and y in N(v) that belong to distict wings

    /* 1. We already enumerated the wings (put the tip of the wings on enumeratedWings and the size gives the number of possible wings)
       2. Save in which wing each neighbor is contained */
    
	vector< vector <int> > whichWings(whiteverticesNonAdjacentAorB.size(),vector<int> (0));//gives index of the wing in enumeratedWings for nonadjacenttoAandB

	
    for (int i = 0;  i < enumerateWings.size(); i++) {
        
        vector<int> enum1;
        vector<int> enum0;
        
        if (enumerateWings[i].size() == 2) {
            enum0.push_back(enumerateWings[i].at(0));
            enum1.push_back(enumerateWings[i].at(1));
			
            
        }
        
        if (enumerateWings[i].size() == 3) {
            enum0.push_back(enumerateWings[i].at(0));
            enum1.push_back(enumerateWings[i].at(1));
            enum1.push_back(enumerateWings[i].at(2));
            
        }
        
        if (enumerateWings[i].size() == 4) {
            enum0.push_back(enumerateWings[i].at(0));
            enum0.push_back(enumerateWings[i].at(1));
            enum1.push_back(enumerateWings[i].at(2));
            enum1.push_back(enumerateWings[i].at(3));

        }
        
        if (enum1.size() == 1 && enum0.size() == 1) {
            for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
                if (C[whiteverticesNonAdjacentAorB.at(j)][enum0.at(0)] && C[whiteverticesNonAdjacentAorB.at(j)][enum1.at(0)]) {
                    
                    whichWings[j].push_back(i);
     
                }
				
            }
			
            
        }
        
        if (enum1.size() == 2 && enum0.size() == 1) {
            for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); i++) {
                if (C[whiteverticesNonAdjacentAorB.at(j)][enum0.at(0)] &&
                    (C[whiteverticesNonAdjacentAorB.at(j)][enum1.at(0)] || C[whiteverticesNonAdjacentAorB.at(j)][enum1.at(1)])) {
                    whichWings[j].push_back(i);
                }
            }
 
        }
    
        if (enum1.size() == 2 && enum0.size() == 2) {
	
            for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); i++) {
                if ((C[whiteverticesNonAdjacentAorB.at(j)][enum0.at(0)] || C[whiteverticesNonAdjacentAorB.at(j)][enum0.at(1)]) && (C[whiteverticesNonAdjacentAorB.at(j)][enum1.at(0)] || C[whiteverticesNonAdjacentAorB.at(j)][enum1.at(1)])) {
                    whichWings[j].push_back(i);
                }
            }
			
        }
	  
    }
    
    //3. sort the vector you created in (2)
    //comment on 1/15
	for (int i = 0; i < whichWings.size(); i++) {
        sort(whichWings[i].begin(), whichWings[i].end());
    }
	
	/*cout << "whiteverticesNonAdjacentAorB size:"<< whiteverticesNonAdjacentAorB.size()<<endl;
	cout << "SECOND TIME whichWingsSssssssssssssssssssssssssssssssss"<<endl;
	PrintVectorVector(whichWings);
	cout<<endl;
    
	cout << "# of wings: "<<endl;;
	for (int i = 0 ; i < whichWings.size(); i++) {
		cout << "size:"<< whichWings[i].size()<<endl;
	}*/
	
	
    //4. Use symmetric difference and myvector.find() to see if two of the neighbor vertices are in different wings
    // Figure out if x and y are in the same wing, that is if they are adjacent to the same black vertex or pair other than V
    vector< vector<int> > nv1(regIIsP,vector<int> (0));
    vector< vector<int> > nv2(regIIsP,vector<int> (0));
	
    for (int i = 0;  i < neighborsV.size(); i++) {
        
        for (int j = 0 ; j < neighborsV[i].size(); j++) {
            for (int k = j+1; k < neighborsV[i].size(); k++){
			
				vector<int>::iterator it;
				vector<int> symmDiffWing(enumerateWings.size());                  // 0  0  0 ... 0 (enumerateWings.size() of them)
                
                
                int x = -1;
                int y = -1;
                for (int l = 0; l < whiteverticesNonAdjacentAorB.size(); l++) {
                    if (neighborsV[i].at(j) == whiteverticesNonAdjacentAorB.at(l)) {
                        x = l;  
				                
					}
					
					
                    if (neighborsV[i].at(k) == whiteverticesNonAdjacentAorB.at(l)) {
						
                        y = l;
                    }
					
					if (x != -1 && y != -1) {
						break;
					}
					
				}
				
				
				if (x!=-1 && y != -1) {
					//Recall whichWings[i] contains indices to the wings that whiteVerticesNonAorB.at(i) belongs to
				
					it=std::set_symmetric_difference (whichWings[x].begin(),whichWings[x].end(), whichWings[y].begin(),whichWings[y].end(), symmDiffWing.begin());
					
					
					//  whitevertices ... 0  0  0  0
					symmDiffWing.resize(it-symmDiffWing.begin());//only white vertices
					
					//cout << "symDiff2:"<<endl;
					//PrintVector(symmDiffWing);
					
					int ww = whichWings[x].size()+ whichWings[y].size();
					
					if (symmDiffWing.size() == ww) {
						//figure put if they are adjacent
						if (!C[whiteverticesNonAdjacentAorB.at(x)][whiteverticesNonAdjacentAorB.at(y)]) {
                        
							bool inNv1x = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(x)) == nv1[i].end();
							bool inNv2x = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(x)) == nv2[i].end();
							
							bool inNv1y = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(y)) == nv1[i].end();
							bool inNv2y = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(y)) == nv2[i].end();
							
							//Case 1. if whiteverticesNonAdjacentAorB.at(x) is n1v or n2v 
							// and whiteverticesNonAdjacentAorB.at(y) is not in the other
							if (inNv1x && !inNv2y) {
								
								nv2[i].push_back(whiteverticesNonAdjacentAorB.at(y));
								
							}
							
							if (inNv2x && !inNv1y) {
								
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(y));
								
							}
							
							//Case 2. if whiteverticesNonAdjacentAorB.at(y) is n1v or n2v 
							//and whiteverticesNonAdjacentAorB.at(x) is not in the other
							if (inNv1y && !inNv2x) {
								
								nv2[i].push_back(whiteverticesNonAdjacentAorB.at(x));
								
							}
							
							if (inNv2y && !inNv1x) {
								
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x));
								
							}
							//Case 3. If neither whiteverticesNonAdjacentAorB.at(x) 
							// or whiteverticesNonAdjacentAorB.at(y) are in n1v or n2v
							if (!inNv1x && !inNv2x && !inNv1y && !inNv2y ) {
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x));
								nv2[i].push_back(whiteverticesNonAdjacentAorB.at(y));
							}
                        }
						else{
						
							// if they are adjacent they can go in the same set of neighbors
							bool inNv1x2 = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(x)) == nv1[i].end();
							bool inNv2x2 = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(x)) == nv2[i].end();
							
							bool inNv1y2 = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(y)) == nv1[i].end();
							bool inNv2y2 = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(y)) == nv2[i].end();

							
							//Case 1. if whiteverticesNonAdjacentAorB.at(x) is n1v or n2v 
							// and whiteverticesNonAdjacentAorB.at(y) is not in the other
							if (inNv1x2 && !inNv2y2) {
								
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(y));
								
							}
							
							if (inNv2x2 && !inNv1y2) {
								
								nv2[i].push_back(whiteverticesNonAdjacentAorB.at(y));
								
							}
							
							//Case 2. if whiteverticesNonAdjacentAorB.at(y) is n1v or n2v 
							//and whiteverticesNonAdjacentAorB.at(x) is not in the other
							if (inNv1y2 && !inNv2x2) {
								
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x));
								
							}
							
							if (inNv2y2 && !inNv1x2) {
								
								nv2[i].push_back(whiteverticesNonAdjacentAorB.at(x));
								
							}
							//Case 3. If neither whiteverticesNonAdjacentAorB.at(x) 
							// or whiteverticesNonAdjacentAorB.at(y) are in n1v or n2v
							if (!inNv1x2 && !inNv2x2 && !inNv1y2 && !inNv2y2 ) {
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x));
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(y));
							}
							
							
						}
						//comment 1/30
					}
					else{
					
						//5. if they are in the same wing they can be in the same group of neighbors
						bool inNv1x3 = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(x)) == nv1[i].end();
						bool inNv2x3 = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(x)) == nv2[i].end();
						
						bool inNv1y3 = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(y)) == nv1[i].end();
						bool inNv2y3 = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(y)) == nv2[i].end();
						
						//Case 1. if whiteverticesNonAdjacentAorB.at(x) is n1v or n2v 
						// and whiteverticesNonAdjacentAorB.at(y) is not in the other
						if (inNv1x3 && !inNv2y3) {
							
							nv1[i].push_back(whiteverticesNonAdjacentAorB.at(y));
							
						}
						
						if (inNv2x3 && !inNv1y3) {
							
							nv2[i].push_back(whiteverticesNonAdjacentAorB.at(y));
							
						}
						
						//Case 2. if whiteverticesNonAdjacentAorB.at(y) is n1v or n2v 
						//and whiteverticesNonAdjacentAorB.at(x) is not in the other
						if (inNv1y3 && !inNv2x3) {
							
							nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x));
							
						}
						
						if (inNv2y3 && !inNv1x3) {
							
							nv2[i].push_back(whiteverticesNonAdjacentAorB.at(x));
							
						}
						//Case 3. If neither whiteverticesNonAdjacentAorB.at(x) 
						// or whiteverticesNonAdjacentAorB.at(y) are in n1v or n2v
						if (!inNv1x3 && !inNv2x3 && !inNv1y3 && !inNv2y3 ) {
							nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x));
							nv1[i].push_back(whiteverticesNonAdjacentAorB.at(y));
						}
                    
					}
					
				}
				
            }
        }
    }
    
	//cout << "rbs before correction:"<< regularII.size() + regularIIPairs.size() <<endl;
    // A regular vertex with empty node class, together with all adjoining white vertices, may be deleted
	for (int i = 0; i < neighborsV.size(); i++) {
        vector<int> tempn1;
        vector<int> tempn2;
        
        for (int j = 0; j < nv1[i].size(); j++) {
            tempn1.push_back(nv1[i].at(j));
        }
        
        for (int j = 0; j < nv2[i].size(); j++) {
            tempn2.push_back(nv2[i].at(j));
        }
        
        if (tempn1.size() == 0 || tempn2.size() == 0) {
            if (i < regularII.size()) {
                regularII.erase(regularII.begin() + i);
                
                for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
                    int du = whiteverticesNonAdjacentAorB.at(j);
                    
                    bool inN = find(neighborsV[i].begin(), neighborsV[i].end(), du) == neighborsV[i].end();
                    if (!inN) {
					
						vector<int>::iterator itw3;
						itw3 = find (whiteInWings.begin(), whiteInWings.end(), whiteverticesNonAdjacentAorB.at(j));
            
						if (itw3 != whiteInWings.end()) {
							whiteInWings.erase(itw3);
						}
                        whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + j);
                        j--;
                    }
                }
                
                neighborsV.erase(neighborsV.begin() + i);
                i--;
            }
            else{
                regularIIPairs.erase(regularIIPairs.begin() + i);
                
                for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
                    int du2 = whiteverticesNonAdjacentAorB.at(j);
                    
                    bool inN = find(neighborsV[i].begin(), neighborsV[i].end(), du2) == neighborsV[i].end();
                    if (!inN) {
						vector<int>::iterator itw4;
						itw4 = find (whiteInWings.begin(), whiteInWings.end(), whiteverticesNonAdjacentAorB.at(j));
            
						if (itw4 != whiteInWings.end()) {
							whiteInWings.erase(itw4);
						}
					
					
                        whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + j);
                        j--;
                    }
                }
                
                neighborsV.erase(neighborsV.begin() + i);
                i--;
            }
        }
        
        
        
        tempn1.clear();
        tempn2.clear();
        
    }
    
    //cout << "regIIsP4: "<< regularII.size() + regularIIPairs.size()<<endl;
    //Construction of the Edmonds' Graph
	//Assume the RBS has N regular vertices ( N = rbsN+2) for the two regular I vertices
	int rbsN = regularII.size() + regularIIPairs.size();
	N = 2*rbsN+2;
    //We form a graph with 2N+2 nodes and N black branches
    // New graph store it in a vectorofvectors in vertex vertex weight format
    //The maping is a pair entry(i)-> entries 2i and 2i+1 weight
	vector< vector<int> > newGraphB;
	
	
	//a^ and b^ will be nodes 0 and 1 respectively
	//x_a^1, x_a^2 and x_b^1, x_b^2 will be 2,3,4,5
	
	
	//cout << "rbsN: "<< rbsN <<endl;
	
	//These two are white edges
	vector<int> tempR;
	tempR.push_back(0);
	tempR.push_back(2);
	tempR.push_back(weights.at(freeVertexA));
	newGraphB.push_back(tempR);
	tempR.clear();
	tempR.push_back(1);
	tempR.push_back(4);
	tempR.push_back(weights.at(freeVertexB));
	newGraphB.push_back(tempR);
	tempR.clear();
	/////////////////Now for xa and xb, they are black edges
	tempR.push_back(2);
	tempR.push_back(3);
	tempR.push_back(weightxA);
	newGraphB.push_back(tempR);
	tempR.clear();
	tempR.push_back(4);
	tempR.push_back(5);
	tempR.push_back(weightxB);
	newGraphB.push_back(tempR);
	tempR.clear();
	
	//cout << "reg II single size:" << regularII.size()<<endl;
	//cout << "rbsN: "<< rbsN<<endl;
	//These are black edges
	for (int i = 3; i < rbsN+3; i++) {
        if (i < regularII.size()+3) {
            vector<int> tempE;
            tempE.push_back(2*i);
            tempE.push_back(2*i+1);
            tempE.push_back(weights.at(regularII.at(i-3)));
            newGraphB.push_back(tempE);
            tempE.clear();
            
        }
        else{
            //This takes care of the contraction of black pairs
            vector<int> tempE2;
            tempE2.push_back(2*i);
            tempE2.push_back(2*i+1);
            int sum1 = regularIIPairs[i-3].at(0);
            int sum2 = regularIIPairs[i-3].at(1);
            int w = weights.at(sum1) + weights.at(sum2);
            tempE2.push_back(w);
            newGraphB.push_back(tempE2);
            tempE2.clear();
            
        }
    }

	//Add first the WHITE edges associated with regular I vertices
	//N1(xa) = freevertexA and N1(xb) = frevertexB
    // for a regular vertex of the first kind put the free vertex into one node and the rest of its neighbors into the other class
    //n2xa and n2xb;
	//Here we only consider n2xa and n2xb, since there will not be IWAPS starting at xa or xb
	//iwapBF(vector<int> &weight, vector<int> whiteSource, vector<int> whiteSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices, vector< vector<int> > irregularPairs)
	
	int numberOfEdgesEdmonds = rbsN+2;
	bool IWAP = false;
	int maxW = 0;
	vector<int> wxa;
	vector<int> whiteReachI;
	for (int i = 0; i < n2xa.size(); i++) {
		vector<int> xaSource;
		xaSource.push_back(n2xa.at(i));
		vector<bool> IWAPBF = iwapBF(wxa,whiteReachI, xaSource, n2xb, whiteInWings, irregular, irregularPairs);
		int counter = 0;
		for (int j = 0; j < IWAPBF.size(); j++) {
			if (IWAPBF.at(j)) {
				IWAP = true;
				
				if (wxa.at(counter) > maxW) {
					maxW = wxa.at(counter);
					counter++;
				}
			}
		}
		
		xaSource.clear();
		wxa.clear();
	}
	
	if (IWAP) {
		vector<int> tempE3;
		//Recall node A is 0 and node B is 1 
		tempE3.push_back(0);
		tempE3.push_back(1);
		tempE3.push_back(maxW);
		newGraphB.push_back(tempE3);
		numberOfEdgesEdmonds++;
		tempE3.clear();
	}
	
	bool IWAP1 = false;
	bool IWAP2 = false;
	int maxW1 = 0;
	int maxW2 = 0;
	vector<int> wxa1, wxa2;
    for (int i = 0; i < n2xa.size(); i++) {
		vector<int> xaaSource;
		xaaSource.push_back(n2xa.at(i));
		for (int j = 0; j < nv1.size(); j++) {
			vector<bool> IWAPBF2 = iwapBF(wxa1,whiteReachI, xaaSource, nv1.at(j), whiteInWings, irregular, irregularPairs);
			int counter1 = 0;
			for (int k = 0; k < IWAPBF2.size(); k++) {
				if (IWAPBF2.at(k)) {
					IWAP1 = true;
					
					if (wxa.at(counter1) > maxW1) {
						maxW1 = wxa.at(counter1);
						counter1++;
					}
				}
			}
			
			if (IWAP1) {
				vector<int> tempE4;
				//Recall node A is 0 and regular vertices xj with N1 are mapped to vertex 2j
				tempE4.push_back(0);
				tempE4.push_back(2*(j+3));
				tempE4.push_back(maxW1);
				newGraphB.push_back(tempE4);
				numberOfEdgesEdmonds++;
				tempE4.clear();
				IWAP1 = false;
			}
			
			wxa1.clear();
		}
		

		for (int j = 0; j < nv2.size(); j++) {
			vector<bool> IWAPBF3 = iwapBF(wxa2, whiteReachI, xaaSource, nv2.at(j), whiteInWings, irregular, irregularPairs);
			int counter2 = 0;
			for (int k = 0; k < IWAPBF3.size(); k++) {
				if (IWAPBF3.at(k)) {
					IWAP2 = true;
					
					if (wxa.at(counter2) > maxW2) {
						maxW2 = wxa.at(counter2);
						counter2++;
					}
				}
			}
			
			if (IWAP2) {
				vector<int> tempE5;
				//Recall node A is 0 and regular vertices xj with N2 are mapped to vertex 2j+1
				tempE5.push_back(0);
				tempE5.push_back(2*(j+3)+1);
				tempE5.push_back(maxW2);
				newGraphB.push_back(tempE5);
				numberOfEdgesEdmonds++;
				tempE5.clear();
				IWAP2 = false;
			}
			wxa2.clear();
		}
		xaaSource.clear();
		
	}
	
	bool IWAP1b = false;
	bool IWAP2b = false;
	int maxW1b = 0;
	int maxW2b = 0;
	vector<int> wxb1, wxb2;
    for (int i = 0; i < n2xb.size(); i++) {
		vector<int> xbSource;
		xbSource.push_back(n2xb.at(i));
		for (int j = 0; j < nv1.size(); j++) {
			vector<bool> IWAPBF4 = iwapBF(wxb1,whiteReachI, xbSource, nv1.at(j), whiteInWings, irregular, irregularPairs);
			int counter1 = 0;
			for (int k = 0; k< IWAPBF4.size(); k++) {
				if (IWAPBF4.at(k)) {
					IWAP1b = true;
					
					if (wxb1.at(counter1) > maxW1b) {
						maxW1 = wxa.at(counter1);
						counter1++;
					}
				}
			}
			
			if (IWAP1b) {
				vector<int> tempE6;
				//Recall node B is 1 and regular vertices xj with N1 are mapped to vertex 2j
				tempE6.push_back(1);
				tempE6.push_back(2*(j+3));
				tempE6.push_back(maxW1b);
				newGraphB.push_back(tempE6);
				numberOfEdgesEdmonds++;
				tempE6.clear();
				IWAP1b = false;
			}
			
			wxb1.clear();
		}
	
		for (int j = 0; j < nv2.size(); j++) {
			vector<bool> IWAPBF5 = iwapBF(wxb2, whiteReachI, xbSource, nv2.at(j), whiteInWings, irregular, irregularPairs);
			int counter2 = 0;
			for (int k = 0; k < IWAPBF5.size(); k++) {
				if (IWAPBF5.at(k)) {
					IWAP2b = true;
					
					if (wxb2.at(counter2) > maxW2b) {
						maxW2b = wxb2.at(counter2);
						counter2++;
					}
				}
			}
			
			if (IWAP2b) {
				vector<int> tempE7;
				//Recall node B is 1 and regular vertices xj with N2 are mapped to vertex 2j+1
				tempE7.push_back(1);
				tempE7.push_back(2*(j+3)+1);
				tempE7.push_back(maxW2b);
				newGraphB.push_back(tempE7);
				numberOfEdgesEdmonds++;
				tempE7.clear();
				IWAP2b = false;
			}
			
			wxb2.clear();
		}
		
		xbSource.clear();
		
	}	
	
	//cout << "size of new graph B1: "<< newGraphB.size()<<endl;
	
    
	//Correct Edmonds
	vector<int> whiteInIrregularWings;
	for (int i = 0; i < whiteverticesNonAdjacentAorB.size(); i++) {
		int countwiW =0;
		for (int j = 0; j < irregular.size(); j++) {
			if (C[whiteverticesNonAdjacentAorB.at(i)][irregular.at(j)]) {
				countwiW++;	
			}
		}
		
		vector<int> tiiW;
		tiiW. push_back(whiteverticesNonAdjacentAorB.at(i));
		for (int j = 0; j < irregularPairs.size(); j++) {
			if (adjacent(tiiW, irregularPairs.at(j))) {
				countwiW++;
			}
		}
		tiiW.clear();
		
		if (countwiW >=2) {
			whiteInIrregularWings. push_back(whiteverticesNonAdjacentAorB.at(i));
		}
	}
	
	
	vector<int>::iterator itWiiW;
	itWiiW = unique(whiteInIrregularWings.begin(), whiteInIrregularWings.end());
	whiteInIrregularWings.resize(distance(whiteInIrregularWings.begin(),itWiiW) );
	
	sort(whiteInIrregularWings.begin(), whiteInIrregularWings.end());
	//Construct W(xi,xj):= Union of all wings reachable by irregular vertices to both xi and xj
	//So far we have white in wing, we want to save white in irregular wings
	
	vector<int> wnv;
	bool IWAPnv = false;
	vector<int> whiteReachIW;
	for (int i = 0; i < nv1.size(); i++) {
		int maxNV = 0;
		for (int j = 0 ; j < nv1[i].size(); j++) {
			vector<int> nv1Source;
			nv1Source.push_back(nv1[i].at(j));
			
			for (int k = 0; k < nv2.size(); k++) {
			
				if (i > k || i <k) {
					vector<bool> IWAPBF8 = iwapBF(wnv, whiteReachIW, nv1Source, nv2.at(k), whiteInWings, irregular, irregularPairs);
					int counternv = 0;
					for (int l = 0; l < IWAPBF8.size(); l++) {
						if (IWAPBF8.at(l)) {
							IWAPnv = true;
							
							if (wnv.at(counternv) > maxNV) {
								maxNV = wnv.at(counternv);
								counternv++;
							}
						}
					}
					
					if (IWAPnv) {
		
						vector<int> tempE8;
						//Recall vertices xi with N1 are mapped to vertex 2i
						//and xk with N2 are mapped to vertex 2k+1
						tempE8.push_back(2*(i+3));
						tempE8.push_back(2*(k+3)+1);
						sort(tempE8.begin(), tempE8.end());
						tempE8.push_back(maxNV);
						newGraphB.push_back(tempE8);
						numberOfEdgesEdmonds++;
						tempE8.clear();
						IWAPnv = false;
					}
					wnv.clear();
					whiteReachIW.clear();

				}
			}
			
		}
	}
	
	sort( newGraphB. begin(), newGraphB. end());
	
	
	int edmondsEdgesCorrected = numberOfEdgesEdmonds;
	//Make correction of Edmonds' graph according to Nakamura & Tamura
	//edmondsCorrection1( vector<int> regNeib, vector<int> whiteWings, vector<int> irreg, vector< vector<int> >irregPairs)
	//Lemma 3.2 correction
	int regCorr = regularII.size() + regularIIPairs.size();
	for (int i = 0; i < regCorr; i++) {
		vector<int> wingsReachThroughIrregxi = edmondsCorrection1( neighborsV.at(i), whiteInWings, irregular, irregularPairs);
		for (int j = i+1; j < regCorr-1; j++) {
			vector<int> wingsReachThroughIrregxj = edmondsCorrection1( neighborsV.at(j), whiteInWings, irregular, irregularPairs);
			vector<int> Wxixj(wingsReachThroughIrregxi.size()+wingsReachThroughIrregxj.size());
			//Get union of all wings that are reachable to both xi and xj
			vector<int>::iterator itInt;
			
			itInt =set_intersection (wingsReachThroughIrregxi.begin(), wingsReachThroughIrregxi.end(), wingsReachThroughIrregxj.begin(),wingsReachThroughIrregxj.end(), Wxixj.begin());
			
			Wxixj.resize(itInt - Wxixj.begin());
			///////////////////////////////////////////////////////////////////////////////////
			vector<int> n1eWxixj(nv1[j].size()+Wxixj.size());
			
			vector<int>::iterator itInt2;
			
			itInt2 = set_intersection( nv1[j].begin(), nv1[j].end(), Wxixj. begin(), Wxixj. end(), n1eWxixj . begin());
			
			n1eWxixj.resize(itInt2 - n1eWxixj.begin());
			
			if (n1eWxixj.size() == nv1[j].size() && nv1[j].size() > 0) {
				//delete xi1xj2 and xi2xj2 (2*(i+3),2*(j+3)+1) and (2*(i+3)+1,2*(j+3)+1)
				for (int k = 0; k < newGraphB.size(); k++) {
					if (newGraphB[k].at(0) == 2*(i+3) && newGraphB[k].at(1) == (2*(j+3)+1)) {
						newGraphB.erase(newGraphB.begin() + k);
						edmondsEdgesCorrected--;
						k--;
					}
					
					if (newGraphB[k].at(0) == (2*(i+3)+1) && newGraphB[k].at(1) == (2*(j+3)+1)) {
						newGraphB.erase(newGraphB.begin() + k);
						edmondsEdgesCorrected--;
						k--;
					}
				}
				
			}
			////////////////////////////////////////////////////////////////////////////////////
			vector<int> n2eWxixj(nv2[j].size()+Wxixj.size());
			
			vector<int>::iterator itInt3;
			
			itInt3 = set_intersection( nv2[j].begin(), nv2[j].end(), Wxixj. begin(), Wxixj. end(), n2eWxixj . begin());
			
			n2eWxixj.resize(itInt3 - n2eWxixj.begin());
			
			if (n2eWxixj.size() == nv2[j].size() && nv2[j].size() > 0) {
				//delete xi1xj1 and xi2xj1 (2*(i+3),2*(j+3)) and (2*(i+3)+1,2*(j+3))
				
				for (int k = 0; k < newGraphB.size(); k++) {
					if (newGraphB[k].at(0) == 2*(i+3) && newGraphB[k].at(1) == 2*(j+3)) {
						newGraphB.erase(newGraphB.begin() + k);
						edmondsEdgesCorrected--;
						k--;
					}
					
					if (newGraphB[k].at(0) == (2*(i+3)+1) && newGraphB[k].at(1) == 2*(j+3)) {
						newGraphB.erase(newGraphB.begin() + k);
						edmondsEdgesCorrected--;
						k--;
					}
				}
			}
			//////////////////////////////////////////////////////////////////////////////////////
			vector<int> n1xieWxixj(nv1[i].size()+Wxixj.size());
			
			vector<int>::iterator itInt4;
			
			itInt4 = set_intersection( nv1[i].begin(), nv1[i].end(), Wxixj. begin(), Wxixj. end(), n1xieWxixj . begin());
			
			n1xieWxixj.resize(itInt4 - n1xieWxixj.begin());
			
			if (n1xieWxixj.size() == nv1[i].size() && nv1[i].size() > 0) {
				//delete xi2xj1 and xi2xj2 (2*(i+3)+1,2*(j+3)) and (2*(i+3)+1,2*(j+3)+1)
				for (int k = 0; k < newGraphB.size(); k++) {
					if (newGraphB[k].at(0) == (2*(i+3)+1) && newGraphB[k].at(1) == 2*(j+3)) {
						newGraphB.erase(newGraphB.begin() + k);
						edmondsEdgesCorrected--;
						k--;
					}
					
					if (newGraphB[k].at(0) == (2*(i+3)+1) && newGraphB[k].at(1) == (2*(j+3)+1)) {
						newGraphB.erase(newGraphB.begin() + k);
						edmondsEdgesCorrected--;
						k--;
					}
				}
			}
			////////////////////////////////////////////////////////////////////////////////////
			vector<int> n2xieWxixj(nv2[i].size()+Wxixj.size());
			
			vector<int>::iterator itInt5;
			
			itInt5 = set_intersection( nv2[i].begin(), nv2[i].end(), Wxixj. begin(), Wxixj. end(), n2xieWxixj . begin());
			
			n2xieWxixj.resize(itInt5 - n1xieWxixj.begin());
			
			if (n2xieWxixj.size() == nv2[i].size() && nv2[i].size() > 0) {
				//delete xi1xj1 and xi1xj2 (2*(i+3),2*(j+3)) and (2*(i+3),2*(j+3)+1)
				for (int k = 0; k < newGraphB.size(); k++) {
					if (newGraphB[k].at(0) == 2*(i+3) && newGraphB[k].at(1) == 2*(j+3)) {
						newGraphB.erase(newGraphB.begin() + k);
						edmondsEdgesCorrected--;
						k--;
					}
					
					if (newGraphB[k].at(0) == 2*(i+3) && newGraphB[k].at(1) == (2*(j+3)+1)) {
						newGraphB.erase(newGraphB.begin() + k);
						edmondsEdgesCorrected--;
						k--;
					}
				}
			}
			
			wingsReachThroughIrregxj.clear();
			
			Wxixj.clear();
		}
		
		wingsReachThroughIrregxi.clear();
	}
	
	//Make correction for lemma 3.4
	
	//cout << "nv1.size:"<<nv1.size()<<endl;
	//cout << "nv2.size:"<<nv2.size()<<endl;
	// Recall size of nv1 and nv2 is the same
	int nv = nv1.size();
	for (int i = 0; i < nv; i++) {
		for (int j = 0; j <nv; j++) {
			for (int l = 0; l < nv1[i].size(); l++) {
				for (int ll = 0; ll < nv2[j].size(); ll++) {
					if (C[nv1[i].at(l)][nv2[j].at(ll)]) {
						//delete xi1xj1 and xi2xj2 (2*(i+3),2*(j+3)) and (2*(i+3)+1,2*(j+3)+1)
						for (int k = 0; k < newGraphB.size(); k++) {
							if (newGraphB[k].at(0) == 2*(i+3) && newGraphB[k].at(1) == 2*(j+3)) {
								newGraphB.erase(newGraphB.begin() + k);
								edmondsEdgesCorrected--;
								k--;
							}
							
							if (newGraphB[k].at(0) == (2*(i+3)+1) && newGraphB[k].at(1) == (2*(j+3)+1)) {
								newGraphB.erase(newGraphB.begin() + k);
								edmondsEdgesCorrected--;
								k--;
							}
						}
					}
				}
				
			}
		}
	}
	
	//Make correction from lemma 3.5
	// edmondsCorrection3(int y11, int yl1, int y12, int yl2, vector<int> whiteWings, vector<int> irreg, vector< vector<int> > irregPairs, zk, wP11i, wP22i, wP11j, wP22j)
	// The function above returns the size of paths y11,z1,y21,...,zl-1,yl1 & y22,z1,y22,...,zl-1,yl2
	for (int i = 0; i < nv; i++) {
		for (int j = 0; j <nv; j++) {
			for (int l = 0; l < nv1[i].size(); l++) {
				for (int l2 = l+1; l2 < nv1[i].size(); l2++) {
					for (int ll = 0; ll < nv2[j].size(); ll++) {
						for (int ll2 = ll+1; ll2 < nv2[j].size(); ll2++) {
							vector<int> zk;
							int wP11i = 0;
							int wP22i = 0; 
							int wP11j = 0; 
							int wP22j = 0;
							if (edmondsCorrection3(nv1[i].at(l), nv1[i].at(l2),nv2[j].at(ll), nv2[j].at(ll2), whiteInWings, irregular, irregularPairs, zk, wP11i, wP22i, wP11j, wP22j) > 1) {
								//delete the four edges xi1xj1, xi2xj2, xi1xj2 and xi2xj1
								//(2*(i+3),2*(j+3)), (2*(i+3)+1,2*(j+3)+1), (2*(i+3),2*(j+3)+1) and (2*(i+3)+1,2*(j+3))
								for (int k = 0; k < newGraphB.size(); k++) {
									if ((newGraphB[k].at(0) == 2*(i+3) && newGraphB[k].at(1) == 2*(j+3)) ||
										(newGraphB[k].at(0) == 2*(i+3)+1 && newGraphB[k].at(1) == 2*(j+3)+1) ||
										(newGraphB[k].at(0) == 2*(i+3) && newGraphB[k].at(1) == 2*(j+3)+1) ||
										(newGraphB[k].at(0) == 2*(i+3)+1 && newGraphB[k].at(1) == 2*(j+3))
									) {
										newGraphB.erase(newGraphB.begin() + k);
										edmondsEdgesCorrected--;
										k--;
							
									}
								}
								
								//add to new vertices z(kk)i and z(kk)j(kk =k on Nakamura's paper)
								//The maping is a pair entry(i)-> entries 2i and 2i+1 weight
								vector <int> tempzk;
								if (zk. size() == 1) {
									tempzk. push_back(2*zk. at(0));
									tempzk. push_back(2*zk. at(0)+1);
									tempzk. push_back(weights. at(zk. at(0)));
								}
								else if (zk. size() == 2) {
								
									tempzk. push_back(2*zk. at(0));
									tempzk. push_back(2*zk. at(0)+1);
									tempzk. push_back(weights. at(zk. at(0)) + weights. at(zk. at(1)));
								
								}
								
								if (tempzk. size() > 0) {
									newGraphB. push_back(tempzk);
									tempzk. clear();
								}
								
								//add four white edges xi1-zki, xi2-zki, xj1-zkj, xj2-zkj
								if (zk. size() > 0 && wP11i >= 0 && wP22i >= 0 && wP11j >= 0 && wP22j >= 0) {
									tempzk. push_back(nv1[i].at(l));
									tempzk. push_back(2*zk. at(0));
									sort(tempzk.begin(), tempzk.end());
									tempzk. push_back(wP11i);
									newGraphB. push_back(tempzk);
									tempzk. clear();
									
									
									tempzk. push_back(nv1[i].at(l2));
									tempzk. push_back(2*zk. at(0));
									sort(tempzk.begin(), tempzk.end());
									tempzk. push_back(wP22i);
									newGraphB. push_back(tempzk);
									tempzk. clear();
									
									tempzk. push_back(nv1[i].at(ll));
									tempzk. push_back(2*zk. at(0)+1);
									sort(tempzk.begin(), tempzk.end());
									tempzk. push_back(wP11j);
									newGraphB. push_back(tempzk);
									tempzk. clear();
									
									
									tempzk. push_back(nv1[i].at(ll2));
									tempzk. push_back(2*zk. at(0)+1);
									sort(tempzk.begin(), tempzk.end());
									tempzk. push_back(wP22j);
									newGraphB. push_back(tempzk);
									tempzk. clear();
								}
								
							}
						}
					}
				}
				
			}
		}
	}

	//cout << "size of new graph B2: "<< newGraphB.size()<<endl;
	edmondsG = newGraphB;
  
    return edmondsG;   
}

void Graph::Minty(vector<int> S){
    
    //This recursion happens while S is not optimal. Later we will determine optimality of S
    //1. S = empty
    vector<int> myS;// We will make a copy of S, which is originally empty
    
    myS.resize( S.size(), '\0');  // '\0' is initial arbitrary value for elements
    copy(S.begin(), S.end(), myS.begin());
    
    //1.1 Obtain white vertices
    vector<int> whiteVertices(n);                      // 0  0  0 ... 0 (n of them)
    vector<int>::iterator it;
    
    sort (myS.begin(),myS.end());
    sort (nodes.begin(),nodes.end());   // Sort for consistency
    
    it=std::set_symmetric_difference (myS.begin(),myS.end(), nodes.begin(),nodes.end(), whiteVertices.begin());
    //  whitevertices ... 0  0  0  0
    whiteVertices.resize(it-whiteVertices.begin());//the white vertices
    
    cout << "The symmetric difference of S and nodes has " << (whiteVertices.size()) << " elements:\n";
    for (it=whiteVertices.begin(); it!=whiteVertices.end(); ++it)
       cout << ' ' << *it;
    cout << '\n';
    
    //2. Clasify the white vertices based off the black vertices
    
    // We have superfree, free and bounded vertices
    vector<int> superFreeVertices, freeVertices, boundedVertices;
    
    for (int i = 0; i < whiteVertices.size(); i++) {
        int numBlack = numberOfBlackNeighbors(whiteVertices.at(i), myS);
        if ( numBlack == 0) {
            superFreeVertices.push_back(whiteVertices.at(i));
        }
        else if (numBlack == 1) {
            freeVertices.push_back(whiteVertices.at(i));
        }
        else if (numBlack == 2) {
            boundedVertices.push_back(whiteVertices.at(i));            
        }
        else {
            cerr << "Graph is not claw-free.\n";
        }
    }
    
    // //3. Find maximum weighted white augmenting path between two distinct free vertices
    vector<int> maximumWWAP = maxWeightWhiteAugPath(superFreeVertices, freeVertices, boundedVertices, myS);
    //vector<int> maximumWWAP;

    //cout<<"Max weight white augmenting path"<<endl;
    //PrintVector(maximumWWAP);
    
    
    // //4. If no such path exists for S, then stop. S is optimal;
    if (maximumWWAP.size() == 0) {
        
        cout<<"#####################max weight co-2-plex########################:"<<endl;
        PrintVector(myS);
        //exit(0);
        
    }
    else{
       //5. S = S-P and return to (2)
        vector<int> myNewS(n); // 0  0  0 ... 0 (n of them)
        vector<int>::iterator itt;
        
        sort (myS.begin(),myS.end());
        sort (maximumWWAP.begin(),maximumWWAP.end());
        
        itt=std::set_symmetric_difference (myS.begin(),myS.end(), maximumWWAP.begin(), maximumWWAP.end(), myNewS.begin());
        //  myNewS ... 0  0  0  0
        myNewS.resize(itt-myNewS.begin());//only white vertices
        
        /*cout << "The symmetric difference of S and P has " << (myNewS.size()) << " elements:\n";
        for (itt=myNewS.begin(); itt!=myNewS.end(); ++itt)
            cout << ' ' << *itt;
        cout << '\n';*/
        
		Minty(myNewS);
    }
}

int main(int argc, char* argv[]){
    //cout << argv[1] << argv[2];
    if(argc != 3) {
        cerr<<" USAGE: ./executable graphAdjacency graphWeights"<<endl;
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
		// need to get vertices of g from G
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