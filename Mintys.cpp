/*
    Minty's Algorithm for the Max Weight Stable Set on a claw-free Graph
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

	struct CCdatagroup;  //might need file that creates this struct
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

Graph::Graph(char* inputFile, char* weightFile){
    ifstream f_in, f_win;
	f_in.open(inputFile);
	
	f_in>>n;//number of vertices
	f_in>>m;//number of edges
    //cout << "Number of Nodes: " << n << "\nNumber of Edges: " << m << "\n";
	
	numMatchingProblem = 0;

    int *x = new int[m];
	int *y = new int[m];

    for(int i = 0; i < m; i++){
		f_in>>x[i]>>y[i];  //Cynthias code stores a third element in a temp variable z
	}
	f_in.close();
    // for(int i = 0; i < m; i++){
	// 	cout << x[i] << " " << y[i] << "\n";
	// }

    int w;
    f_win.open(weightFile);
    for (int i = 0; i < n; i++) {
        f_win>>w; // file just has weights (use f_win>>z>>w if it has nodes and weights)
        //nodes.push_back(z);
        weights.push_back(w);
    }
    // for(int i = 0; i < n; i++){
	// 	cout << nodes[i] << " " << weights[i] << "\n";
	// }
    
    f_win.close();
	
	//Create the attribute C matrix
	C = new int*[n];
	for(int i = 0; i < n; i++) {
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
        //S.push_back(i); // if you want to see the nodes pushed (not supposed to be in S)
    }
    //PrintVector(S);
    Minty(S);
}

bool Graph::adjacent(vector<int> A, vector<int> B){
	// If a node in A is adjacent to a node in B, return true
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
        if (C[vertex][blackVertices.at(i)]) {
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
	
	sort(blackVertices.begin(), blackVertices.end());
	//printVector(blackVertices)

    //2.1 Generate all white alternating paths of length 0 or 2
    
    //Save SF vertex with maximum weight only
    int maxWeightl0 = 0; // Changed these to 0 (from -1) so that an empty vector is returned if all paths are nonpositive
	int maxWeightln = 0;
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

	//Generate WAPs of length 2 and save it only if its weight is greater than the current max weight
    vector<int> weightWhite;
    vector<int> weightBlack;
    for (int i = 0; i < blackVertices.size(); i++) {
        int adjacent1 = -1;
        int adjacent2 = -1;
        weightBlack.push_back(weights.at(blackVertices.at(i)));
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
				//xa and xb are the black vertices adjacent to a and b, respectively
				int xa = -1;
				int xb = -1;
                for (int i = 0; i < blackVertices.size(); i++) {
                    if (C[a][blackVertices.at(i)] && xa == -1) {
						xa = blackVertices.at(i);
                    }
                    
                    if (C[b][blackVertices.at(i) && xb == -1]) {
						xb = blackVertices.at(i);
                    }
					if (xa != -1 && xb != -1) {
						break;
					}
				}
                
				if (xa != xb) {
					//Find MWWAP
					vector< vector<int> > EdmondsG;
					int N = 0; // number of nodes in the Edmonds graph
					EdmondsG = Edmonds(F.at(a), F.at(b), xa, xb, blackVertices, B, N);

					sort(EdmondsG.begin(), EdmondsG.end()); //doesn't rearrange each component vertex, right?
                    					
					if (N > EdmondsG.size()) {
						double matzeit = 0.0;
						double genzeit = 0.0;
						int ncount, ecount;
						//long l;
						time_t l; //l must be of this type
						
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
                                //Negative since we are solving for the max weighted stable set and perfect_match 
								//solves for min weight perfect matching
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
						
                        
						/*Returns 0 if it worked and 1 otherwise (for example, when one
						of the mallocs failed). The nodes in the graph should be named 
						0 through #nodes - 1 */
						// Generalize to maximum matching instead of perfect matching
						// Perfect match: every vertex is adjacent to 1 edge
						if (perfect_match (ncount, NULL, ecount, &elist, &elen, &ourMatch, &ourWeights,
						   blossom_file, match_file, just_frac, no_frac, 
						   use_all_trees, partialprice, &matzeit)) {
							fprintf(stderr, "perfect_match failed\n");
			
						}

						//TODO: M delta M* contains a maximum weight alternating path between a^ and b^
						//Need to extract this path and store the vertices appropriately
						//M is the set of black edges in the Edmonds graph
						//M* is the maximum weight matching found
						
						//The mapping is a pair entry(i)-> entries 2i and 2i+1 weight
						//Entries 2i and 2i+1 correspond to node i in the original graph
						for (int i = 0; i < ncount; i++) {
							if (ourMatch[i] % 2 == 0) {
								maxWWAPln.push_back(ourMatch[i]/2);
							}
							else {
								maxWWAPln.push_back((ourMatch[i]-1)/2);
							}

						}
						
						int tempMaxln = 0;
						for (int i = 0; i < maxWWAPln.size(); i++) {
							tempMaxln += weights.at(maxWWAPln.at(i));
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

	//if (maxWWAPln > maxWWAPl0) {
	if (maxWeightln > maxWeightl0) {
		return maxWWAPln;
	}
	
    return maxWWAPl0;
}

vector<bool> Graph::iwapBF(vector<int> &weight, vector<int> &whiteReachIrreg, vector<int> whiteSource, vector<int> whiteSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices){

	// weight must be empty and will contain the weight of IWAPS 
	//if the first IWAP between the source and a whiteSink is whiteSink.at(2) (so third entry)
	//then the first entry of weight will be the entry (weight?) of such a path
	//i.e., weight only stores the weights of the IWAPs, not for each path
	vector<bool> IWAP( whiteSinks.size(), false);
	vector<int> vertexSet;
	//Construct a DAG (directed acyclic graph) to solve the problem
	vector< vector<int> > level;
	level.push_back(whiteSource);
	vertexSet.push_back(whiteSource.at(0));//We're doing contraction
	bool moreLevelsNeeded = true;
	bool white = false;//search through white vertices if true
	
	vector<int> edgeSet;
	vector<int> edgeWeights; //Bellman-Ford finds shortest path
	
	//only save two vertices edgeSet even if they are a pair, but do consider the adjacencies of the whole thing
	
	vector< vector<int> > irregularVP; // can maybe make this vector<int> since we don't have pairs anymore
	for (int i = 0; i < irregularVertices.size(); i++) {
		vector<int> t;
		
		t.push_back(irregularVertices.at(i));
		
		irregularVP.push_back(t);
		t.clear();
	}
	
	// for (int i = 0 ; i < irregularPairs.size(); i++) {
	// 	irregularVP.push_back(irregularPairs.at(i));
	// }
	
	vector<int> whiteIndices;
	while (moreLevelsNeeded) {
		if (white) { // not executed at fist since white is initialized to false
		
			int counter = 0;
			vector< vector<int> > tempLevel;
			for (int i = 0 ; i < level.size(); i++) {
				//level contains all the vertices in the current level, each as a vector
				bool inWhiteLevel = false;
				// Go through the eligible wings
				for (int j = 0; j < whiteVerticesInWings.size(); j++) {
					vector<int> tempW;
					tempW.push_back(whiteVerticesInWings.at(j));
					//because level[i] and tempW are vectors of one vertex, can do C[level[i].at(0)][tempW.at(0)]
					//should whiteInWings store all the vertices in each wing then?
					//at least one vertex in level is adjacent to the one vertex in tempW
					if (adjacent(level[i],tempW)) {
						counter++;
						inWhiteLevel = true;
						edgeSet.push_back(level[i].at(0));
						edgeSet.push_back(tempW.at(0)); //original code only pushes one vertex as well
						vertexSet.push_back(tempW.at(0));
						// edgeWeights.push_back(-weights.at(tempW.at(0)));
						//The weight should be negative because it's a white vertex and we're finding shortest path
						edgeWeights.push_back(-weights.at(tempW.at(0)));
						
					}
					
					
					if (inWhiteLevel) {
					
						for (int k = 0; k < whiteSinks.size(); k++) {
							if (whiteSinks.at(k) == whiteVerticesInWings.at(j)) {
								IWAP.at(k) = true;
								whiteIndices.push_back(vertexSet.size()-1); //because 0-based indexing
								break;
							}
						}
						
						//erase entry j from possible white vertices
						//can do so because a white vertex cannot be adjacent to two irregular vertices in the
						//same level and won't appear in the next level
						whiteVerticesInWings.erase(whiteVerticesInWings.begin() + j);
						//When you remove you decrease the size of the vector
						j--;
						inWhiteLevel = false;
						
						tempLevel.push_back(tempW);
						
					}
					tempW.clear();
				}
				
			}
			if (counter == 0) {
				moreLevelsNeeded = false;
			}
			
			level.clear();
			level = tempLevel;
			tempLevel.clear();
			white = false;
		}
		else {
			int counter2 = 0;
			vector< vector<int> > tLevel;
			
			for (int i = 0; i < level.size(); i++) {
				bool inBlackLevel = false;
				for (int j = 0 ; j < irregularVP.size(); j++) {
					//because level[i] and irregularVP[j] are vectors of one vertex, can do C[level[i].at(0)][irregularVP[j].at(0)]
					if (adjacent(level[i], irregularVP[j])) { // adjacent returns true if any node in the first is adjacent to any node in the second
						counter2++;
						inBlackLevel = true;
						edgeSet.push_back(level[i].at(0));
						edgeSet.push_back(irregularVP[j].at(0)); //original code only pushes one vertex as well
						vertexSet.push_back(irregularVP[j].at(0));
						// int w2 = 0;
						// for (int jj = 0; jj < irregularVP[j].size(); jj++) {
						// 	w2 -= weights.at(irregularVP[j].at(jj));
						// }
						// edgeWeights.push_back(w2);

						//There should only be only element in each irregularVP[j]
						//Weight should be positive because we're finding shortest path
						edgeWeights.push_back(weights.at(irregularVP[j].at(0)));
					}
					
					if (inBlackLevel) {
						tLevel.push_back(irregularVP.at(j));//this pushes a vector (of one vertex) into tLevel
					
						irregularVP.erase(irregularVP.begin() + j);
						j--;
						
						inBlackLevel = false;
						
						
					}
				}
				
			}
			
			if (counter2 == 0) {
				moreLevelsNeeded = false;
			}
			
			level.clear();
			level = tLevel;
			tLevel.clear();
			white = true;
			
		}

		
	}
	
	
	//Now do bellman ford
	
	// Step 1: initialize graph
   /*for each vertex v in vertices:
       if v is source then distance[v] := 0
       else distance[v] := inf
       predecessor[v] := null*/
	   
	vector<int> distance;
	vector<int> predecessor;
	
	for (int i = 0; i < vertexSet.size(); i++) {
		if (vertexSet.at(i) == whiteSource.at(0)) {
			distance.push_back(0);
			//null = -1 (for practical purposes)
			predecessor.push_back(-1);
		}
		else {
			//infinity = 1e7 (for practical purposes)
			distance.push_back(10000000);
			//null = -1 (for practical purposes)
			predecessor.push_back(-1);
		}

	}
	
	
	// Step 2: relax edges repeatedly
   /*for i from 1 to size(vertices)-1:
       for each edge (u, v) with weight w in edges:
           if distance[u] + w < distance[v]:
               distance[v] := distance[u] + w
               predecessor[v] := u*/
	for (int i = 0; i < vertexSet.size()-1; i++) {
		for (int j = 0; j < edgeSet.size()-2; j+=2) {
			int entry1, entry2;	
			int countEntries = 0;
			for (int k = 0; k < vertexSet.size(); k++) {
				if (edgeSet.at(j) == vertexSet.at(k)) {
					entry1 = k;
					countEntries++;
				}
				if (edgeSet.at(j+1) == vertexSet.at(k)) {
					entry2 = k;
					countEntries++;
				}
				if (countEntries ==2) {
					break;
				}
				
			}
				
				
			if (distance.at(entry1) + edgeWeights.at(j/2) < distance.at(entry2)) {
				
				distance.at(entry2) = distance.at(entry1) + edgeWeights.at(j/2);
				
				//predecessor.at(entry2) = edgeSet.at(entry1); //edgeSet.at(j) ? or vertexSet.at(entry1)?
				predecessor.at(entry2) = edgeSet.at(j);
			}
		}
	}
	
	// Step 3: check for negative-weight cycles
	//NO need to check this since we have a DAG
			   

	//Used saved index of white vertices to figure out the weight of the IWAPs
	
	for (int i = 0; i < whiteIndices.size(); i++) {
	
		//weight.push_back(distance.at(whiteIndices.at(i)));
		weight.push_back(-distance.at(whiteIndices.at(i)) + weights.at(whiteSource.at(0)));
	}
	
	return IWAP;
}

vector<int> Graph::edmondsCorrection1(vector<int> regNeib, vector<int> whiteWings, vector<int> irreg){
	vector<int> whiteReachableWingsThroughIrreg;
	bool moreLevels = false;
	
	vector<int> blackLevel;
	bool blackLevelFound = false;
	for (int i = 0; i < irreg.size(); i++) {
		for (int j = 0; j < regNeib.size(); j++) {
			if (C[regNeib.at(j)][irreg.at(i)]) {
				blackLevelFound = true;
				whiteReachableWingsThroughIrreg.push_back(regNeib.at(j));
								
			}
		}

		if (blackLevelFound) {
			moreLevels = true;
			blackLevel.push_back(irreg.at(i));
			//erase entry i from possible black vertices
			irreg.erase(irreg.begin() + i);
			//When you remove you decrease the size of the vector
			i--; //doesn't actually matter because we break
			//A black level can only have one vertex (because irregular)
			break;
		}
		
	}
	
	// bool moreLevels = true;

	/********** Need to remove white vertices in wings already visited *********/
	/* Can pass whichWings and remove vertices by wings, or a set difference*/
	//Case: whiteWings is vector<int>
	//vector<int> whiteWingsDiff(whiteWings.size() + regNeib.size()); //0 0 ... 0
	vector<int>::iterator it;
	sort(whiteWings.begin(), whiteWings.end());
	sort(regNeib.begin(), regNeib.end());
	it = std::set_difference(whiteWings.begin(), whiteWings.end(), regNeib.begin(), regNeib.end(), whiteWings.begin());
	whiteWings.resize(it-whiteWings.begin());

	//Case: whiteWings is vector< vector<int> >
	//whiteWings = set(whiteWings) - set(regNeib);


	vector<int> whiteLevel;
	while (moreLevels) {
		
		for (int i = 0; i < whiteWings.size(); i++) {
			//for (int j = 0; j < whiteWings[i].size(); j++) {
			if (C[blackLevel.at(0)][whiteWings.at(i)]) {
				whiteLevel.push_back(whiteWings.at(i));
				
				//erase entry i from possible white vertices
				whiteWings.erase(whiteWings.begin() + i);
				//When you remove you decrease the size of the vector
				i--;
			}

				// if (C[blackLevel.at(0)][whiteWings[i].at(j)]) {
				// 	whiteLevel.push_back(whiteWings[i].at(j));
					
				// 	//erase entry j from possible white vertices
				// 	whiteWings.erase(whiteWings[i].begin() + j);
				// 	//When you remove you decrease the size of the vector
				// 	j--;
				// }

			//}
		}
		
		blackLevel.clear();
		
		for (int i = 0; i < whiteLevel.size(); i++) {
			for (int j = 0; j < irreg.size(); j++) {
				if (C[whiteLevel.at(i)][irreg.at(j)]) {
					whiteReachableWingsThroughIrreg. push_back(whiteLevel.at(i));
					
					blackLevel. push_back(irreg.at(j));
					//erase entry j from possible black vertices
					irreg.erase(irreg.begin() + j);
					//When you remove you decrease the size of the vector
					j--;
					
				}
			}
			
		}

		/********** Need to remove white vertices in wings already visited *********/
		//Case: whiteWings is vector<int>
		//vector<int> whiteWingsDiff(whiteWings.size() + regNeib.size()); //0 0 ... 0
		vector<int>::iterator it;
		sort(whiteWings.begin(), whiteWings.end());
		sort(whiteLevel.begin(), whiteLevel.end());
		it = std::set_difference(whiteWings.begin(), whiteWings.end(), whiteLevel.begin(), whiteLevel.end(), whiteWings.begin());
		whiteWings.resize(it-whiteWings.begin());

		//Case: whiteWings is vector< vector<int> >
		// CODE
		
		whiteLevel.clear();

		if (blackLevel.size() == 0) {
			moreLevels = false;
		}
	}

	
	// vector< vector<int> > blackLevel;
	// for (int i = 0; i < regNeib.size(); i++) {
	// 	for (int j = 0; j < irreg.size(); j++) {
	// 		if (C[regNeib.at(i)][irreg.at(j)]) {
	// 			whiteReachableWingsThroughIrreg.push_back(regNeib.at(i));
	// 			vector<int> tempBL;
				
	// 			tempBL. push_back(irreg.at(j));
	// 			blackLevel. push_back(tempBL);
	// 			tempBL.clear();
	// 			//erase entry j from possible black vertices
	// 			irreg.erase(irreg.begin() + j);
	// 			//When you remove you decrease the size of the vector
	// 			j--;
				
	// 		}
	// 	}
		
	// 	// vector<int> tempNbr;
	// 	// tempNbr.push_back(regNeib.at(i));
	// 	// for (int j = 0; j < irregPairs.size(); j++) {
			 
	// 	// 	if (adjacent(tempNbr, irregPairs.at(j))) {
			
	// 	// 		whiteReachableWingsThroughIrreg. push_back(regNeib.at(i));
	// 	// 		vector<int> tempBL2;
				
	// 	// 		tempBL2. push_back(irreg.at(j));
	// 	// 		blackLevel. push_back(tempBL2);
	// 	// 		tempBL2.clear();
				
	// 	// 		//erase entry j from possible black pairs
	// 	// 		irregPairs.erase(irregPairs.begin() + j);
	// 	// 		//When you remove you decrease the size of the vector
	// 	// 		j--;
				
	// 	// 	}
	// 	// }
	// 	// tempNbr.clear();
	// }
	
	// bool moreLevels = true;

	// vector<int> whiteLevel;
	// while (moreLevels) {
		
	// 	for (int i = 0; i < blackLevel.size() ; i++) {
	// 		for (int j = 0; j < whiteWings.size(); j++) {
	// 			vector<int> tWhite;
	// 			tWhite.push_back(whiteWings.at(j));
				
	// 			if (adjacent(blackLevel.at(i), tWhite)) {
	// 				whiteLevel.push_back(whiteWings.at(j));
					
	// 				//erase entry j from possible white vertices
	// 				whiteWings.erase(whiteWings.begin() + j);
	// 				//When you remove you decrease the size of the vector
	// 				j--;
	// 			}
	// 			tWhite.clear();
	// 		}	
	// 	}
		
	// 	blackLevel.clear();
		
	// 	for (int i = 0; i < whiteLevel.size(); i++) {
	// 		for (int j = 0; j < irreg.size(); j++) {
	// 			if (C[whiteLevel.at(i)][irreg.at(j)]) {
	// 				whiteReachableWingsThroughIrreg. push_back(whiteLevel.at(i));
	// 				vector<int> tempBL2;
					
	// 				tempBL2. push_back(irreg.at(j));
	// 				blackLevel. push_back(tempBL2);
	// 				tempBL2.clear();
	// 				//erase entry j from possible black vertices
	// 				irreg.erase(irreg.begin() + j);
	// 				//When you remove you decrease the size of the vector
	// 				j--;
					
	// 			}
	// 		}
			
			
	// 		// vector<int> tempWht;
	// 		// tempWht. push_back(whiteLevel.at(i));
	// 		// for (int j = 0; j < irregPairs.size(); j++) {
				 
	// 		// 	if (adjacent(tempWht, irregPairs.at(j))) {
				
	// 		// 		whiteReachableWingsThroughIrreg. push_back(whiteLevel.at(i));
	// 		// 		vector<int> tempBL2;
					
				
	// 		// 		blackLevel. push_back(irregPairs.at(j));
					
					
	// 		// 		//erase entry j from possible black pairs
	// 		// 		irregPairs.erase(irregPairs.begin() + j);
	// 		// 		//When you remove you decrease the size of the vector
	// 		// 		j--;
					
	// 		// 	}
	// 		// }
	// 		// tempWht.clear();
	// 	}
	// 	whiteLevel. clear();
		
	// 	if (blackLevel.size() == 0) {
	// 		moreLevels = false;
	// 	}
	// }
	
	vector<int>::iterator itrW;
	itrW = unique( whiteReachableWingsThroughIrreg.begin(), whiteReachableWingsThroughIrreg.end());
	whiteReachableWingsThroughIrreg.resize( std::distance(whiteReachableWingsThroughIrreg.begin(),itrW));	

	return whiteReachableWingsThroughIrreg;
}

int Graph::edmondsCorrection3(int y11, int yl1, int y12, int yl2, vector<int> whiteWings, vector<int> irreg, vector< vector<int> > irregPairs, vector<int> &k, int &wP11i, int &wP22i, int &wP11j, int &wP22j){
	int zk = 0;
	
	bool moreLevels = true;
	bool zkFound = false;
	

	vector<int> whiteLevelL1;
	vector<int> whiteLevelL2;
	
	vector< vector<int> > blackLevelL1;
	vector< vector<int> > blackLevelL2;
	
	whiteLevelL1. push_back(y11);
	//wP11i += weights. at(y11);
	whiteLevelL2. push_back(y12);
	//wP22i += weights. at(y12);
	while (moreLevels) {
	
		if (!zkFound && (blackLevelL1. size() > 0 || blackLevelL2. size() > 0)) {
		
			vector<int> bLevelL1;
			vector<int> bLevelL2;
			for (int i = 0; i < blackLevelL1. size(); i++) {
				for (int j = 0; j < blackLevelL1[i]. size(); j++) {
					bLevelL1. push_back(blackLevelL1[i]. at(j));
				}
			
			}
			
			for (int i = 0; i < blackLevelL2. size(); i++) {
				for (int j = 0; j < blackLevelL2[i]. size(); j++) {
					bLevelL2. push_back(blackLevelL2[i]. at(j));
				}
			
			}
	
			vector<int> symmBlack(n);	// 0  0  0 ... 0 (n of them)
			vector<int>::iterator it;
			sort(bLevelL1. begin(), bLevelL1. end());
			sort(bLevelL2. begin(), bLevelL2. end());
			it=std::set_symmetric_difference (bLevelL1 .begin(),bLevelL1.end(), bLevelL2 .begin(),bLevelL2 .end(), symmBlack .begin());
			symmBlack. resize(it - symmBlack. begin());
			
			if (symmBlack. size() == 1) {
				zkFound = true;
				k. push_back( symmBlack. at(0));
			}
			else if ( symmBlack. size() == 2 && C[symmBlack. at(0)][symmBlack. at(1)] ) {
				zkFound = true;
				k. push_back( symmBlack. at(0));
				k. push_back( symmBlack. at(1));
			
			}
			else if ( symmBlack. size() > 1 ) {
				zkFound = true;
				k. push_back( symmBlack. at(0));
				
				for (int i = 1; i < symmBlack. size(); i++) {
					if (C[symmBlack. at(0)][ symmBlack. at(i)] ) {
						k. push_back( symmBlack. at(i));
						break;
					}
				}
				
			
			}
		
		}
	
		int irrSz = irreg.size() + irregPairs.size();
		vector<bool> bL1( irrSz, false);
		vector<bool> bL2( irrSz, false);
		
		vector<bool> wL1( whiteWings.size(), false);
		vector<bool> wL2( whiteWings.size(), false);
		
		int maxWhiteWeightL1 = 0;
		int maxWhiteWeightL2 = 0;
		for (int i = 0; i < whiteLevelL1.size(); i++) {
			vector<int> tempWhiteL1;
			tempWhiteL1. push_back(whiteLevelL1.at(i));
			for (int j = 0; j < irreg.size(); j++) {
				vector<int> tBlack;
				tBlack. push_back(irreg.at(j));
				if (adjacent(tempWhiteL1, tBlack)) {
					blackLevelL1. push_back(tBlack);
					
					//To erase entry j from possible white vertices
					bL1.at(j) = true;
					
					//For correction 3
					if (maxWhiteWeightL1 < weights. at(tempWhiteL1. at(0))) {
						maxWhiteWeightL1 = weights. at(tempWhiteL1. at(0));
					}
					
				}
				
				tBlack.clear();
				
			}
			
			for (int j = 0; j < irregPairs.size(); j++) {
				if (adjacent(tempWhiteL1, irregPairs.at(j))) {
					blackLevelL1. push_back(irregPairs.at(j));
					bL1.at(irreg.size()+j) = true;
				}
			}
			
			
			tempWhiteL1.clear();
		}
		
		//For correction 3
		if (maxWhiteWeightL1 > 0 && !zkFound) {
			wP11i += maxWhiteWeightL1;
		}
		else if(maxWhiteWeightL1 > 0){
			wP11j += maxWhiteWeightL1;
		}

		
		for (int i = 0; i < whiteLevelL2.size(); i++) {
			vector<int> tempWhiteL2;
			tempWhiteL2. push_back(whiteLevelL2.at(i));
			for (int j = 0; j < irreg.size(); j++) {
				vector<int> t2Black;
				t2Black. push_back(irreg.at(j));
				if (adjacent(tempWhiteL2, t2Black)) {
					blackLevelL2. push_back(t2Black);
					
					//To erase entry j from possible white vertices
					bL2.at(j) = true;
					
					
					//For correction 3
					if (maxWhiteWeightL2 < weights. at(tempWhiteL2. at(0))) {
						maxWhiteWeightL2 = weights. at(tempWhiteL2. at(0));
					}
					
				}
				t2Black.clear();
			}
			
			for (int j = 0; j < irregPairs.size(); j++) {
				if (adjacent(irregPairs.at(i), tempWhiteL2)) {
					blackLevelL2. push_back(irregPairs. at(j));
					bL2.at(irreg.size() +j) = true;
				}
			}
			tempWhiteL2. clear();
		}
		
		//For correction 3
		if (maxWhiteWeightL2 > 0 && !zkFound) {
			wP22i += maxWhiteWeightL1;
		}
		else if ( maxWhiteWeightL2 > 0 ) {
			wP22j += maxWhiteWeightL1;
		}
		
		for (int i = 0; i < irrSz; i++) {
			if (i < irreg.size() && (bL1.at(i) == true || bL2.at(i) == true)) {
				irreg. erase(irreg.begin() + i);
				i--;
			}
			
			if (i >= irreg.size() && (bL1.at(i) == true || bL2.at(i) == true)) {
				irregPairs. erase(irregPairs.begin() + i - irreg.size());
				i--;
			}
		}
		whiteLevelL1. clear();
		whiteLevelL2. clear();
		
		
		int minBlackWeightL1 = 0;
		int minBlackWeightL2 = 0;
		
		if (blackLevelL1. size() > 0) {
		
			for (int i = 0; i < blackLevelL1[0].size(); i++) {
				minBlackWeightL1 += weights. at(blackLevelL1[0]. at(i));
			}
		}
		
		if (blackLevelL2. size() > 0) {
		
			for (int i = 0; i < blackLevelL2[0].size(); i++) {
				minBlackWeightL2 += weights. at(blackLevelL2[0]. at(i));
			}
		}
		
		
		
		for (int i = 0; i < blackLevelL1.size(); i++) {
			for (int j = 0; j < whiteWings.size(); j++) {
				vector<int> tWhite;
				tWhite. push_back(whiteWings.at(j));
				
				if (adjacent(blackLevelL1.at(i), tWhite)) {
					whiteLevelL1. push_back(whiteWings.at(j));
					
					//To erase entry j from possible white vertices
					wL1.at(j) = true;
					
					//For correction 3
					int weightCorr3 = 0;
					
					for (int k = 0; k < blackLevelL1[i].size(); k++) {
						weightCorr3 += weights. at(blackLevelL1[i]. at(k));
					}
					
					if (minBlackWeightL1 > weightCorr3) {
						//For correction 3
						minBlackWeightL1 = weightCorr3;
						
					}
					
				}
				tWhite.clear();
			}
		}
		
		//For correction 3
		if (minBlackWeightL1 > 0 && !zkFound) {
			wP11i -= minBlackWeightL1;
		}
		else if ( minBlackWeightL1 > 0 ) {
			wP11j -= minBlackWeightL1;
		}

		
		for (int i = 0; i < blackLevelL2.size(); i++) {
			for (int j = 0; j < whiteWings.size(); j++) {
				vector<int> t2White;
				t2White. push_back(whiteWings.at(j));
				
				if (adjacent(blackLevelL2.at(i), t2White)) {
					whiteLevelL2. push_back(whiteWings.at(j));
					
					//To erase entry j from possible white vertices
					wL2.at(j) = true;
					
					//For correction 3
					int weightCorr3 = 0;
					
					for (int k = 0; k < blackLevelL2[i].size(); k++) {
						weightCorr3 += weights. at(blackLevelL2[i]. at(k));
					}
					
					if (minBlackWeightL2 > weightCorr3) {
						//For correction 3
						minBlackWeightL2 = weightCorr3;
						
					}
					
					
				}
				t2White.clear();
			}
		}
		
		//For correction 3
		if (minBlackWeightL2 > 0 && !zkFound) {
			wP22i -= minBlackWeightL2;
		}
		else if ( minBlackWeightL2 > 0 ) {
			wP22j -= minBlackWeightL2;
		}
		
		
		
		
		for (int i = 0; i < whiteWings.size(); i++) {
		
			if (wL1.at(i) == true || wL2.at(i) == true) {
				whiteWings. erase(whiteWings.begin() + i);
				i--;
			}
			
		}
		
		if (whiteLevelL1.size() == 0 || whiteLevelL2.size() == 0) {
			moreLevels = false;
		}else {
			//Check if adjacent to yl1 or yl2 and then clear
			vector<int> tyl1, tyl2;
			tyl1. push_back(yl1);
			tyl2. push_back(yl2);
			int i = 0;
			int j = 0;
			for (i = 0; i < blackLevelL1.size(); i++) {
				if (adjacent(tyl1, blackLevelL1.at(i))) {
					break;
				}
			}
			
			for (j = 0; j < blackLevelL2.size(); j++) {
				if (adjacent(tyl2, blackLevelL2.at(j))) {
					break;
				}
			}
			
			if (i == j && i > 0 && j > 0) {
				//Then we have found the path we want and that path has length i=j
				moreLevels = false;
				zk = i;
			}
			
		
			blackLevelL1.clear();
			blackLevelL2.clear();
		}


	}
	
	
	return zk;
}

vector < vector<int> > Graph::Edmonds(int freeVertexA, int freeVertexB, int xa, int xb, vector<int> blackVertices, vector<int> boundedVertices, int &N){
    
    vector< vector<int> > edmondsG;
	
	int weightxA = weights.at(xa);
	int weightxB = weights.at(xb);
	
    /*We want to construct the reduced basic structure
     1.1 Ignore all the SF vertices (not passed to function)
     1.2 Ignore free vertices except a and b (not passed to function)
     1.3 Ignore all white vertices adjacent to a or b, since they will never appear in an alternating path
     */
    
	vector <int> whiteverticesNonAdjacentAorB;

	whiteverticesNonAdjacentAorB.resize(boundedVertices.size(),'\0'); // does not contain A or B
	copy(boundedVertices.begin(),boundedVertices.end(),whiteverticesNonAdjacentAorB.begin());
    sort(whiteverticesNonAdjacentAorB.begin(), whiteverticesNonAdjacentAorB.end());
    
    //cout<<"#########################################################"<<endl;
    //PrintVector(whiteverticesNonAdjacentAorB);
   
   // remove white vertices adjacent to A or B
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
    
    /* Now we form the reduced basic structure that is also the weight function. However, we dont have to get
	a copy of the weights of the nodes since we have global access to them (only for black vertices) <- ??? */
	vector<int> whiteRBS; //Contains the white vertices in the RBS, including A and B
    
    whiteRBS.resize(whiteverticesNonAdjacentAorB.size(), '\0');
    copy(whiteverticesNonAdjacentAorB.begin(), whiteverticesNonAdjacentAorB.end(), whiteRBS.begin());
    whiteRBS.push_back(freeVertexA);
    whiteRBS.push_back(freeVertexB);
    
	vector<int> blackRBS; //Contains the black vertices in the RBS

	blackRBS.resize(blackVertices.size(),'\0');
	copy(blackVertices.begin(), blackVertices.end(), blackRBS.begin());
    
    
    /***********************CLASSIFICATION OF BLACK VERTICES*************************/
    
    // A nonempty set of all bounded vertices, which are adjacent to the same two black vertices x and y, is called a wing
    // First figure out which black vertices are adjacent to bounded vertices (or vertices in whiteverticesNonAdjacentAorB)
    
    vector< int > numWings(blackRBS.size()); //keeps a counter to the number of wings each black vertex is adjacent to. They are 0 originally.
    
    vector< vector<int> > enumerateWings;//Keeps track of tip of the wings as pairs (ie, the two black vertices). Its size also gives how many different wings we have
	vector<int> whiteInWings; //Will store only one vertex that is in each distinct wing
	// enumerateWings and whiteInWings are the same size, and the i-th component of whiteInWings
	// is one vertex that is in the wing corresponding to the i-th component of enumerateWings
	vector< vector<int> > whiteInWingsAll; //Store all white vertices in each wing
	
    for (int i = 0; i < blackRBS.size(); i++) {
        int numi = numWings.at(i);
        for (int j = i+1; j < blackRBS.size(); j++) {
            int numj = numWings.at(j);
			bool wingAdded = false;
			vector<int> tempWhites;
            for (int k = 0; k < whiteverticesNonAdjacentAorB.size(); k++) {
                if (C[whiteverticesNonAdjacentAorB.at(k)][blackRBS.at(j)] && C[whiteverticesNonAdjacentAorB.at(k)][blackRBS.at(i)]) {
					
					if (!wingAdded) {
						numi++;
						numj++;
						vector<int> tempWing;
						tempWing.push_back(blackRBS.at(i));
						tempWing.push_back(blackRBS.at(j));
						sort(tempWing.begin(), tempWing.end());
						enumerateWings.push_back(tempWing);
						tempWing.clear();
						wingAdded = true;
					}

					tempWhites.push_back(whiteverticesNonAdjacentAorB.at(k));

                    // whiteInWings.push_back(whiteverticesNonAdjacentAorB.at(k));
                    // vector<int> tempWing;
                    // tempWing.push_back(blackRBS.at(i));
                    // tempWing.push_back(blackRBS.at(j));
                    // sort(tempWing.begin(), tempWing.end());
                    // enumerateWings.push_back(tempWing);
                    // tempWing.clear();
					// // the break statement here stops the search for other white vertices in the wing
					// // if you want to add those, need to make whiteInWings type vector< vector<int> >
                    // break;//as long as the set = wing is nonempty then we know a particular vertex is part of a wing
    
                }
            }
			if (wingAdded) {
				whiteInWingsAll.push_back(tempWhites);
				tempWhites.clear();
			}

            numWings.at(j) = numj;
        }

        numWings.at(i) = numi;
    }
    
    
	// Does this sorting affect the index correspondence of enumerateWings and whiteInWings?
	// Later, for a white vertex adjacent to a black vertex, find() is used to find it in whiteInWings
	// sort(enumerateWings.begin(), enumerateWings.end());
	// sort(whiteInWings.begin(), whiteInWings.end());

	// Only sort white vertices within each wing
	for (int i = 0; i < whiteInWingsAll.size(); i++) {
		sort(whiteInWingsAll[i].begin(), whiteInWingsAll[i].end());
	}
	
    //cout<<"***************************************************************"<<endl;
    //PrintVector(numWings);
    
	//Regular I: Vertices xa and xb
		//Stored as xa and xb
    //Regular II: A black vertex adjacent to 3 or more wings
    vector<int> regularII;
    //Irregular: Adjacent to exactly 2 wings
    vector<int> irregular;
	//useless otherwise
	vector<int> useless;
	
    for (int i = 0; i < blackRBS.size(); i++) {
		if (blackRBS.at(i) == xa || blackRBS.at(i) == xb) {
			continue;
		}
        if (numWings.at(i) > 2) {
            regularII.push_back(blackRBS.at(i));
			//cout << "regularII"<< blackRBS.at(i)<<endl;
        }
        if (numWings.at(i) == 2) {
            irregular.push_back(blackRBS.at(i));
			//cout << "irregular"<< blackRBS.at(i)<<endl;
        }

        // A useless vertex along with all adjoining white vertices can be deleted from blackRBS
		// as they will not appear any WAP between a and b.
		// I think it's fine since blackVertices still contains all the black vertices and
		// whiteRBS still contains all white vertices in the RBS, including A and B.
		
		// Actually, when doing so, other black vertices may then be misclassified.
        // Eg, You remove black vertex j. Black vertex i was irregular but is now useless, and it'll remain
        // irregular if it was already classified.
        // Best to just not remove them for now, but label as useless.
        if (numWings.at(i) == 1) {
			useless.push_back(blackRBS.at(i));
			//cout << "useless"<< blackRBS.at(i)<<endl;
		}
		//If not adjacent to any wings, remove it from blackRBS
		if (numWings.at(i) == 0) {
			blackRBS.erase(blackRBS.begin() + i);
			numWings.erase(numWings.begin() + i);
			i--;
		}
            
		// if (numWings.at(i) < 2) {
		// 	useless.push_back(blackRBS.at(i));
		// 	//cout << "useless"<< blackRBS.at(i)<<endl;

			// //If you want to remove all useless vertices and their adjacent white vertices
            // for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
            //     if (C[blackRBS.at(i)][whiteverticesNonAdjacentAorB.at(j)]) {
				
			// 		vector<int>::iterator itw;
			// 		itw = find(whiteInWings.begin(), whiteInWings.end(), whiteverticesNonAdjacentAorB.at(j));
            
			// 		if (itw != whiteInWings.end()) { // ie, if white vertex "j" is in whiteInWings
			// 			whiteInWings.erase(itw);
			// 		}
					
			// 		whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + j);
            //         j--;
	
            //     }
            // }
            // blackRBS.erase(blackRBS.begin() + i);
            // numWings.erase(numWings.begin() + i);
            // i--;
        // }
    }
    
    
    //whiteRBS becomes whiteverticesNonAdjacentAorB U freevertexA U freeVertexB
	//To update whiteRBS
	// whiteRBS.clear();
	// whiteRBS.resize(whiteverticesNonAdjacentAorB.size(), '\0');
    // copy(whiteverticesNonAdjacentAorB.begin(), whiteverticesNonAdjacentAorB.end(), whiteRBS.begin());
    // whiteRBS.push_back(freeVertexA);
    // whiteRBS.push_back(freeVertexB);

    //N1(xa) = freevertexA and N1(xb) = freevertexB
    //For each regular I vertex, put the free vertex into one class and the rest of its neighbors into the other class
    vector<int> n2xa;
    vector<int> n2xb;
    for(int i = 0 ; i < whiteverticesNonAdjacentAorB.size(); i++) {

        if (C[xa][whiteverticesNonAdjacentAorB.at(i)]) {
			n2xa.push_back(whiteverticesNonAdjacentAorB.at(i));    
        }
		
        if (C[xb][whiteverticesNonAdjacentAorB.at(i)]) {
			n2xb.push_back(whiteverticesNonAdjacentAorB.at(i));  
        }
		
    }
    
	//Store the neighbors of each regularII vertex
	vector< vector<int> > neighborsV;
    vector<int> tempRegII;
    for (int i = 0;  i < regularII.size(); i++) {    
		for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
			if (C[regularII.at(i)][whiteverticesNonAdjacentAorB.at(j)]) {
				tempRegII.push_back(whiteverticesNonAdjacentAorB.at(j));
			}
		}
		
		neighborsV.push_back(tempRegII); //neigborsV.at(i) then contains a vector of neighbors of regularII node i, which may empty
		tempRegII.clear();
    }
	
    
    /* Fact 2.2: For any regularII vertex, N(v) is uniquely partitioned into N1(v) and N2(v) so that for any x and y in distinct wings
	adjacent to v, x is not adjacent to y if and only if one of x and y is N1(v) and the other is in N2(v) */
    //Now we want to find all the pairs x and y in N(v) that belong to distinct wings

    /* 1. We already enumerated the wings (put the tip - ie, black vertices - of the wings on enumeratedWings and the size gives the number of distinct wings)
       2. Save in which wing each white vertex (excluding a and b) is contained */
    
	vector< vector <int> > whichWings(whiteverticesNonAdjacentAorB.size(), vector<int> (0)); //the 0 means empty
	//Gives index of the wing in enumeratedWings for which each vertex in whiteverticesNonAdjacentAorB belongs to
	//whichWings could perhaps be vector<int> since each white vertex should only be in one wing; would require lots of edits later
	
    for (int i = 0;  i < enumerateWings.size(); i++) {
        int enum0 = enumerateWings[i].at(0);
        int enum1 = enumerateWings[i].at(1);
        
		for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
			if (C[whiteverticesNonAdjacentAorB.at(j)][enum0] && C[whiteverticesNonAdjacentAorB.at(j)][enum1]) {
				whichWings[j].push_back(i); //each white vertex (excluding A and B) will be in one wing
			}
			
		}
    }

	/*//White vertices should be in exactly one wing
	for (int i=0; i < whichWings.size(); i++) {
		if (whichWings[i].size() != 1) {
			cout << "whiterverticesNonAdjacentAorB vertex" << i << "is in" << whichWings[i].size() << "wings." << endl;
		}
	} */
	
	/*cout << "whiteverticesNonAdjacentAorB size:"<< whiteverticesNonAdjacentAorB.size()<<endl;
	cout << "SECOND TIME whichWingsSssssssssssssssssssssssssssssssss"<<endl;
	PrintVectorVector(whichWings);
	cout<<endl;
    
	cout << "# of wings: "<<endl;;
	for (int i = 0 ; i < whichWings.size(); i++) {
		cout << "size:"<< whichWings[i].size()<<endl;
	}*/
	
	
    //4. Determine if two neighbor vertices, x and y, are in different wings and assign them to a class accordingly
	vector< vector<int> > nv1(regularII.size(), vector<int> (0));
    vector< vector<int> > nv2(regularII.size(), vector<int> (0));

    for (int i = 0;  i < neighborsV.size(); i++) {

        for (int j = 0 ; j < neighborsV[i].size(); j++) {
            for (int k = j+1; k < neighborsV[i].size(); k++){
			
                int x = -1;
                int y = -1;
                for (int l = 0; l < whiteverticesNonAdjacentAorB.size(); l++) {
                    if (neighborsV[i].at(j) == whiteverticesNonAdjacentAorB.at(l)) {
                        x = l;  // there is a 1-1 correspondence between whichWings and whiteverticeNonAdjacentAorB
				                // x is an index
					}
					
                    if (neighborsV[i].at(k) == whiteverticesNonAdjacentAorB.at(l)) {
                        y = l; // y is an index
                    }
					
					if (x != -1 && y != -1) {
						break;
					}
					
				}
				
				
				if (x!=-1 && y != -1) {
					// //Recall whichWings[i] contains the index to the wing that whiteverticesNonAdjacentAorB.at(i) belongs to
					// //Don't need to sort first because they should only contain one element
					// vector<int>::iterator it;
					// vector<int> symmDiffWing(enumerateWings.size()); // 0  0  0 ... 0 (enumerateWings.size() of them)
					// it=std::set_symmetric_difference(whichWings[x].begin(), whichWings[x].end(), whichWings[y].begin(), whichWings[y].end(), symmDiffWing.begin());
					
					// symmDiffWing.resize(it-symmDiffWing.begin());
					// //symmDiffWing should contain the indices of the wings that x and y are in, excluding those they're both in
					
					// //cout << "symDiff2:"<<endl;
					// //PrintVector(symmDiffWing);
					
					// int ww = whichWings[x].size() + whichWings[y].size();

					
					//Recall whichWings[i] contains the index to the wing that whiteverticesNonAdjacentAorB.at(i) belongs to
					// if (symmDiffWing.size() == ww) { // means x and y are not in any mutual wings
					if (whichWings[x].at(0) != whichWings[y].at(0)) {
						//x and y are in different wings
						if (!C[whiteverticesNonAdjacentAorB.at(x)][whiteverticesNonAdjacentAorB.at(y)]) {
							//x and y are nonadjacent

							// Original code had == instead of !=
							bool inNv1x = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(x) ) != nv1[i].end();
							bool inNv2x = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(x) ) != nv2[i].end();
							
							bool inNv1y = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(y) ) != nv1[i].end();
							bool inNv2y = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(y) ) != nv2[i].end();

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
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x)); //w.l.o.g
								nv2[i].push_back(whiteverticesNonAdjacentAorB.at(y));
							}
                        }
						else{
							//x and y are adjacent, so they go in the same set of neighbors (eg, N1(v) or N2(v))
							// Original code had == instead of !=
							bool inNv1x2 = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(x) ) != nv1[i].end();
							bool inNv2x2 = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(x) ) != nv2[i].end();
							
							bool inNv1y2 = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(y) ) != nv1[i].end();
							bool inNv2y2 = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(y) ) != nv2[i].end();

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
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x)); //w.l.o.g
								nv1[i].push_back(whiteverticesNonAdjacentAorB.at(y));
							}
							
						}
						
					}
					// I don't think the following is true; it doesn't hold for the example in the Minty paper
					// else{
					// 	//5. if they are in the same wing they can be in the same group of neighbors
					
					//  // Original code had == instead of !=
					// 	bool inNv1x3 = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(x)) != nv1[i].end();
					// 	bool inNv2x3 = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(x)) != nv2[i].end();
						
					// 	bool inNv1y3 = find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(y)) != nv1[i].end();
					// 	bool inNv2y3 = find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(y)) != nv2[i].end();
						
					// 	//Case 1. if whiteverticesNonAdjacentAorB.at(x) is n1v or n2v 
					// 	// and whiteverticesNonAdjacentAorB.at(y) is not in the other
					// 	if (inNv1x3 && !inNv2y3) {
					// 		nv1[i].push_back(whiteverticesNonAdjacentAorB.at(y));
					// 	}
						
					// 	if (inNv2x3 && !inNv1y3) {
					// 		nv2[i].push_back(whiteverticesNonAdjacentAorB.at(y));
					// 	}
						
					// 	//Case 2. if whiteverticesNonAdjacentAorB.at(y) is n1v or n2v 
					// 	//and whiteverticesNonAdjacentAorB.at(x) is not in the other
					// 	if (inNv1y3 && !inNv2x3) {
					// 		nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x));
					// 	}
						
					// 	if (inNv2y3 && !inNv1x3) {
					// 		nv2[i].push_back(whiteverticesNonAdjacentAorB.at(x));
					// 	}

					// 	//Case 3. If neither whiteverticesNonAdjacentAorB.at(x) 
					// 	// or whiteverticesNonAdjacentAorB.at(y) are in n1v or n2v
					// 	if (!inNv1x3 && !inNv2x3 && !inNv1y3 && !inNv2y3 ) {
					// 		nv1[i].push_back(whiteverticesNonAdjacentAorB.at(x)); //w.l.o.g
					// 		nv1[i].push_back(whiteverticesNonAdjacentAorB.at(y));
					// 	}
                    
					// }
					
				}
				
            }
        }
    }

	/*//Check it every vertex v is in exactly one class
	for (int j=0; j < whiteverticesNonAdjacentAorB.size(); j++) {
		int numClassesIn= 0;
		for (int i=0; i < regularII.size(); i++) {
			bool inN1 = false;
			bool inN2 = false;
			if (find(nv1[i].begin(), nv1[i].end(), whiteverticesNonAdjacentAorB.at(j)) != nv1[i].end()) {
				numClassesIn++;
				inN1 = true;
			}
			if (find(nv2[i].begin(), nv2[i].end(), whiteverticesNonAdjacentAorB.at(j)) != nv2[i].end()) {
				numClassesIn++;
				inN2 = true;
			}
			if (inN1 && in N2) {
				cout << "White vertex" << j << "is in N1 and N2 classes for white vertex" << i << endl;
			}
		}

		if (numClassesIn != 1) {
			cout << "White vertex" << j << "is in" << numClassesIn << "classes." << endl;
		}
	}*/

	/* Can also check if nv1(i).at(i) union nv2.at(i) equals neighborsV.at(i)*/
    
	//cout << "rbs before correction:"<< regularII.size() <<endl;
    // A regularII vertex with an empty node class, together with all adjoining white vertices, may be deleted
	// I'm guessing they never appear in a path

	//Going to not delete them for now as deleting them may mess up with the classification of black vertices
	// for (int i = 0; i < neighborsV.size(); i++) {
        
	// 	if (nv1[i].size() == 0 || nv2[i].size() == 0) {
            
	// 		regularII.erase(regularII.begin() + i);

	// 		for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
				
	// 			bool inN = find(neighborsV[i].begin(), neighborsV[i].end(), whiteverticesNonAdjacentAorB.at(j)) != neighborsV[i].end();
	// 			if (inN) { //if white vertex j is a neighbor of regularII node i	
	// 				//Would need to modify for whiteInWingsAll
	// 				vector<int>::iterator itw3;
	// 				itw3 = find(whiteInWings.begin(), whiteInWings.end(), whiteverticesNonAdjacentAorB.at(j));
		
	// 				if (itw3 != whiteInWings.end()) {
	// 					whiteInWings.erase(itw3);
	// 				}
	// 				whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + j);
	// 				j--;
	// 			}
	// 		}
			
	// 		neighborsV.erase(neighborsV.begin() + i);
	// 		nv1.erase(nv1.begin() + i);
	// 		nv2.erase(nv2.begin() + i);
	// 		i--;
	// 	}
    // }

	//cout << "regularII size:"<< regularII.size()<<endl;
	// Update whiteRBS
	// whiteRBS.clear();
	// whiteRBS.resize(whiteverticesNonAdjacentAorB.size(), '\0');
    // copy(whiteverticesNonAdjacentAorB.begin(), whiteverticesNonAdjacentAorB.end(), whiteRBS.begin());
    // whiteRBS.push_back(freeVertexA);
    // whiteRBS.push_back(freeVertexB);

    
	/***************** Construction of the Edmonds' Graph *****************/
	
	// R is the number of black vertices in the Edmonds graph
	int R = regularII.size() + 2; // 2 for the regularI vertices
	N = 2*R + 2; // number of vertices in the Edmonds graph (+2 for a^ and b^)
	//cout << "R:  "<< R<<endl;
    //We form a graph with 2R+2 vertices and R black branches
	//Each inner vector corresponds to an edge and consists of vertex1InEdge, vertex2InEdge, weightOfEdge
	vector< vector<int> > newGraphB;
	
	//a^ and b^ will be nodes 0 and 1 respectively
	//x_a^1, x_a^2 and x_b^1, x_b^2 will be 2,3,4,5
	
	//These two are white edges (a, x_a^1), (b, x_b^1)
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
	
	//cout << "regularII size:" << regularII.size()<<endl;
	//These are black edges
	for (int i = 0; i < regularII.size(); i++) {
		vector<int> tempE;
		tempE.push_back(2*(i+3)); //x_i^1 regularII vertex
		tempE.push_back(2*(i+3)+1); //x_i^2 regularII vertex
		tempE.push_back(weights.at(regularII.at(i)));
		newGraphB.push_back(tempE);
		tempE.clear();
    }

	int numberOfEdgesEdmonds = R + 2; // regularII + 2 regularI + 2 white; will increase as white edges are added

	//Add first the edges associated with regular I vertices (Note: (a,x_a^1) and (b,x_b^1) were added above)
	//N1(xa) = freeVertexA and N1(xb) = freeVertexB
    //For a regular vertex of the first kind put the free vertex into one class and the rest of its neighbors into the other class
	//Here we only consider n2xa and n2xb, since there will not be IWAPS starting at x_a^1 or x_b^1
	//iwapBF(vector<int> &weight, vector<int> &whiteReachIrreg, vector<int> whiteSource,
		// vector<int> whiteSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices);
	
	//Add edges between xa^2 and xb^2
	bool IWAP = false;
	int maxW = 0;
	vector<int> wxa;
	vector<int> whiteReachI;
	bool maxSet = false;
	for (int i = 0; i < n2xa.size(); i++) {
		vector<int> xaSource;
		xaSource.push_back(n2xa.at(i));
		vector<bool> IWAPBF = iwapBF(wxa, whiteReachI, xaSource, n2xb, whiteInWings, irregular);
		// int counter = 0;
		int tempMaxW;
		if (wxa.size() > 0) {
			IWAP = true;
			tempMaxW = *max_element(wxa.begin(), wxa.end());
			if (!maxSet) {
				maxW = tempMaxW;
				maxSet = true;
			}
			else if (tempMaxW > maxW) {
				maxW = tempMaxW;
			}
			// for (int j = 0; j < IWAPBF.size(); j++) {
			// 	if (IWAPBF.at(j)) {
			// 		IWAP = true;
			// 		// Only weights for IWAPs are stored
			// 		if (wxa.at(counter) > maxW) {
			// 			maxW = wxa.at(counter);
			// 			counter++;
			// 		}
			// 	}
			// }
		}
		
		xaSource.clear();
		wxa.clear();
	}
	
	if (IWAP) {
		vector<int> tempE3;
		//Recall node x_a^2 is 3 and x_b^2 is 5 
		tempE3.push_back(3);
		tempE3.push_back(5);
		tempE3.push_back(maxW);
		newGraphB.push_back(tempE3);
		numberOfEdgesEdmonds++;
		tempE3.clear();
	}
		

	// Add edges between x_a^2 and regularII nodes x_j
	bool IWAP1 = false;
	bool IWAP2 = false;
	vector<int> wxa1, wxa2;
	// First we add between x_a^2 and x_j^1
    for (int j = 0; j < nv1.size(); j++) {
		int maxW1 = 0;
		bool maxSet = false;
		for (int i = 0; i < n2xa.size(); i++){
			vector<int> xaaSource;
			xaaSource.push_back(n2xa.at(i));
			vector<bool> IWAPBF2 = iwapBF(wxa1,whiteReachI, xaaSource, nv1.at(j), whiteInWings, irregular);
			//int counter1 = 0;
			int tempMaxW;
			if (wxa1.size() > 0) {
				IWAP1 = true;
				tempMaxW = *max_element(wxa1.begin(), wxa1.end());
				if (!maxSet) {
					maxW1 = tempMaxW;
					maxSet = true;
				}
				else if (tempMaxW > maxW1) {
					maxW1 = tempMaxW;
				}
			}

			// for (int k = 0; k < IWAPBF2.size(); k++) {
			// 	if (IWAPBF2.at(k)) {
			// 		IWAP1 = true;
					
			// 		if (wxa1.at(counter1) > maxW1) {
			// 			maxW1 = wxa1.at(counter1);
			// 			counter1++;
			// 		}
			// 	}
			// }
			xaaSource.clear();
			wxa1.clear();
		}

		if (IWAP1) {
			vector<int> tempE4;
			//Recall node x_a^2 is 3 and regular vertices xj with N1 are mapped to vertex 2j (+6)
			tempE4.push_back(3);
			tempE4.push_back(2*(j+3));
			tempE4.push_back(maxW1);
			newGraphB.push_back(tempE4);
			numberOfEdgesEdmonds++;
			tempE4.clear();
			IWAP1 = false;
		}
	}

	// Next we add edges between x_a^2 and x_j^2
	for (int j = 0; j < nv2.size(); j++) {
		int maxW2 = 0;
		bool maxSet = false;
		for (int i = 0; i < n2xa.size(); i++){
			vector<int> xaaSource;
			xaaSource.push_back(n2xa.at(i));
			vector<bool> IWAPBF3 = iwapBF(wxa2, whiteReachI, xaaSource, nv2.at(j), whiteInWings, irregular);
			//int counter2 = 0;
			int tempMaxW;
			if (wxa2.size() > 0) {
				IWAP2 = true;
				tempMaxW = *max_element(wxa2.begin(), wxa2.end());
				if (!maxSet) {
					maxW2 = tempMaxW;
					maxSet = true;
				}
				else if (tempMaxW > maxW2) {
					maxW2 = tempMaxW;
				}
			}

			// for (int k = 0; k < IWAPBF3.size(); k++) {
			// 	if (IWAPBF3.at(k)) {
			// 		IWAP2 = true;
					
			// 		if (wxa2.at(counter2) > maxW2) {
			// 			maxW2 = wxa2.at(counter2);
			// 			counter2++;
			// 		}
			// 	}
			// }
			xaaSource.clear();
			wxa2.clear();
		}

		if (IWAP2) {
			vector<int> tempE5;
			//Node x_a^2 is 3 and regular vertices xj with N2 are mapped to vertex 2j+1 (+6)
			tempE5.push_back(3);
			tempE5.push_back(2*(j+3)+1);
			tempE5.push_back(maxW2);
			newGraphB.push_back(tempE5);
			numberOfEdgesEdmonds++;
			tempE5.clear();
			IWAP2 = false;
		}		
	}

	
	// Add edges between x_b^2 and regularII nodes x_j
	bool IWAP1b = false;
	bool IWAP2b = false;
	vector<int> wxb1, wxb2;
	// First we add between x_b^2 and x_j^1
	for (int j = 0; j < nv1.size(); j++) {
		int maxW1b = 0;
		bool maxSet = false;
    	for (int i = 0; i < n2xb.size(); i++) {
			vector<int> xbSource;
			xbSource.push_back(n2xb.at(i));
			vector<bool> IWAPBF4 = iwapBF(wxb1,whiteReachI, xbSource, nv1.at(j), whiteInWings, irregular);
			//int counter1 = 0;
			int tempMaxW;
			if (wxb1.size() > 0) {
				IWAP1b = true;
				tempMaxW = *max_element(wxb1.begin(), wxb1.end());
				if (!maxSet) {
					maxW1b = tempMaxW;
					maxSet = true;
				}
				else if (tempMaxW > maxW1b) {
					maxW1b = tempMaxW;
				}
			}

			// for (int k = 0; k< IWAPBF4.size(); k++) {
			// 	if (IWAPBF4.at(k)) {
			// 		IWAP1b = true;
					
			// 		if (wxb1.at(counter1) > maxW1b) {
			// 			maxW1b = wxb1.at(counter1);
			// 			counter1++;
			// 		}
			// 	}
			// }
			xbSource.clear();
			wxb1.clear();
		}

		if (IWAP1b) {
			vector<int> tempE6;
			//Recall node x_b^2 is 5 and regular vertices xj with N1 are mapped to vertex 2j (+6)
			tempE6.push_back(5);
			tempE6.push_back(2*(j+3));
			tempE6.push_back(maxW1b);
			newGraphB.push_back(tempE6);
			numberOfEdgesEdmonds++;
			tempE6.clear();
			IWAP1b = false;
		}
	}

	// Next we add between x_b^2 and x_j^2
	for (int j = 0; j < nv2.size(); j++) {
		int maxW2b = 0;
		bool maxSet = false;
		for (int i = 0; i < n2xb.size(); i++){
			vector<int> xbSource;
			xbSource.push_back(n2xb.at(i));
			vector<bool> IWAPBF5 = iwapBF(wxb2, whiteReachI, xbSource, nv2.at(j), whiteInWings, irregular);
			//int counter2 = 0;
			int tempMaxW;
			if (wxb2.size() > 0) {
				IWAP2b = true;
				tempMaxW = *max_element(wxb2.begin(), wxb2.end());
				if (!maxSet) {
					maxW2b = tempMaxW;
					maxSet = true;
				}
				else if (tempMaxW > maxW2b) {
					maxW2b = tempMaxW;
				}
			}

			// for (int k = 0; k < IWAPBF5.size(); k++) {
			// 	if (IWAPBF5.at(k)) {
			// 		IWAP2b = true;
					
			// 		if (wxb2.at(counter2) > maxW2b) {
			// 			maxW2b = wxb2.at(counter2);
			// 			counter2++;
			// 		}
			// 	}
			// }
			
			xbSource.clear();
			wxb2.clear();
		}
		
		if (IWAP2b) {
			vector<int> tempE7;
			//Recall node x_b^2 is 5 and regular vertices xj with N2 are mapped to vertex 2j+1 (+6)
			tempE7.push_back(5);
			tempE7.push_back(2*(j+3)+1);
			tempE7.push_back(maxW2b);
			newGraphB.push_back(tempE7);
			numberOfEdgesEdmonds++;
			tempE7.clear();
			IWAP2b = false;
		}
		
	}	
	
	//cout << "size of new graph B: "<< newGraphB.size()<<endl;
	
	// Add edges between regularII vertices x_i and x_j ---> four cases
	// Edges between x_i^1 and x_j^2
	vector<int> wnv12;
	bool IWAPnv12 = false;
	vector<int> whiteReachIW12;
	for (int i = 0; i < nv1.size(); i++) {
		for (int k = 0; k < nv2.size(); k++) {
			int maxNV12 = 0;
			bool maxSet = false;
			for (int j = 0 ; j < nv1[i].size(); j++) {
				vector<int> nv1Source;
				nv1Source.push_back(nv1[i].at(j));
				if (i != k) {
					vector<bool> IWAPBF8 = iwapBF(wnv12, whiteReachIW12, nv1Source, nv2.at(k), whiteInWings, irregular);
					//int counternv12 = 0;
					int tempMaxW;
					if (wnv12.size() > 0) {
						IWAPnv12 = true;
						tempMaxW = *max_element(wnv12.begin(), wnv12.end());
						if (!maxSet) {
							maxNV12 = tempMaxW;
							maxSet = true;
						}
						else if (tempMaxW > maxNV12) {
							maxNV12 = tempMaxW;
						}
					}

					// for (int l = 0; l < IWAPBF8.size(); l++) {
					// 	if (IWAPBF8.at(l)) {
					// 		IWAPnv12 = true;
							
					// 		if (wnv12.at(counternv12) > maxNV12) {
					// 			maxNV12 = wnv12.at(counternv12);
					// 			counternv12++;
					// 		}
					// 	}
					// }
					wnv12.clear();
					whiteReachIW12.clear(); //this doesn't get altered in iwapBF() at this moment
				}
				nv1Source.clear();
			}

			if (IWAPnv12) {
				vector<int> tempE8;
				//Recall vertices xi with N1 are mapped to vertex 2i (+6)
				//and xk with N2 are mapped to vertex 2k+1 (+6)
				tempE8.push_back(2*(i+3));
				tempE8.push_back(2*(k+3)+1);
				sort(tempE8.begin(), tempE8.end());
				tempE8.push_back(maxNV12);
				newGraphB.push_back(tempE8);
				numberOfEdgesEdmonds++;
				tempE8.clear();
				IWAPnv12 = false;
			}

		}
		
	}
	
	// Edges between x_i^2 and x_j^1
	vector<int> wnv21;
	bool IWAPnv21 = false;
	vector<int> whiteReachIW21;
	for (int i = 0; i < nv2.size(); i++) {
		for (int k = 0; k < nv1.size(); k++) {
			int maxNV21 = 0;
			bool maxSet = false;
			for (int j = 0 ; j < nv2[i].size(); j++) {
				vector<int> nv2Source;
				nv2Source.push_back(nv2[i].at(j));
				if (i != k) {
					vector<bool> IWAPBF9 = iwapBF(wnv21, whiteReachIW21, nv2Source, nv1.at(k), whiteInWings, irregular);
					//int counternv21 = 0;
					int tempMaxW;
					if (wnv21.size() > 0) {
						IWAPnv21 = true;
						tempMaxW = *max_element(wnv21.begin(), wnv21.end());
						if (!maxSet) {
							maxNV21 = tempMaxW;
							maxSet = true;
						}
						else if (tempMaxW > maxNV21) {
							maxNV21 = tempMaxW;
						}
					}

					// for (int l = 0; l < IWAPBF9.size(); l++) {
					// 	if (IWAPBF9.at(l)) {
					// 		IWAPnv21 = true;
							
					// 		if (wnv21.at(counternv21) > maxNV21) {
					// 			maxNV21 = wnv21.at(counternv21);
					// 			counternv21++;
					// 		}
					// 	}
					// }
					wnv21.clear();
					whiteReachIW21.clear(); //this doesn't get altered in iwapBF() at this moment
				}
				nv2Source.clear();
			}

			if (IWAPnv21) {
				vector<int> tempE9;
				//Recall vertices xi with N2 are mapped to vertex 2i+1 (+6)
				//and xk with N1 are mapped to vertex 2k (+6)
				tempE9.push_back(2*(i+3)+1);
				tempE9.push_back(2*(k+3));
				sort(tempE9.begin(), tempE9.end());
				tempE9.push_back(maxNV21);
				newGraphB.push_back(tempE9);
				numberOfEdgesEdmonds++;
				tempE9.clear();
				IWAPnv21 = false;
			}

		}
		
	}

	// Edges between x_i^1 and x_j^1
	vector<int> wnv11;
	bool IWAPnv11 = false;
	vector<int> whiteReachIW11;
	for (int i = 0; i < nv1.size(); i++) {
		for (int k = 0; k < nv1.size(); k++) {
			int maxNV11 = 0;
			bool maxSet = false;
			for (int j = 0 ; j < nv1[i].size(); j++) {
				vector<int> nv1Source;
				nv1Source.push_back(nv1[i].at(j));
				if (i != k) {
					vector<bool> IWAPBF10 = iwapBF(wnv11, whiteReachIW11, nv1Source, nv1.at(k), whiteInWings, irregular);
					//int counternv11 = 0;
					int tempMaxW;
					if (wnv11.size() > 0) {
						IWAPnv11 = true;
						tempMaxW = *max_element(wnv11.begin(), wnv11.end());
						if (!maxSet) {
							maxNV11 = tempMaxW;
							maxSet = true;
						}
						else if (tempMaxW > maxNV11) {
							maxNV11 = tempMaxW;
						}
					}

					// for (int l = 0; l < IWAPBF10.size(); l++) {
					// 	if (IWAPBF10.at(l)) {
					// 		IWAPnv11 = true;
							
					// 		if (wnv11.at(counternv11) > maxNV11) {
					// 			maxNV11 = wnv11.at(counternv11);
					// 			counternv11++;
					// 		}
					// 	}
					// }
					wnv11.clear();
					whiteReachIW11.clear(); //this doesn't get altered in iwapBF() at this moment
				}
				nv1Source.clear();
			}
			if (IWAPnv11) {
				vector<int> tempE10;
				//Recall vertices xi with N1 are mapped to vertex 2i (+6)
				//and xk with N1 are mapped to vertex 2k (+6)
				tempE10.push_back(2*(i+3));
				tempE10.push_back(2*(k+3));
				sort(tempE10.begin(), tempE10.end());
				tempE10.push_back(maxNV11);
				newGraphB.push_back(tempE10);
				numberOfEdgesEdmonds++;
				tempE10.clear();
				IWAPnv11 = false;
			}

		}
		
	}

	// Edges between x_i^2 and x_j^2
	vector<int> wnv22;
	bool IWAPnv22 = false;
	vector<int> whiteReachIW22;
	for (int i = 0; i < nv2.size(); i++) {
		for (int k = 0; k < nv2.size(); k++) {
			int maxNV22 = 0;
			bool maxSet = false;
			for (int j = 0 ; j < nv2[i].size(); j++) {
				vector<int> nv2Source;
				nv2Source.push_back(nv2[i].at(j));
				if (i != k) {
					vector<bool> IWAPBF11 = iwapBF(wnv22, whiteReachIW22, nv2Source, nv2.at(k), whiteInWings, irregular);
					//int counternv22 = 0;
					int tempMaxW;
					if (wnv22.size() > 0) {
						IWAPnv22 = true;
						tempMaxW = *max_element(wnv22.begin(), wnv22.end());
						if (!maxSet) {
							maxNV22 = tempMaxW;
							maxSet = true;
						}
						else if (tempMaxW > maxNV22) {
							maxNV22 = tempMaxW;
						}
					}

					// for (int l = 0; l < IWAPBF11.size(); l++) {
					// 	if (IWAPBF11.at(l)) {
					// 		IWAPnv22 = true;
							
					// 		if (wnv22.at(counternv22) > maxNV22) {
					// 			maxNV22 = wnv22.at(counternv22);
					// 			counternv22++;
					// 		}
					// 	}
					// }
					wnv22.clear();
					whiteReachIW22.clear(); //this doesn't get altered in iwapBF() at this moment
				}
				nv2Source.clear();
			}
			if (IWAPnv22) {
				vector<int> tempE11;
				//Recall vertices xi with N2 are mapped to vertex 2i+1 (+6)
				//and xk with N2 are mapped to vertex 2k+1 (+6)
				tempE11.push_back(2*(i+3)+1);
				tempE11.push_back(2*(k+3)+1);
				sort(tempE11.begin(), tempE11.end());
				tempE11.push_back(maxNV22);
				newGraphB.push_back(tempE11);
				numberOfEdgesEdmonds++;
				tempE11.clear();
				IWAPnv22 = false;
			}

		}
		
	}

	sort(newGraphB.begin(), newGraphB.end());


	//Correct Edmonds
	vector<int> whiteInIrregularWings;
	for (int i = 0; i < whiteverticesNonAdjacentAorB.size(); i++) {
		int countwiW =0;
		for (int j = 0; j < irregular.size(); j++) {
			if (C[whiteverticesNonAdjacentAorB.at(i)][irregular.at(j)]) {
				countwiW++;	
			}
		}
		
		// vector<int> tiiW;
		// tiiW. push_back(whiteverticesNonAdjacentAorB.at(i));
		// for (int j = 0; j < irregularPairs.size(); j++) {
		// 	if (adjacent(tiiW, irregularPairs.at(j))) {
		// 		countwiW++;
		// 	}
		// }
		// tiiW.clear();
		
		//Has to be adjacent to two irregular vertices to be in an irregular wing (can't be adjacent two more than 2 either)
		if (countwiW >=2) {
			whiteInIrregularWings.push_back(whiteverticesNonAdjacentAorB.at(i));
		}
	}
	
	//There shouldn't any duplicates, but can keep for now
	vector<int>::iterator itWiiW;
	itWiiW = unique(whiteInIrregularWings.begin(), whiteInIrregularWings.end());
	whiteInIrregularWings.resize(distance(whiteInIrregularWings.begin(),itWiiW) );
	
	sort(whiteInIrregularWings.begin(), whiteInIrregularWings.end());
	//Construct W(xi,xj):= Union of all wings reachable by irregular vertices to both xi and xj
	//So far we have white in wing, we want to save white in irregular wings

	
	int edmondsEdgesCorrected = numberOfEdgesEdmonds;
	//Make correction of Edmonds' graph according to Nakamura & Tamura
	//Lemma 3.2 correction
	//int regCorr = regularII.size() + regularIIPairs.size();
	int regCorr = regularII.size();
	for (int i = 0; i < regCorr; i++) {
		vector<int> wingsReachThroughIrregxi = edmondsCorrection1(neighborsV.at(i), whiteInWings, irregular);
		for (int j = i+1; j < regCorr-1; j++) {
			vector<int> wingsReachThroughIrregxj = edmondsCorrection1(neighborsV.at(j), whiteInWings, irregular);
			vector<int> Wxixj( wingsReachThroughIrregxi.size() + wingsReachThroughIrregxj.size() );
			//Get union of all wings that are reachable to both xi and xj
			vector<int>::iterator itInt;
			
			itInt =set_intersection (wingsReachThroughIrregxi.begin(), wingsReachThroughIrregxi.end(), wingsReachThroughIrregxj.begin(),wingsReachThroughIrregxj.end(), Wxixj.begin());
			
			Wxixj.resize(itInt - Wxixj.begin());
			///////////////////////////////////////////////////////////////////////////////////
			vector<int> n1eWxixj(nv1[j].size()+Wxixj.size());
			
			vector<int>::iterator itInt2;
			
			itInt2 = set_intersection(nv1[j].begin(), nv1[j].end(), Wxixj.begin(), Wxixj.end(), n1eWxixj.begin());
			
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
			
			itInt3 = set_intersection(nv2[j].begin(), nv2[j].end(), Wxixj.begin(), Wxixj.end(), n2eWxixj.begin());
			
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
			
			itInt4 = set_intersection(nv1[i].begin(), nv1[i].end(), Wxixj.begin(), Wxixj.end(), n1xieWxixj.begin());
			
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
			
			itInt5 = set_intersection(nv2[i].begin(), nv2[i].end(), Wxixj.begin(), Wxixj.end(), n2xieWxixj.begin());
			
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

	// Can use a while loop?
	int nv = nv1.size();
	for (int i = 0; i < nv; i++) {
		for (int j = 0; j < nv; j++) {
			for (int l = 0; l < nv1[i].size(); l++) {
				for (int ll = 0; ll < nv2[j].size(); ll++) {
					// Is this a valid way to check that paths P11 and P22 have no irregular vertices?
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
    //1. Copy stable set S into myS (S is intially empty)
    vector<int> myS;
    
    myS.resize(S.size(), '\0');  // '\0' is initial arbitrary value for elements
    copy(S.begin(), S.end(), myS.begin());
    
    //1.1 Obtain white vertices
    vector<int> whiteVertices(n);                      // 0  0  0 ... 0 (n of them)
    vector<int>::iterator it;
    
    sort(myS.begin(),myS.end());
    sort(nodes.begin(),nodes.end()); //Need to sort before taking symmetric difference
    
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
    
    //3. Find maximum weighted white augmenting path between two distinct free vertices
    vector<int> maximumWWAP = maxWeightWhiteAugPath(superFreeVertices, freeVertices, boundedVertices, myS);

    //cout<<"Max weight white augmenting path"<<endl;
    //PrintVector(maximumWWAP);
    
    //4. If no such path exists for S, then stop. S is optimal;
    if (maximumWWAP.size() == 0) {
        
        cout<<"#####################max weight co-2-plex########################:"<<endl;
        PrintVector(myS);
		//Should write myS into a file, along with its weight
        //exit(0);
        
    }
    else{
       //5. S = S delta P and return to (1)
        vector<int> myNewS(n); // 0  0  0 ... 0 (n of them)
        vector<int>::iterator itt;
        
        sort(myS.begin(),myS.end());
        sort(maximumWWAP.begin(), maximumWWAP.end());
        
        itt=std::set_symmetric_difference(myS.begin(),myS.end(), maximumWWAP.begin(), maximumWWAP.end(), myNewS.begin());
        myNewS.resize(itt-myNewS.begin());
        
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