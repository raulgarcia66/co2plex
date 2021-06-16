/*
 *  Algorithm to find Maximum weighted co-2-plex
 *
 *
 *
 *  Created by Cynthia Wood on ?/2014
 *
 */

//#include<assert.h>
#include "modcppstd.hh"

//#include "time.h"


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

Graph::Graph(char* InputFile, char* weightFile){
	
	ifstream f_in, f_win;
	f_in.open(InputFile);
    
	
	f_in>>n;//number of vertices
	f_in>>m;//number of edges
   
    numMatchingProblem = 0;
	
	int *x = new int[m];
	int *y = new int[m];
    
	int z,w;
	
	for(int i = 0; i < m; i++){
		f_in>>x[i]>>y[i]>>z;
	}
	f_in.close();
    
    f_win.open(weightFile);
    for (int i = 0; i < n; i++) {
        f_win>>z>>w;
        weights.push_back(w);
    }
    
    f_win.close();
	
	//Create matrix
	C = new int*[n];
	for( int i = 0; i < n;i++) {
		C[i] = new int[n];
		
	}
	
	for(int i = 0; i < n; i++){
		for(int j = i; j < n; j++) {
			
			C[i][j] = 0;
			C[j][i] = 0;
		}
	}
	
	for(int i = 0; i < m; i++){
		C[x[i]][y[i]] = 1;
		C[y[i]][x[i]] = 1;
	}
	
    //PrintMatrix(C,n);
	delete[] x;
	delete[] y;
}

Graph::~Graph(){
	
	for(int i = 0; i < n; i++){
		delete[] C[i];
	}
	delete[] C;
}


//Only checks the number of neighbors upto 4
int Graph::numberOfBlackNeighbors(int vertex, vector<int> blackVertices){
    
    int numOfN = 0;
    
    for (int i = 0; i < blackVertices.size(); i++) {
        if ( C[vertex][blackVertices.at(i)]) {
            
            numOfN++;
            
        }
        
        if (numOfN == 4) {
            break;
        }
        
    }
    
    return numOfN;
}


void Graph::initMinty(){
    
    //1. S is originally empty
    vector<int> S;
   // All vertices are white on the first step
    
    for (int i = 0; i < n; i++) {
        nodes.push_back(i);
        
    }
    
    
    //PrintVectorVector(white);
    
    Minty(S);
    
}


int Graph::weightOfaPath(vector<int> weightOfBlackPath, vector<int> weightofWhitePath){
    
    int weightWhite = accumulate(weightofWhitePath.begin(),weightofWhitePath.end(),0);
    
    int weightBlack = accumulate(weightOfBlackPath.begin(), weightOfBlackPath.end(), 0);
    
    int weightPath = weightWhite-weightBlack;
    
    //cout<<"weightWhite: "<<weightWhite<<endl;
    //cout<<"weightBlack: "<<weightBlack<<endl;
    
    return weightPath;
}

//Check if a white vertex adjacent to a black vertex violates the co-2-plex property to know if it is FT1 or SFT2
//Return true if violates co-2-plex prop
bool Graph::violateCo2PlexProp( int vertex, vector<int> blackVertices){
    
    int neighbor = -1;
   //Find the neighbor of vertex
    for (int i = 0; i < blackVertices.size(); i++) {
        if (C[vertex][blackVertices.at(i)]) {
            neighbor = blackVertices.at(i);
            break;
        }
    }
    
    //Determine if the black neighbor has a neighbor within the set of black vertices
    if (neighbor > -1) {
        for (int i = 0; i < blackVertices.size(); i++) {
            if (C[neighbor][blackVertices.at(i)]) {
                return true;
            }
        }
    }
    
    return false;
    
}


//Violates co2plex prop, but gives a co2plex whenever the symmetric difference is taken, for this the black vertices adjacent to vertex, must not be adjacent to any other vertex in the black vertex set
//Returns true if co2plex after symm diff
bool Graph::co2PlexAfterSymmDiff( int vertex, vector<int> blackVertices){
    
    int neighbor1 = -1;
    int neighbor2 = -1;
    //Find the neighbor of vertex
    for (int i = 0; i < blackVertices.size(); i++) {
        if (neighbor1 == -1 && C[vertex][blackVertices.at(i)]) {
            neighbor1 = blackVertices.at(i);
        }
        if (neighbor1 > -1 && C[vertex][blackVertices.at(i)]) {
            neighbor2 = blackVertices.at(i);
        }
        
        if (neighbor1 != neighbor2 && neighbor2 > -1 ) {
            break;
        }
    }
    
    if (neighbor1 > -1 && neighbor2 >-1) {
        for (int i = 0; i < blackVertices.size(); i++) {
            if (C[neighbor1][blackVertices.at(i)] || C[neighbor2][blackVertices.at(i)]) {
                return false;
            }
        }
    }
    
    return true;
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




vector<int> Graph::maxWeightWhiteAugPath(vector<int> SFT1, vector<int> SFT2, vector<int> FT1, vector<int> FT2, vector<int> blackVertices, vector<int> BT1, vector<int> BT2, vector<int> BT3){
    sort (blackVertices.begin(),blackVertices.end());
    vector<int> mergeFT12(FT1.size()+FT2.size()+SFT2.size());//Contains FT1, FT2 & SFT2
    vector<int> mergeF(FT1.size()+FT2.size());
    merge (FT1.begin(),FT1.end(),FT2.begin(),FT2.end(),mergeF.begin());
    merge (mergeF.begin(),mergeF.end(),SFT2.begin(),SFT2.end(),mergeFT12.begin());
    //cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
    //PrintVector(mergeFT12);
    //2.1 Generate all white alternating pathsof length 0 or 2
    
    //Save SF vertex with maximum weight only
    int maxWeightl0 = -1;
	int maxWeightln = -1;
    vector<int> maxWWAPl0;
	vector<int> maxWWAPln;
    for (int i = 0; i < SFT1.size(); i++) {
        int tempMax = weights.at(SFT1.at(i));
        if (tempMax > maxWeightl0) {
            maxWWAPl0.clear();
            maxWWAPl0.push_back(SFT1.at(i));
            maxWeightl0 = tempMax;
            
        }
    }
    
    
    //SFT2 vertices can only be added to the set S when there is no WAP of positive weight, otherwise they may bound vertices that should not be bounded, and we we add them if necessary later in this funtion
    
    
    //Now we will generate WAP of length 2 and save it only if its weight is greater than the current max weight
    vector<int> blackPairs;//Will always be even. An odd entry will always be paired with the even entry that preceeds it
    vector<int> sortBlackPairs;
    for (int i = 0; i < blackVertices.size(); i++) {
        for (int j = i+1; j < blackVertices.size(); j++) {
            if (C[blackVertices.at(i)][blackVertices.at(j)]) {
                blackPairs.push_back(blackVertices.at(i));
                blackPairs.push_back(blackVertices.at(j));
                sortBlackPairs.push_back(blackVertices.at(i));
                sortBlackPairs.push_back(blackVertices.at(j));
                break;// Since it is only adjacent to one other black vertex
            }
        }
    }
    
    //Black pairs are not contained in white alternating paths of length 2
    vector<int> blackSingle(n);                      // 0  0  0 ... 0 (n of them)
    vector<int>::iterator it;
    
	//sort (blackVertices.begin(), blackVertices.end());you sort them at the beginning of the function
    sort (sortBlackPairs.begin(),sortBlackPairs.end());
    
    it=std::set_symmetric_difference (blackVertices.begin(),blackVertices.end(), sortBlackPairs.begin(),sortBlackPairs.end(), blackSingle.begin());
    //  blackvertices ... 0  0  0  0
    blackSingle.resize(it-blackSingle.begin());//only black vertices
    vector<int> weightWhite;
    vector<int> weightBlack;
    for (int i = 0; i < blackSingle.size(); i++) {
        int adjacent1 = -1;
        int adjacent2 = -1;
        weightBlack.push_back(weights.at(blackSingle.at(i)));
        for (int j = 0; j < mergeFT12.size(); j++) {
            if (adjacent1 == -1 && C[blackSingle.at(i)][mergeFT12.at(j)]) {
                adjacent1 = mergeFT12.at(j);
                for (int k = j+1; k < mergeFT12.size(); k++) {
                    if (C[blackSingle.at(i)][mergeFT12.at(k)]  && !C[mergeFT12.at(j)][mergeFT12.at(k)] ) {
                        adjacent2 = mergeFT12.at(k);
                        weightWhite.push_back(weights.at(mergeFT12.at(j)));
                        weightWhite.push_back(weights.at(mergeFT12.at(k)));
                        
                        int tempw = weightOfaPath(weightBlack, weightWhite);
                        if (tempw > maxWeightl0) {
                            maxWWAPl0.clear();
                            maxWWAPl0.push_back(mergeFT12.at(j));
                            maxWWAPl0.push_back(blackSingle.at(i));
                            maxWWAPl0.push_back(mergeFT12.at(k));
                            maxWeightl0 = tempw;
                        }
                        weightWhite.clear();
                    }
                    
                }
            }
        }
        weightBlack.clear();
    }
    
    //2.2 For each pair of nonadjacent free vertices and super free vertices of type 2 a and b do
    
    //For each pair of nonadjacent Free VERTICES a & b do
    for (int a = 0; a < mergeFT12.size(); a++) {
        for (int b = a+1; b < mergeFT12.size(); b++) {
            if (!C[mergeFT12.at(a)][mergeFT12.at(b)]) {
                //Let xa and xb be the black vertices andjacent to a and b respectively (at most 2)
                vector<int> xa;// black vertex or vertices adjacent to a
                vector<int> xb;// black vertex or vertices adjacent to b
                for (int i = 0; i < blackVertices.size(); i++) {
                    if (C[a][blackVertices.at(i)]) {
                        xa.push_back(blackVertices.at(i));
                    }
                    
                    if (C[b][blackVertices.at(i)]) {
                        xb.push_back(blackVertices.at(i));
                    }
                }
                
                bool xEqual = false;
                
                for (int i = 0; i < xa.size(); i++) {
                    for (int j = 0; j < xb.size(); j++) {
                        //Want ot make sure the black vertices to wich they are joined are distinct
                        if (xa.at(i) == xb.at(j)) {
                            xEqual = true;
                            break;
                        }
                    }
                    
                    if (xEqual) {
                        break;
                    }
                }
                
                //if xa != xb
                if (!xEqual) {
				
					vector< vector<int> > EdmondsG;
                    //find MWWAP
					//N = 2*rbsN+2
					int N = 0;
                    EdmondsG = Edmonds(mergeFT12.at(a), mergeFT12.at(b), xa, xb, blackSingle, blackPairs, BT1, BT2, BT3,N);
					
					sort(EdmondsG.begin(), EdmondsG.end());
                    
                   // cout<<"LINE 460\n";
					
					if ( N > EdmondsG.size()) {
						double matzeit = 0.0;
						double genzeit = 0.0;
						int ncount, ecount;
						//long l;
						time_t l; // my fix
						//include the right header for this one otherwise it gives you a warning( "time.h")
						seed = time (&l); //stores the current time in l
						
						
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
						
						//The maping is a pair entry(i)-> entries 2i and 2i+1 weight
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
    
    if (maxWeightl0 < 0) {
        for (int i = 0; i < SFT2.size(); i++) {
            int tempMax2 = weights.at(SFT2.at(i));
            if (tempMax2 > maxWeightl0) {
                maxWWAPl0.clear();
                maxWWAPl0.push_back(SFT2.at(i));
                maxWeightl0 = tempMax2;
                
            }
        }
    }
    
	if (maxWWAPln > maxWWAPl0) {
		return maxWWAPln;
	}
	
    return maxWWAPl0;
}

vector<int> Graph::edmondsCorrection1(vector<int> regNeib, vector<int> whiteWings, vector<int> irreg, vector< vector<int> >irregPairs){
	vector<int> whiteReachableWingsThroughIrreg;
	
	
	vector< vector<int> > blackLevel;
	for (int i = 0; i < regNeib.size(); i++) {
		for (int j = 0; j < irreg.size(); j++) {
			if (C[regNeib.at(i)][irreg.at(j)]) {
				whiteReachableWingsThroughIrreg. push_back(regNeib.at(i));
				vector<int> tempBL;
				
				tempBL. push_back(irreg.at(j));
				blackLevel. push_back(tempBL);
				tempBL.clear();
				//erase entry j from possible black vertices
				irreg.erase(irreg.begin() + j);
				//When you remove you decrease the size of the vector
				j--;
				
			}
		}
		
		vector<int> tempNbr;
		tempNbr. push_back(regNeib.at(i));
		for (int j = 0; j < irregPairs.size(); j++) {
			 
			if (adjacent(tempNbr, irregPairs.at(j))) {
			
				whiteReachableWingsThroughIrreg. push_back(regNeib.at(i));
				vector<int> tempBL2;
				
				tempBL2. push_back(irreg.at(j));
				blackLevel. push_back(tempBL2);
				tempBL2.clear();
				
				//erase entry j from possible black pairs
				irregPairs.erase(irregPairs.begin() + j);
				//When you remove you decrease the size of the vector
				j--;
				
			}
		}
		tempNbr.clear();
	}
	
	bool moreLevels = true;

	vector<int> whiteLevel;
	while (moreLevels) {
		
		for (int i = 0; i < blackLevel.size() ; i++) {
			for (int j = 0; j < whiteWings.size(); j++) {
				vector<int> tWhite;
				tWhite. push_back(whiteWings.at(j));
				
				if (adjacent(blackLevel.at(i), tWhite)) {
					whiteLevel. push_back(whiteWings.at(j));
					
					//erase entry j from possible white vertices
					whiteWings.erase(whiteWings.begin() + j);
					//When you remove you decrease the size of the vector
					j--;
				}
				tWhite.clear();
			}	
		}
		
		blackLevel. clear();
		
		for (int i = 0; i < whiteLevel.size(); i++) {
			for (int j = 0; j < irreg.size(); j++) {
				if (C[whiteLevel.at(i)][irreg.at(j)]) {
					whiteReachableWingsThroughIrreg. push_back(whiteLevel.at(i));
					vector<int> tempBL2;
					
					tempBL2. push_back(irreg.at(j));
					blackLevel. push_back(tempBL2);
					tempBL2.clear();
					//erase entry j from possible black vertices
					irreg.erase(irreg.begin() + j);
					//When you remove you decrease the size of the vector
					j--;
					
				}
			}
			
			
			vector<int> tempWht;
			tempWht. push_back(whiteLevel.at(i));
			for (int j = 0; j < irregPairs.size(); j++) {
				 
				if (adjacent(tempWht, irregPairs.at(j))) {
				
					whiteReachableWingsThroughIrreg. push_back(whiteLevel.at(i));
					vector<int> tempBL2;
					
				
					blackLevel. push_back(irregPairs.at(j));
					
					
					//erase entry j from possible black pairs
					irregPairs.erase(irregPairs.begin() + j);
					//When you remove you decrease the size of the vector
					j--;
					
				}
			}
			tempWht.clear();
		}
		whiteLevel. clear();
		
		if (blackLevel.size() == 0) {
			moreLevels = false;
		}
	}
	
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

vector<bool> Graph::iwapBF(vector<int> &weight, vector<int> &whiteReachIrreg, vector<int> whiteSource, vector<int> whiteSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices, vector< vector<int> > irregularPairs){

	// weight must be empty and will contain the weight of IWAPS 
	//if the first IWAP between the source and a whiteSink is whiteSink.at(2)
	//then the first entry of weight will be the entry of such a path
	vector<bool> IWAP( whiteSinks.size(), false);
	vector<int> vertexSet;
	//Construct a DAG to solve the problem
	vector< vector<int> > level;
	level.push_back(whiteSource);
	vertexSet.push_back(whiteSource.at(0));//We're doing contraction
	bool moreLevelsNeeded = true;
	bool white = false;//search through white vertices if true
	
	vector<int> edgeSet;
	vector<int> edgeWeights;
	
	//only save two vertices even if they are a pair, but do consider the adjacencies of the whole thing
	
	vector< vector<int> > irregularVP;
	for (int i = 0; i < irregularVertices.size(); i++) {
		vector<int> t;
		
		t.push_back(irregularVertices.at(i));
		
		irregularVP.push_back(t);
		t.clear();
	}
	
	for (int i = 0 ; i < irregularPairs.size(); i++) {
		irregularVP.push_back(irregularPairs.at(i));
	}
	
	vector<int> whiteIndices;
	while (moreLevelsNeeded) {
		if (white) {
		
			int counter = 0;
			vector< vector<int> > tempLevel;
			for (int i = 0 ; i < level.size(); i++) {
				bool inWhiteLevel = false;
				
				for (int j = 0; j < whiteVerticesInWings.size(); j++) {
					vector<int> tempW;
					tempW.push_back(whiteVerticesInWings.at(j));
					if (adjacent(level[i],tempW)) {
						counter++;
						inWhiteLevel = true;
						edgeSet.push_back(level[i].at(0));
						edgeSet.push_back(tempW.at(0));
						vertexSet.push_back(tempW.at(0));
						edgeWeights.push_back(-weights.at(tempW.at(0)));
						
					}
					
					
					if (inWhiteLevel) {
					
						for (int k = 0; k < whiteSinks.size(); k++) {
							if (whiteSinks.at(k) == whiteVerticesInWings.at(j)) {
								IWAP.at(k) = true;
								whiteIndices.push_back(vertexSet.size()-1);
								break;
							}
						}
						
						//erase entry j from possible white vertices
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
					if (adjacent(level[i], irregularVP[j])) {
						counter2++;
						inBlackLevel = true;
						edgeSet.push_back(level[i].at(0));
						edgeSet.push_back(irregularVP[j].at(0));
						vertexSet.push_back(irregularVP[j].at(0));
						int w2 = 0;
						for (int jj = 0; jj < irregularVP[j].size(); jj++) {
							w2 -= weights.at(irregularVP[j].at(jj));
						}
						edgeWeights.push_back(w2);
					}
					
					if (inBlackLevel) {
						tLevel.push_back(irregularVP.at(j));
					
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
				
				predecessor.at(entry2) = edgeSet.at(entry1);
			}
		}
	}
	
	// Step 3: check for negative-weight cycles
	//NO need to check this since we have a DAG
			   

	//Used saved index of white vertices to figure out the weight of the IWAPs
	
	for (int i = 0; i < whiteIndices.size(); i++) {
	
		weight.push_back(distance.at(whiteIndices.at(i)));
	}
	
	return IWAP;
}


vector<int> Graph::bellmanFordVariant(vector<int> regularSource, vector< vector<int> > regularSinks, vector<int> whiteVerticesInWings, vector<int> irregularVertices, vector< vector<int> >irregularPairs) {
	
	vector<int> vertexSet;
	//Construct a DAG to solve the problem
	vector< vector<int> > level;
	level.push_back(regularSource);
	vertexSet.push_back(regularSource.at(0));//We're doing contraction
	bool moreLevelsNeeded = true;
	bool white = true;//search through white vertices
	
	vector<int> edgeSet;
	vector<int> edgeWeights;
	
	//only save two vertices even if they are a pair, but do consider the adjacencies of the whole thing
	
	vector< vector<int> > irregularVP;
	for (int i = 0; i < irregularVertices.size(); i++) {
		vector<int> t;
		
		t.push_back(irregularVertices.at(i));
		
		irregularVP.push_back(t);
		t.clear();
	}
	
	for (int i = 0 ; i < irregularPairs.size(); i++) {
		irregularVP.push_back(irregularPairs.at(i));
	}
	
	vector<int> whiteIndices;
	while (moreLevelsNeeded) {
		if (white) {
		
			int counter = 0;
			vector< vector<int> > tempLevel;
			for (int i = 0 ; i < level.size(); i++) {
				bool inWhiteLevel = false;
				
				for (int j = 0; j < whiteVerticesInWings.size(); j++) {
					vector<int> tempW;
					tempW.push_back(whiteVerticesInWings.at(j));
					if (adjacent(level[i],tempW)) {
						counter++;
						inWhiteLevel = true;
						edgeSet.push_back(level[i].at(0));
						edgeSet.push_back(tempW.at(0));
						whiteIndices.push_back(vertexSet.size());
						vertexSet.push_back(tempW.at(0));
						edgeWeights.push_back(-weights.at(tempW.at(0)));
						
					}
					
					
					if (inWhiteLevel) {
						//erase entry j from possible white vertices
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
					if (adjacent(level[i], irregularVP[j])) {
						counter2++;
						inBlackLevel = true;
						edgeSet.push_back(level[i].at(0));
						edgeSet.push_back(irregularVP[j].at(0));
						vertexSet.push_back(irregularVP[j].at(0));
						int w2 = 0;
						for (int jj = 0; jj < irregularVP[j].size(); jj++) {
							w2 -= weights.at(irregularVP[j].at(jj));
						}
						edgeWeights.push_back(w2);
					}
					
					if (inBlackLevel) {
						tLevel.push_back(irregularVP.at(j));
					
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
		if (vertexSet.at(i) == regularSource.at(0)) {
			distance.push_back(0);
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
		for (int j = 0; j < edgeSet.size(); j+=2) {
			if (distance.at(j) + edgeWeights.at(j/2) < distance.at(j+1)) {
				distance.at(j+1) = distance.at(j) + edgeWeights.at(j);
				
				predecessor.at(j+1) = edgeSet.at(j);
			}
		}
	}
	
	// Step 3: check for negative-weight cycles
	//NO need to check this since we have a DAG
			   

	
	vector<int> IWAPweights;
	//Used saved index of white vertices to figure out which regularSinks must be joined
	
	for (int i = 0; i < whiteIndices.size(); i++) {
		vector<int> tempA;
		tempA.push_back(vertexSet.at(i));
		for (int j = 0; j < regularSinks.size(); j++) {
			if (adjacent(tempA, regularSinks.at(j))) {
				
				IWAPweights.push_back(-distance.at(i));
				
				
				
			}
			else {
				IWAPweights.push_back(-1);
			}

		}
		
		tempA.clear();
	}
	
	//returns the max weight of the IWAP found between two distinct regular vertices
	//-1 if no path to the given source
	// the weight of the IWAP otherwise
	return IWAPweights;
}


vector < vector<int> > Graph::Edmonds(int freeVertexA, int freeVertexB, vector<int> xa, vector<int> xb, vector<int> blackSetSingle,vector<int> blackSetPairs, vector<int> boundedVerticesT1, vector<int> boundedVerticesT2, vector<int> boundedVerticesT3, int &N){
    
    vector< vector<int> > edmondsG;
   
    sort(xa.begin(),xa.end());
    sort(xb.begin(),xb.end());
    sort(blackSetSingle.begin(),blackSetSingle.end());//Don't really need to do this anymore, but leave it for now
	int weightxA = 0;
	int weightxB = 0;
	
	for (int i = 0 ; i < xa.size(); i++) {
		weightxA += weights.at(xa.at(i));
	}
	
	for (int i = 0 ; i < xb.size(); i++) {
		weightxB += weights.at(xb.at(i));
	}
	
    /*We want to construct the reduced basic structure
     1.1 Ignore all the SF objects (just don't pass them to the function)
     1.2 Ignore F  vertices except a and b, we should also ignore SFT2 ( don't have to pass them either)
     1.3 Ignore all white vertices/objects adjacent to a and b, since they will never appear in an alternating path between
     */
    
    vector<int> whiteverticesNonAdjacentAorB( boundedVerticesT1.size() + boundedVerticesT2.size() + boundedVerticesT3.size());
    
    vector<int> whitev2(boundedVerticesT1.size()+ boundedVerticesT2.size());
    merge (boundedVerticesT1.begin(),boundedVerticesT1.end(),boundedVerticesT2.begin(),boundedVerticesT2.end(),whitev2.begin());
    merge (whitev2.begin(),whitev2.end(),boundedVerticesT3.begin(),boundedVerticesT3.end(),whiteverticesNonAdjacentAorB.begin());
    
    sort(whiteverticesNonAdjacentAorB.begin(), whiteverticesNonAdjacentAorB.end());
    
    //cout<<"#########################################################"<<endl;
    //PrintVector(whiteverticesNonAdjacentAorB);
   
    
    for (int i = 0; i < whiteverticesNonAdjacentAorB.size(); i++) {
        //Remove A & B from free, since we don't want repetitions
        
        if (C[whiteverticesNonAdjacentAorB.at(i)][freeVertexA] ||C[whiteverticesNonAdjacentAorB.at(i)][freeVertexB] || whiteverticesNonAdjacentAorB.at(i) == freeVertexA || whiteverticesNonAdjacentAorB.at(i) == freeVertexB) {
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
    vector<int> whiteRBS;
    
    whiteRBS.resize( whiteverticesNonAdjacentAorB.size(), '\0');
    copy(whiteverticesNonAdjacentAorB.begin(), whiteverticesNonAdjacentAorB.end(), whiteRBS.begin());
    whiteRBS.push_back(freeVertexA);
    whiteRBS.push_back(freeVertexB);
    
    vector<int> blackSingleRBS;
    
    blackSingleRBS.resize( blackSetSingle.size(), '\0');
    copy(blackSetSingle.begin(), blackSetSingle.end(), blackSingleRBS.begin());
    
    vector< vector<int> > blackPairsRBS;
    vector<int> weightsBlackPairs;
    
    for (int i = 0; i < blackPairsRBS.size(); i++) {
        //The original black set comes sorted, which implies entry[i] < entry[i+1]
        if (i%2 == 0) {
            vector<int> temp;
            temp.push_back(blackSetPairs.at(i));
            temp.push_back(blackSetPairs.at(i+1));
            sort(temp.begin(),temp.end());
            
            if (temp != xa && temp != xb) {
                blackPairsRBS.push_back(temp);
                int weight = 0;
                weight = weights.at(blackSetPairs.at(i)) + weights.at(blackSetPairs.at(i+1)) ;
                weightsBlackPairs.push_back(weight);
            }
            temp.clear();
        }
    }
    
    
    /***********************CLASIFICATION OF BLACK VERTICES*************************/
    
    // A nonempty set of all bounded vertices, which are adjacent to the same two black objects x and y is called a wing
    // First figure out which black vertices and objects are adjacent to bounded vertices
    
    int sinPairSize = blackSingleRBS.size() + blackPairsRBS.size();
    vector< int > numWings(sinPairSize); //keeps a counter to the number of wings each black vertex or pair is adjacent. They are 0 originally.
    
    vector< vector<int> > enumerateWings;//Keeps track of tip of the wings as pairs. Its size also gives how many different wings we have
    
	vector<int> whiteInWings;
	
    for (int i = 0; i < blackSingleRBS.size(); i++) {
        int numi = numWings.at(i);
        for (int j = i+1; j < blackSingleRBS.size(); j++) {
            
            int numj = numWings.at(j);
            for (int k = 0; k < whiteverticesNonAdjacentAorB.size(); k++) {
                if (C[whiteverticesNonAdjacentAorB.at(k)][blackSingleRBS.at(j)] && C[whiteverticesNonAdjacentAorB.at(k)][blackSingleRBS.at(i)]) {
                    numi++;
                    numj++;
					
                    whiteInWings.push_back(whiteverticesNonAdjacentAorB.at(k));
                    vector<int> tempWing;
                    tempWing.push_back(blackSingleRBS.at(i));
                    tempWing.push_back(blackSingleRBS.at(j));
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
    
    for (int i = 0; i < blackSingleRBS.size(); i++) {
        int numi = numWings.at(i);
        for (int j = 0; j < blackPairsRBS.size(); j++) {
            
            int numj = numWings.at(blackSingleRBS.size() + j);
            for (int k = 0; k < whiteverticesNonAdjacentAorB.size(); k++) {
                if (C[whiteverticesNonAdjacentAorB.at(k)][blackSingleRBS.at(j)] && (C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[j].at(0)] || C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[j].at(1)])) {
                    numi++;
                    numj++;
					whiteInWings.push_back(whiteverticesNonAdjacentAorB.at(k));
                    vector<int> tempWing1;
                    tempWing1.push_back(blackSingleRBS.at(i));
                    sort(blackPairsRBS[j].begin(), blackPairsRBS[j].end());
                    tempWing1.push_back(blackPairsRBS[j].at(0));
                    tempWing1.push_back(blackPairsRBS[j].at(1));
                    enumerateWings.push_back(tempWing1);
                    tempWing1.clear();
                    
                    
                    
                    break;//as long as the set = wing is nonempty then we know a particular vertex is part of a wing
                }
            }
            numWings.at(blackSingleRBS.size() + j) = numj;
        }
        
        numWings.at(i) = numi;
    }
    
    
    for (int i = 0; i < blackPairsRBS.size(); i++) {
        int numi = numWings.at(blackSingleRBS.size() + i);
		sort(blackPairsRBS[i].begin(), blackPairsRBS[i].end());
        for (int j = i+1; j < blackPairsRBS.size(); j++) {
            
            int numj = numWings.at(blackSingleRBS.size() + j);
            for (int k = 0; k < whiteverticesNonAdjacentAorB.size(); k++) {
                if ((C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[i].at(0)] || C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[i].at(1)]) && (C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[j].at(0)] || C[whiteverticesNonAdjacentAorB.at(k)][blackPairsRBS[j].at(1)])) {
                    numi++;
                    numj++;
					
					whiteInWings.push_back(whiteverticesNonAdjacentAorB.at(k));
                    vector<int> tempWing2;
                    
					sort(blackPairsRBS[j].begin(), blackPairsRBS[j].end());
                    tempWing2.push_back(blackPairsRBS[i].at(0));
                    tempWing2.push_back(blackPairsRBS[i].at(1));
                    tempWing2.push_back(blackPairsRBS[j].at(0));
                    tempWing2.push_back(blackPairsRBS[j].at(1));
                    enumerateWings.push_back(tempWing2);
                    tempWing2.clear();
                    break;//as long as the set = wing is nonempty then we know a particular vertex is part of a wing
                }
            }
            numWings.at(blackSingleRBS.size() + j) = numj;
        }
        
        numWings.at(blackSingleRBS.size() + i)= numi;
    }
    
	sort(enumerateWings.begin(), enumerateWings.end());
	sort(whiteInWings.begin(), whiteInWings.end());
	vector<int>::iterator itWiW;
	itWiW = unique( whiteInWings.begin(), whiteInWings.end());
	whiteInWings.resize( std::distance(whiteInWings.begin(),itWiW));
	
	
    //cout<<"***************************************************************"<<endl;
    //PrintVector(numWings);
    
    //Regular I: vectors xa and xb. Make sure you don't count this in black pairs
    //Regular II: A black vertex adjacent to 3 or more wings
    vector<int> regularII;
    vector< vector<int> > regularIIPairs;
    //Irregular: Adjacent to exactly 2 wings
    //useless otherwise
    vector<int> irregular;
    vector< vector<int> > irregularPairs;
	int irregularPairsSize =0;
	
    for (int i = 0; i < blackSingleRBS.size(); i++) {
        if (numWings.at(i) > 3) {
			//cout << "regular"<<endl;
            regularII.push_back(blackSingleRBS.at(i));
        }
        if (numWings.at(i) == 2) {
			//cout << "irregular"<<endl;
            irregular.push_back(blackSingleRBS.at(i));
        }
        
        //A useless vertex along with all adjoining white vertices can be deleted
        if (numWings.at(i) <2) {
            
            for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
                if (C[blackSingleRBS.at(i)][whiteverticesNonAdjacentAorB.at(j)]) {
				
					vector<int>::iterator itw;
					itw = find (whiteInWings.begin(), whiteInWings.end(), whiteverticesNonAdjacentAorB.at(j));
            
					if (itw != whiteInWings.end()) {
						whiteInWings.erase(itw);
					}
					
					whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + j);
                    j--;
	
					
                }
            }
            blackSingleRBS.erase(blackSingleRBS.begin() + i);
            numWings.erase(numWings.begin() + i);
            i--;
        }
    }
    
    int bsp = blackSingleRBS.size() + blackPairsRBS.size();
    for (int i = blackSingleRBS.size(); i < bsp; i++) {
        
        if (numWings.at(i) > 3) {
			cout << "regular pair"<<endl;
            regularIIPairs.push_back(blackPairsRBS[i - blackSingleRBS.size()]);
            
        }
        if (numWings.at(i) == 2) {
			cout << "irregular pair"<<endl;
            irregularPairs.push_back(blackPairsRBS[i - blackSingleRBS.size()]);
			irregularPairsSize++;
        }
        
        //A useless object along with all adjoining white vertices can be deleted
        // don't include them
        if (numWings.at(i) <2) {
            
            for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
                if (C[blackPairsRBS[i - blackSingleRBS.size()].at(0)][whiteverticesNonAdjacentAorB.at(j)] || C[blackPairsRBS[i - blackSingleRBS.size()].at(1)][whiteverticesNonAdjacentAorB.at(j)]) {
                    vector<int>::iterator itw2;
					itw2 = find (whiteInWings.begin(), whiteInWings.end(), whiteverticesNonAdjacentAorB.at(j));
            
					if (itw2 != whiteInWings.end()) {
						whiteInWings.erase(itw2);
					}
	
					whiteverticesNonAdjacentAorB.erase(whiteverticesNonAdjacentAorB.begin() + j);
                    j--;
                }
            }
            
            //just keep them empty here, then you don't have to decrease i, change this to remove it for real 10.15.14
            blackPairsRBS[i - blackSingleRBS.size()].clear();
            blackPairsRBS.erase(blackPairsRBS.begin() + i - blackSingleRBS.size());
            i--;
        }
       
    }
    
    //NEW  white RBS whiteverticesNonAdjacentAorB U freevertexA U freeVertexB instead of white RBS
    
    //N1(xa) = freevertexA and N1(xb) = frevertexB
    // for a regular vertex of the first kind put the free vertex into one node and the rest of its neighbors into the other class

    vector<int> n2xa;
    vector<int> n2xb;
    for(int i = 0 ; i < whiteverticesNonAdjacentAorB.size(); i++)
    {
		
		vector<int> tempw;
		tempw.push_back(whiteverticesNonAdjacentAorB.at(i));
        if (xa.size() == 2) {
			
			if (adjacent(xa, tempw)){
				//cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
				n2xa.push_back(whiteverticesNonAdjacentAorB.at(i));
				//PrintVector(n2xa);
			}
        }
		
		if (xb.size() == 2) {
			//cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			if (adjacent(xb, tempw))n2xb.push_back(whiteverticesNonAdjacentAorB.at(i));
    
        }
		
		tempw.clear();
        if (xa.size() == 1) {
            if(C[xa.at(0)][whiteverticesNonAdjacentAorB.at(i)]) n2xa.push_back(whiteverticesNonAdjacentAorB.at(i));
            
        }
		if (xb.size() == 1) {
            if(C[xb.at(0)][whiteverticesNonAdjacentAorB.at(i)]) n2xb.push_back(whiteverticesNonAdjacentAorB.at(i));
            
        }
		
    }
    
    
    //For any regularII black vertex or pair v  or {v,w} get the neighborhood of v or {v,w}
    int regIIsP = regularII.size() + regularIIPairs.size();
    vector< vector<int> > neighborsV;//These are the neighbors of regularII vertices and pairs
    
    vector<int> tempRegII;
    for (int i = 0;  i < regIIsP; i++) {
        
        if (i < regularII.size()) {
            
            for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
                if (C[regularII.at(i)][whiteverticesNonAdjacentAorB.at(j)]) {
                    tempRegII.push_back(whiteverticesNonAdjacentAorB.at(j));
                }
            }
            
            neighborsV.push_back(tempRegII);
            tempRegII.clear();
           
        }
        else{
            
            for (int j = 0; j < whiteverticesNonAdjacentAorB.size(); j++) {
                if (C[regularIIPairs[i - regularII.size()].at(0)][whiteverticesNonAdjacentAorB.at(j)] || C[regularIIPairs[i - regularII.size()].at(1)][whiteverticesNonAdjacentAorB.at(j)]) {
                    tempRegII.push_back(whiteverticesNonAdjacentAorB.at(j));
                }
            }
            
            neighborsV.push_back(tempRegII);
            tempRegII.clear();
            
        }
        
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
    
    myS.resize( S.size(), '\0');
    copy(S.begin(), S.end(), myS.begin());
    
    //1.1 Obtain white vertices
    vector<int> whiteVertices(n);                      // 0  0  0 ... 0 (n of them)
    vector<int>::iterator it;
    
    sort (myS.begin(),myS.end());
    sort (nodes.begin(),nodes.end());   // Sort for consistency
    
    it=std::set_symmetric_difference (myS.begin(),myS.end(), nodes.begin(),nodes.end(), whiteVertices.begin());
    //  whitevertices ... 0  0  0  0
    whiteVertices.resize(it-whiteVertices.begin());//only white vertices
    
    //cout << "The symmetric difference of S and nodes has " << (whiteVertices.size()) << " elements:\n";
    //for (it=whiteVertices.begin(); it!=whiteVertices.end(); ++it)
      //  cout << ' ' << *it;
    //cout << '\n';
    
    //2. Clasify the white vertices based off the black vertices
    
    //A white vertex on a {claw, bull}-free graph is adjacent to at most 4 vertices
    vector<int> superFreeVerticesT1,superFreeVerticesT2, freeVerticesT1, freeVerticesT2, boundedVerticesT1, boundedVerticesT2, boundedVerticesT3;
    
    for (int i = 0; i < whiteVertices.size(); i++) {
        
        int numBlack = numberOfBlackNeighbors(whiteVertices.at(i), myS);
        if ( numBlack == 0) {
            superFreeVerticesT1.push_back(whiteVertices.at(i));
        }
        //SFT2 or FT1
        if (numBlack == 1 && !violateCo2PlexProp(whiteVertices.at(i), myS)) {
            superFreeVerticesT2.push_back(whiteVertices.at(i));
            
        }
        if (numBlack == 1 && violateCo2PlexProp(whiteVertices.at(i), myS)) {
            freeVerticesT1.push_back(whiteVertices.at(i));
        }
        //FT2 or BT1
        if (numBlack == 2) {
            
            bool co2pAfSymmDiff = co2PlexAfterSymmDiff(whiteVertices.at(i), myS);
            
            if (co2pAfSymmDiff) {
                freeVerticesT2.push_back(whiteVertices.at(i));
            }
            else{
                boundedVerticesT1.push_back(whiteVertices.at(i));
            }
            
        }
   
                                     
                                     
        if (numBlack == 3) {
          
            boundedVerticesT2.push_back(whiteVertices.at(i));
        }
        if (numBlack == 4) {
            
            boundedVerticesT3.push_back(whiteVertices.at(i));
        }
    }
    
    //3. Find maximum weighted white augmenting path between two distinct free vertices
    vector<int> maximumWWAP = maxWeightWhiteAugPath(superFreeVerticesT1, superFreeVerticesT2, freeVerticesT1, freeVerticesT2, myS, boundedVerticesT1, boundedVerticesT2, boundedVerticesT3);
    
    //cout<<"Max weight white augmenting path"<<endl;
    //PrintVector(maximumWWAP);
    
    
    //4. If no such path exists for S, the stop. S is optimal;
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
	
    if(argc != 3) {
        
        cerr<<" USAGE: ./executable graphAdjacency graphWeights"<<endl;
        exit(0);
        
    }
    
    clock_t t1,t2;
    t1=clock();
    
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
    
    
    ostringstream fname;
	fname << "verticeCo2PlexResults.txt";
	//fname << ".txt";
	ofstream f_out;
	f_out.open( fname.str().c_str(),ofstream::app);
    
    t2=clock();
    float diff (((float)t2-(float)t1)/CLOCKS_PER_SEC);
    
    
	f_out<<"time for: "<< argv[1]  << " was: "<< diff <<" s"<<endl;
	
    
    //cout<<"time for: "<< argv[1]  << " was: "<< diff <<" s"<< endl;
    
    return 0;	
	
}

