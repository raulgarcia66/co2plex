#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

int main(){
    // vector<int> test1 = {2,2,0,0,0};
    // vector<int> test2 = {2,2,0,0,0,0};
    // if (test1 > test2){
    //     cout << "It's sum.";
    // }
    // else {
    //     cout << "It's cardinality if tied at each element.";//winner
    // }

    //int myints[] = {5,6,7,8,9};
    //std::vector<int> myvector (myints,myints+4);
    vector<int> myvector = {5,6,7,11,9};
    vector<int> othervector = {5,4,10,7,6};

    // vector<int> newvector(myvector.size() + othervector.size()); //0 0 ... 0
	vector<int>::iterator it;
	sort(myvector.begin(), myvector.end());
	sort(othervector.begin(), othervector.end());
	it = std::set_difference(myvector.begin(), myvector.end(), othervector.begin(), othervector.end(), myvector.begin());
	myvector.resize(it-myvector.begin());
    
    for (int i = 0; i < myvector.size(); i++) {
        cout << myvector.at(i) << endl;
    }
}