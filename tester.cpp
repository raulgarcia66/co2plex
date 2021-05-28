#include <iostream>
#include <vector>

using namespace std;

int main(){
    vector<int> test1 = {2,2,0,0,0};
    vector<int> test2 = {2,2,0,0,0,0};
    if (test1 > test2){
        cout << "It's sum.";
    }
    else {
        cout << "It's cardinality if tied at each element.";//winner
    }
}