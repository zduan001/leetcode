//
//  main.cpp
//  lol
//
//  Created by Duan, David on 12/20/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <vector>

using namespace std;

int main(int argc, const char * argv[]) {
    
    int i;
    i = 5;
    //cout<<i<<endl;
    vector<int> res;
    for(int i = 150;i> 29;i++)
    {
        res.push_back(i);
        cout<<"I am in the first loop"<<endl;
    }
    
    for(int i = 150;i>29;i++)
    {
        cout<<res[i]<<endl;
        cout<<"I am in the second loop and never has chance to run"<<endl;
    }
    
    // insert code here...
    
    return 81;
}
