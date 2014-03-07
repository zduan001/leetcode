//
//  container.cpp
//  LeetCode
//
//  Created by Duan, David on 3/4/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include "container.h"
#include <vector>
#include <forward_list>
#include <list>
using namespace std;

container::container(){
    
}

container::~container(){
    
}

void container::testMath()
{
    vector<int> vec;
    vec.push_back(1);
    vec.push_back(2);
    vec.push_back(3);
    
    cout << vec[2]<<endl;
    cout <<vec.at(2)<<endl;
    
    int *p = &vec[1];
    
    vector<int> vec2;
    vec2.swap(vec);
    
    for(vector<int>::iterator it = vec2.begin(); it != vec2.end();it++){
        cout<< *it<<endl;
        cout<<"sizeof it  "<<sizeof(it)<<endl;
        cout<<sizeof(*it)<<endl;
    }
    
    list<int> *myList = new list<int>();
    myList->push_back(1);
    myList->push_front(2);
    
    for(list<int>::iterator it = myList->begin(); it != myList->end(); it++){
        cout<<*it<<endl;
    }
}
