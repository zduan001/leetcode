//
//  main.cpp
//  threadxcode
//
//  Created by Duan, David on 3/11/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <thread>
#include <mutex>

using namespace std;

mutex mu;

void function1(){
    thread::id threadId = this_thread::get_id();
    lock_guard<mutex> guard(mu);
    cout<<"From Function1 "<< threadId<<endl;
}

class A{
public:
    void operator()(string& n){
        lock_guard<mutex> guard(mu);
        cout<<"From Functor " << n << endl;
        n = "new string from functor";
    }
};

int main(int argc, const char * argv[])
{
    //A a;
    string n = "haha";
    thread t1(function1);
    thread t2((A()), ref(n));
    thread::id threadId = this_thread::get_id();
    t2.join();
    cout << "Hello, World! "<< threadId<<endl;
    t1.join();
    cout<<n<<endl;

    return 0;
}

