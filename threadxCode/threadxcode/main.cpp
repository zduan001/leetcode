//
//  main.cpp
//  threadxcode
//
//  Created by Duan, David on 3/11/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <deque>
#include <future>
#include <thread>

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

class LogFile{
    mutex _mu;
    mutex _openMu;
    once_flag _flag;
    ofstream _f;
public:
    LogFile(){
        _f.open("log.txt");
    }
    ~LogFile(){
        _f.close();
    }
    void shared_print(string id, int value){
        //_mu.lock();
        //{
        //    unique_lock<mutex> locker2(_openMu);
        //    if(!_f.is_open()){
        //        _f.open("log.txt");
        //    }
        //}
        
        call_once(_flag, [&](){_f.open("log.txt");});
    
        unique_lock<mutex> guard(mu);
        _f<< "From: "<< id << ": " << value <<endl;
    }
    
};

deque<int> q;
condition_variable cond;

void function_a(){
    int count = 10;
    while(count >0){
        unique_lock<mutex> locker(mu);
        q.push_front(count);
        locker.unlock();
        //cond.notify_one();
        cond.notify_all();
        this_thread::sleep_for(chrono::seconds(1));
        count--;
    }
}

void function_b(){
    int data = 0;
    while(data!= 1){
        unique_lock<mutex> locker(mu);
        //if(!q.empty()){
        cond.wait(locker, [](){return !q.empty();}); // spurious wake.
        data = q.back();
        q.pop_back();
        locker.unlock();
        cout<< " t2 got value from t1 " << data <<endl;

    }
}


int function_factorial(int n){
    int res = 1;
    for(int i = n;i>1;i--){
        res *=i;
    }
    this_thread::sleep_for(chrono::seconds(5));
    return res;
}

int function_factorial_fu(future<int>& f){
    int res = 1;
    int n =f.get();
    for(int i = n;i>1;i--){
        res *=i;
    }
    this_thread::sleep_for(chrono::seconds(5));
    return res;
}


int main(int argc, const char * argv[])
{
/*    //A a;
    string n = "haha";
    thread t1(function1);
    thread t2((A()), ref(n));
    thread::id threadId = this_thread::get_id();
    t2.join();
    cout << "Hello, World! "<< threadId<<endl;
    t1.join();
    cout<<n<<endl;
*/
/*
    thread t1(function_a);
    thread t2(function_b);
    t1.join();
    t2.join();
 */
    
    future<int> fu = async(function_factorial, 4);
    cout<<fu.get()<<endl;
    
    promise<int> p;
    future<int>f = p.get_future();
    fu = async(launch::async, function_factorial_fu, ref(f));
    p.set_value(5);
    cout<<fu.get()<<endl;
    
    
    return 0;
}

