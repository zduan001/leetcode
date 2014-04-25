//
//  main.cpp
//  lettcodeII
//
//  Created by Duan, David on 4/24/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>

using namespace::std;

struct ListNode {
    int val;
    ListNode *next;
    ListNode(int x) : val(x), next(NULL) {}
};

int findKthElement(int A[], int m, int B[], int n, int k)
{
    if(m>n)
        return findKthElement(B, n, A,m, k);
    if(m== 0)
        return B[k-1];
    if(k == 1)
        return min(A[0], B[0]);
    
    int ia = min(m, k/2);
    int ib = k-ia;
    
    if(B[ib-1] > A[ia-1])
    {
        return findKthElement( A + ia, m-ia, B, n, k-ia);
    }
    else if(B[ib-1] < A[ia-1])
    {
        return findKthElement( A, m, B+ib, n-ib, k- ib);
    }
    else
    {
        return A[ia-1];
    }
}

double findMedianSortedArrays(int A[], int m, int B[], int n) {
    if((n+m)%2 == 0){
        return (findKthElement(A, m, B, n, (n+m)/2) +
                findKthElement(A, m, B, n, (n+m)/2 +1)) /2.0;
    }else{
        return (findKthElement(A, m, B, n, (n+m)/2+1));
    }
}

// assume there exactly one solution
vector<int> twoSum(vector<int> &n, int target) {
    unordered_map<int,int> tracker;
    vector<int> res;
    for(int i = 0;i< (int)n.size();i++){
        tracker.insert(make_pair(n.at(i), i));
    }
    
    for(int i = 0;i<(int)n.size();i++){
        auto it = tracker.find(target-n.at(i));
        if(it != tracker.end())
        {
            if(i != it->second)
            {
                res.push_back(i+1);
                res.push_back(it->second+1);
                break;
            }

        }
    }
    return res;
}

int lengthOfLongestSubstring(string s) {
    if(s.length() == 0) return 0;
    int max = 1;
    int length = 1;
    int start = 0;
    int end = 1;
    bool tracker[256];
    for(int i = 0;i< 256;i++){
        tracker[i] = false;
    }
    tracker[s.at(0)] = true;
    
    while(end != s.length()){
        if(!tracker[s.at(end)])
        {
            length++;
            max = max > length? max : length;
            tracker[s.at(end++)] = true;
        } else {
            while(s.at(start) != s.at(end)){
                tracker[s.at(start++)] = false;
                length--;
            }
            start++;
            end++;
        }
    }
    return max;
}



int main(int argc, const char * argv[])
{
    /*
    int A[] = {100};
    int B[] = {101};
    
    double res = findMedianSortedArrays(A, 1, B, 1);
    cout<<res;
    */
    //vector<int> input = {3,2,4};
    //vector<int> res = twoSum(input, 6);
    
    string s = "wlrbbmqbhcdarzowkkyhiddqscdxrjmowfrxsjybldbefsarcbynecdyggxxpklorellnmpapqfwkhopkmco";
    int res = lengthOfLongestSubstring(s);
    cout<<res;
    
    
    return 0;
}

