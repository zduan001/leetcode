//
//  main.cpp
//  LeetCodeIV
//
//  Created by Duan, David on 7/21/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <stack>
#include <assert.h>
#include <numeric>
#include <unordered_set>
#include <queue>

using namespace::std;

struct ListNode {
    int val;
    ListNode *next;
    ListNode(int x) : val(x), next(NULL) {}
};

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode(int x) : val(x), left(NULL), right(NULL) {}
};

struct Point {
    int x;
    int y;
    Point() : x(0), y(0) {}
    Point(int a, int b) : x(a), y(b) {}
};

struct Interval {
    int start;
    int end;
    Interval() : start(0), end(0) {}
    Interval(int s, int e) : start(s), end(e) {}
};

struct TreeLinkNode {
    int val;
    TreeLinkNode *left, *right, *next;
    TreeLinkNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
};

double findKthElement(int k, int A[], int m, int B[], int n)
{
    if(m <=0)
        return B[k-1];
    if(n <=0)
        return A[k-1];
    
    if(k ==1)
        return min(A[0], B[0]);
    
    int a_key = k/2 > m ? INT_MAX : A[k/2-1];
    int b_key = k/2 > n ? INT_MAX : B[k/2-1];
    
    if(a_key < b_key)
        return findKthElement(k - k/2, A + k/2, m - k/2, B, n);
    else
        return findKthElement(k - k/2,  A, m, B+ k/2, n - k/2);
        
}

double findMedianSortedArrays(int A[], int m, int B[], int n) {
    if((n+m)%2==1)
    {
        return findKthElement((n+m)/2+1, A, m, B, n);
    }
    else
    {
        return (double)(findKthElement((n+m)/2, A, m, B, n) + findKthElement((n+m)/2+1, A, m, B, n))/2;
    }
}

bool isMatch(const char *s, const char *p) {
    if(*p== '\0') return *s == '\0';
    
    if(*(p+1) == '*')
    {
        if ((*p == '.' && *s != '\0' ) || *p == *s)
        {
            return isMatch(s+1, p) || isMatch(s, p+2);
        }
        else
        {
            return isMatch(s, p+2);
        }
        
    }
    else if(*p == '.' && *s != '\0')
    {
        return isMatch(s+1, p+1);
    }
    else if(*p == *s)
    {
        return isMatch(s+1, p+1);
    }
    else
    {
        return false;
    }
}


int main(int argc, const char * argv[])
{
/*    int A[] = {1000};
    int B[] = {1001};
    double res = findMedianSortedArrays(A, 1, B, 1);
    cout<<res<<endl;
 */
    cout<<isMatch("aaa", "a*a")<<endl;
    cout<<isMatch("aa","a")<<endl;
    cout<<isMatch("aa","aa")<<endl;
    cout<<isMatch("aaa","aa")<<endl;
    cout<<isMatch("aa", "a*")<<endl;
    cout<<isMatch("aa", ".*")<<endl;
    cout<<isMatch("ab", ".*")<<endl;
    cout<<isMatch("aab", "c*a*b")<<endl;

    return 0;
}

