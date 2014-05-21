//
//  main.cpp
//  LeetCodeIII
//
//  Created by Duan, David on 5/20/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <stack>
#include <assert.h>

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

int findKthElement(int A[], int m , int B[], int n, int k)
{
    if(m > n)
    {
        return findKthElement(B, n, A, m, k);
    }
    if(m == 0)
    {
        return B[k-1];
    }
    if(k == 1)
    {
        return min(A[0], B[0]);
    }
    
    int ia = min(m,k/2);
    int ib = k - ia;
    
    if(A[ia-1] < B[ib-1])
    {
        return findKthElement(A+ia, m-ia, B, n, k-ia);
    }
    else if(B[ib-1] < A[ia-1])
    {
        return findKthElement(A, m, B + ib, n-ib, k-ib);
    }
    else
    {
        return A[ia-1];
    }
}


double findMedianSortedArrays(int A[], int m, int B[], int n)
{
    if((m+n)%2 ==1)
    {
        return findKthElement(A, m, B, n, (m+n)/2+1);
    }
    else
    {
        return (findKthElement(A, m, B, n, (m+n)/2) +
                findKthElement(A, m, B, n, (m+n)/2+1))/2.0;
    }
}

bool isMatchI(const char *s, const char *p)
{
    if( *p == '\0') return *s == '\0';
    if(*(p+1) == '*')
    {
        if((*p == '.' && *s != '\0') || *s == *p)
        {
            return isMatchI(s, p+2) || isMatchI((s+1), p);
        }
        else
        {
            return isMatchI(s, p+2);
        }
    }
    else if((*p == '.' && *s != '\0') || *p == *s)
    {
        return isMatchI(s+1, p+1);
    }
    else
    {
        return false;
    }
}

bool isMatchII(const char *s, const char *p) {
    if (*p == '\0') return *s == '\0';
    // next char is not '*', then must match current character
    if (*(p + 1) != '*') {
        if (*p == *s || (*p == '.' && *s != '\0'))
            return isMatchII(s + 1, p + 1);
        else
            return false;
    } else { // next char is '*'
        while (*p == *s || (*p == '.' && *s != '\0')) {
            if (isMatchII(s, p + 2))
                return true;
            s++;
        }
        return isMatchII(s, p + 2);
    }
}

void swap(vector<int>& num, int i, int j)
{
    int tmp = num[i];
    num[i] = num[j];
    num[j] = tmp;
}

void reverse(vector<int>& num, int i, int j)
{
    while(i<j)
    {
        swap(num, i++,j--);
    }
}

void nextPermutation(vector<int> &num) {
    int i = (int)num.size() -1;
    while(i > 0 && num[i] < num[i-1])
    {
        i --;
    }

    if(i==0)
    {
        reverse(num, 0, (int)num.size()-1);
        return;
    }
    
    int k = (int)num.size()-1;
    while(num[k] <= num[i-1]){
        k--;
    }
    
    swap(num, i-1, k);
    reverse(num, i, (int)num.size()-1);
}


void swap(int num[], int i, int j)
{
    int tmp = num[i];
    num[i] = num[j];
    num[j] = tmp;
}

int firstMissingPositive(int A[], int n) {
    for(int i = 0;i< n;i++)
    {
        while(A[i] != i+1)
        {
            if(A[i] <= 0 || A[i] >= n || A[i] == A[A[i]-1])
            {
                break;
            }
            swap(A, i, A[i]-1);
        }
    }
    
    for(int i = 0;i< n;i++)
    {
        if(A[i] != i+1)
        {
            return i+1;
        }
    }
    return n+1;
}

/*
 '?' Matches any single character.
 '*' Matches any sequence of characters (including the empty sequence).
 
 The matching should cover the entire input string (not partial).
 
 The function prototype should be:
 bool isMatch(const char *s, const char *p)
 
 Some examples:
 isMatch("aa","a") → false
 isMatch("aa","aa") → true
 isMatch("aaa","aa") → false
 isMatch("aa", "*") → true
 isMatch("aa", "a*") → true
 isMatch("ab", "?*") → true
 isMatch("aab", "c*a*b") → false
 */
bool isMatch(const char *s, const char *p) {
    if(*p == '\0' ) return *s == '\0';
    if(*s == '\0')
    {
        while(*p == '*')
        {
            p = p+1;
        }
        if(*p != '\0')
            return false;
        else
            return true;
    }
    
    if(*p == '*' && *s != '\0')
    {
        return isMatch(s+1, p) || isMatch(s, p+1);
    }
    else if(*s == *p || (*p == '?' && *s != '\0'))
    {
        return isMatch(s+1, p+1);
    }
    return false;
}

/**
 search for rotation array II
*/
bool search(int A[], int n, int target) {

    int left = 0;
    int right = n-1;
    while(left<=right)
    {
        int mid = left + (right-left) /2;
        if(A[mid] == target || A[left] == target || A[right] == target)
            return true;
        
        if(A[left] < A[mid])
        {
            if(target > A[left] && target < A[mid])
            {
                right = mid -1;
            }
            else
            {
                left = mid + 1;
            }
        }else if( A[mid] < A[left])
        {
            if(target > A[mid] && target < A[right])
            {
                left = mid + 1;
            }
            else
            {
                right = mid -1;
            }
        }
        else
        {
            left ++;
        }
    }
    return false;
}

int largestRectangleArea(vector<int> &height) {
    stack<int> s;
    height.push_back(0);
    int result = 0;
    for (int i = 0; i < height.size(); )
    {
        if (s.empty() || height[i] > height[s.top()])
        {
            s.push(i++);
        }
        else
        {
            int tmp = s.top();
            s.pop();
            result = max(result,height[tmp] * (s.empty() ? i : i - s.top() - 1));
        }
    }
    
    return result;
}

int main(int argc, const char * argv[])
{
    /*
    int A[] = {3};
    int B[] = {1,2};
    
    double res = findMedianSortedArrays(A, 1, B, 2);
    cout<<res<<endl;
    */
    
    /*
    string s = "acaabbaccbbacaabbbb";
    string p = "a*.*b*.*a*aa*a*";
    bool res = isMatch(s.c_str(), p.c_str());
    cout<<res<<endl;
    */
    
    //vector<int> num;
    //num.push_back(1);
    //num.push_back(2);
    
    //nextPermutation(num);
    
    /*
    string s = "aaaabaaaabbbbaabbbaabbaababbabbaaaababaaabbbbbbaabbbabababbaaabaabaaaaaabbaabbbbaababbababaabbbaababbbba";
    string p = "*****b*aba***babaa*bbaba***a*aaba*b*aa**a*b**ba***a*a*";
    bool res = isMatch(s.c_str(), p.c_str());
    cout<<res<<endl;
    */
    
    vector<int> height;
    height.push_back(2);
    height.push_back(3);
    
    int res = largestRectangleArea(height);
    cout <<res<<endl;
    
    // insert code here...
    //std::cout << "Hello, World!\n";
    return 0;
}

