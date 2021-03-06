//
//  main.cpp
//  projectIII
//
//  Created by Duan, David on 6/19/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <assert.h>
#include <numeric>
#include <queue>
#include <math.h>
#include <stdlib.h>
#include <chrono>
#include <thread>
//#incldue <sstream>

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

struct RandomListNode {
    int label;
    RandomListNode *next, *random;
    RandomListNode(int x) : label(x), next(NULL), random(NULL) {}
};

struct DoubleLinkedListNode{
    int key;
    int val;
    DoubleLinkedListNode *prev;
    DoubleLinkedListNode *next;
    DoubleLinkedListNode(int x, int y) : key(x), val(y), next(NULL), prev(NULL) {}
};

struct GraphNode{
    int val;
    vector<GraphNode*> neighbors;
    bool visited;
    bool inTree;
    GraphNode* parent;
    GraphNode(int x) : val(x), visited(false), inTree(false), parent(NULL) {}
};

struct TrieNode{
    TrieNode* val[256];
    int count[256];
    bool end[256];
    TrieNode()
    {
        fill_n(&count[0], 256, 0);
        fill_n(&val[0], 256, nullptr);
        fill_n(&end[0], 256, false);
    }
};


double findKthElement(int k, int A[], int m, int B[], int n)
{
    if(m<=0)
        return B[k-1];
    if(n <=0)
        return A[k-1];
    
    if(k == 1)
        return min(A[0], B[0]);
    
    int a_key = k/2 - 1 >= m ? INT_MAX: A[k/2-1];
    int b_key = k/2 - 1 >= n ? INT_MAX: B[k/2-1];
    
    if(a_key <b_key)
        return findKthElement(k - k/2, A+k/2, m - k/2, B, n);
    else
        return findKthElement(k - k/2, A, m, B+ k/2, n-k/2);
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
    if(*p == '\0') return *s == '\0';
    
    if(*(p+1) == '*')
    {
        if(*p == *s || (*p == '.' && *s != '\0'))
        {
            return isMatch(s, p+2) || isMatch(s+1, p);
        }
        else
        {
            return isMatch(s, p+2);
        }
    }
    else if((*p == '.' && *s != '\0') || *p == *s)
    {
        return isMatch(s+1, p+1);
    }
    else
    {
        return false;
    }
}

void swap(vector<int>& num, int i , int j)
{
    int tmp = num[i];
    num[i] = num[j];
    num[j] = tmp;
}

void reverse(vector<int>& num , int i, int j)
{
    while(i < j)
        swap(num, i++, j--);
}

void nextPermutation(vector<int> &num) {

    int i = (int)num.size() -1;
    while(i>0)
    {
        if(num[i-1] >= num[i])
            i--;
        else
            break;
    }
    
    if(i ==0)
    {
        //reverse(num.begin(), num.end());
        reverse(num, 0, (int)num.size() -1);
        return;
    }
    i--;
    
    int j = (int)num.size() -1;
    while(j>i)
    {
        if(num[i] >= num[j]) j --;
        else break;
    }
    swap(num, i, j);
    
    reverse(num, i+1, (int)num.size()-1);
    
}

void swap(int A[], int i, int j)
{
    int tmp = A[i];
    A[i] = A[j];
    A[j] = tmp;
}

int firstMissingPositive(int A[], int n) {
    int i = 0;
    while(i< n)
    {
        if((A[i] < 1 || A[i] > n) || A[i] == A[A[i]-1])
            i++;
        else
            swap(A, i, A[i] -1);
    }
    
    for(int i = 0;i< n;i++)
    {
        if(i+1 != A[i])
            return i+1;
    }
    return n+1;
}

bool isMatchII(const char *s, const char *p) {
    if(*p == '\0') return *s == '\0';
    
    if(*p == '*')
    {
        if(*s != '\0')
        {
            return isMatchII(s, p+1) || isMatchII(s+1, p);
        }
        else
        {
            return isMatchII(s, p+1);
        }
    }
    else if((*p == '?' && *s != '\0') || *p == *s)
    {
        return isMatchII(s+1, p+1);
    }
    
    return false;
        
}

bool search(int A[], int n, int target) {
    int left = 0;
    int right = n-1;
    
    while(left <= right)
    {
        int mid = left + (right-left)/2;
        if(target == A[mid] || target == A[left] || target == A[right])
            return true;
        if(A[left] < A[mid])
        {
            if(target > A[left] && target < A[mid])
            {
                right = mid-1;
            }
            else
            {
                left = mid +1;
                
            }
        }
        else if(A[mid] < A[right] )
        {
            if(target > A[mid] && target < A[right])
            {
                left = mid+1;
            }
            else
            {
                right = mid-1;
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
    if (height.size() == 0) return 0;
    height.push_back(0);
    stack<int> s;
    int res = 0;
    for(int i =0;i< height.size();i++)
    {
        if(s.empty() || height[i] >= height[s.top()])
        {
            s.push(i);
        }
        else
        {
            while(!s.empty() && height[i] < height[s.top()])
            {
                int tmp = s.top();
                s.pop();
                if(!s.empty())
                    res = max(res, height[tmp] * (i - s.top()-1));
                else
                    res = max(res, height[tmp] * i);
            }
            s.push(i);
        }
    }
    return res;
}

int maximalRectangle(vector<vector<char> > &matrix) {
    if(matrix.size() == 0 || matrix[0].size() == 0) return 0;
    int res = 0;
    vector<int> q(matrix[0].size(), 0);

    for(int i = 0;i< matrix.size();i++)
    {
        for(int j = 0;j< matrix[i].size();j++)
        {
            if(matrix[i][j] == '1')
                q[j]++;
            else
                q[j] = 0;
        }
        
        res = max(res, largestRectangleArea(q));
    }
    return res;
}

bool isInterleave(string s1, string s2, string s3) {
    if(s1.length() + s2.length() != s3.length()) return false;
    if(s1.length() == 0 && s2.length() == 0  && s3.length() == 0) return true;

    if(s1[0] == s3[0] && s2[0] == s3[0])
    {
        return isInterleave(s1.substr(1, s1.length()-1), s2,s3.substr(1, s3.length()-1)) ||
        isInterleave(s1, s2.substr(1, s2.length()-1), s3.substr(1, s3.length()-1));
    }
    else if(s1[0] == s3[0])
    {
        return isInterleave(s1.substr(1,s1.length()-1), s2, s3.substr(1, s3.length()-1));
    }
    else if(s2[0] == s3[0])
    {
        return isInterleave(s1, s2.substr(1, s2.length()-1), s3.substr(1, s3.length()-1));
    }
    else
        return false;
}

bool isInterleaveDp(string s1, string s2, string s3){
    if(s1.length() + s2.length() != s3.length()) return false;
    vector<vector<bool>> tracker(s1.length()+1, vector<bool>(s2.length() +1, true));
    
    for(int i = 1;i< s1.length()+1;i++)
    {
        tracker[i][0] = s1[i-1] == s3[i-1] && tracker[i-1][0];
    }
    
    for(int j = 1;j< s2.length()+1;j++)
    {
        tracker[0][j] = s2[j-1] == s3[j-1] && tracker[0][j-1];
    }
    
    for(int i = 1;i< s1.length()+1;i++)
    {
        for(int j = 1;j<s2.length();j++)
        {
            tracker[i][j] = (s1[i-1] == s3[i+j-1] && tracker[i-1][j]) ||
            (s2[j-1] == s3[i+j-1] && tracker[i][j-1]);
            
        }
    }
    return tracker[s1.length()][s2.length()];
    
}

string longestPalindrome(string s)
{
    string res = "";
    for (int i = 0;i< s.length();i++)
    {
        int left = i;
        int right = i;
        while(left >=0 && right < s.length() && s[left] == s[right] )
        {
            if((right - left +1)>res.length())
            {
                res = s.substr(left, right-left+1);
            }
            left --;
            right ++;
        }
        left = i;
        right = i+1;
        while(left >=0 && right < s.length() && s[left] == s[right] )
        {
            if((right - left +1)>res.length())
            {
                res = s.substr(left, right-left+1);
            }
            left --;
            right ++;
        }
    }
    return res;
}

ListNode* reverse(ListNode* head)
{
    ListNode* newHead=NULL;
    ListNode* tmp;
    while(head)
    {
        tmp = head;
        head = head->next;
        tmp->next = newHead;
        newHead = tmp;
    }
    return newHead;
}

ListNode *reverseKGroup(ListNode *head, int k) {
    ListNode* dummyHead = new ListNode(0);
    dummyHead->next = head;
    
    ListNode* prevTail = dummyHead;
    ListNode* currentHead = NULL;
    ListNode* currentTail = NULL;
    ListNode* nextHead =NULL;
    
    while(prevTail)
    {
        currentHead = prevTail->next;
        currentTail = currentHead;
        
        for(int i = 0;i<k;i++)
        {
            currentTail = currentTail->next;
            if(!currentTail)
            {
                return dummyHead;
            }
        }
        nextHead = currentTail->next;
        currentTail->next = NULL;
        
        prevTail->next = reverse(currentHead);
        prevTail = currentHead;
        currentHead->next = nextHead;
        currentHead = NULL;
        currentTail = NULL;
        
    }
    return dummyHead->next;
}

char *strStr(char *haystack, char *needle) {
    int m = (int)strlen(needle);
    int n = (int) strlen(haystack);
    if(m>n) return NULL;
        
    int tracker[m];
    tracker[0] = -1;
    int j = tracker[0];
    for(int i = 1;i< strlen(needle);i++)
    {
        while(j != -1 && needle[j+1] == needle[i] )
            j = tracker[j];
        if(needle[j+1] == needle[i])
            j++;
        tracker[i] = j;
    }
    
    j = -1;
    int i = 0;
    while(i<n)
    {
        while(j != -1 && needle[j+1] == haystack[i])
            j = tracker[j];
        if(needle[j+1] == haystack[i])
            j++;
        if(j == m-1)
            return haystack+i-j;
        i++;
    }
    return NULL;
}

int divide(int dividend, int divisor) {
    
    if(divisor == 1) return dividend;
    if(divisor == -1) return -dividend;
    
    
    int sign = 1;
    if (dividend < 0)
    {
        sign *= -1;
        dividend *= -1;
    }
    if (divisor < 0)
    {
        sign *= -1;
        divisor *= -1;
    }
    
    int res = 0;
    int b = divisor;
    
    while(divisor< dividend)
    {
        b= divisor;
        int i = 0;
        while(b<=dividend)
        {
            res += 1<<i;
            dividend -= b;
            b = b<<1;
            i++;
        }
        
    }
    return res;
}

int longestValidParentheses(string s) {
    int res = 0;
    int last = -1;
    stack<char> tracker;
    for(int i = 0;i< s.length();i++)
    {
        if(s[i] == '(')
            tracker.push(i);
        if(s[i] == ')')
        {
            if(tracker.empty())
            {
                last = i;
            }
            else
            {
                tracker.pop();
                if(tracker.empty())
                {
                res = max(res, i - last);
                }
                else
                {
                   res = max(res, i - tracker.top());
                }
            }
        }
    }
    return res;
}

vector<int> searchRange(int A[], int n, int target) {
    vector<int> res;
    if(n<1) return res;
    int left = 0;
    int right = n -1;
    int mid= 0;
    while(left<=right)
    {
        mid = left + (right - left) / 2;
        if(A[mid] < target)
        {
            left = mid + 1;
        }
        else if(A[mid] > target)
        {
            right = mid - 1;
        }
        else
        {
            break;
        }
    }
    
    
    left=0;
    right = mid;
    int tmp;
    while(left <= right)
    {
        tmp = left + (right-left) /2;
        if(tmp == 0 || (A[tmp] == target && A[tmp-1] < target))
        {
            if(A[tmp] == target)
                break;
            else
                left = tmp+1;
        }
        else if(A[tmp] == target)
        {
            right = tmp-1;
        }
        else
        {
            left = tmp+1;
        }
    }
    if(left>right)
    {
        res.push_back(-1);
        res.push_back(-1);
        return res;
    }
    
    res.push_back(tmp);
    left = mid;
    right = n-1;
    
    while(left <= right)
    {
        tmp = left + (right-left) /2;
        if(tmp == n-1 || (A[tmp] == target && A[tmp+1] > target))
        {
            if(A[tmp] == target)
                break;
            else
                right = tmp -1;
        }
        else if(A[tmp] == target)
        {
            left = tmp+1;
        }
        else
        {
            right = tmp-1;
        }
    }
    res.push_back(tmp);
    
    return res;
}

////////////////////////////////////////////////
//
//
// Combination II
//
//
////////////////////////////////////////////////
void worker(vector<int>& num, int gap, vector<int>& tmp, int level, vector<vector<int>>& res)
{
    if(0 == gap)
    {
        res.push_back(tmp);
        return;
    }
    else
    {
        int prev = -1;
        for(int i = level;i< num.size();i++)
        {
            if(num[i] == prev)
            {
                continue;
            }
            if(gap < num[i])
            {
                return;
            }
            prev = num[i];
            tmp.push_back(num[i]);
            worker(num, gap - num[i], tmp, i+1, res);
            tmp.pop_back();
        }
    }
}
vector<vector<int> > combinationSum2II(vector<int> &num, int target) {
    
    sort(num.begin(), num.end());
    vector<vector<int>> res;
    vector<int> tmp;
    worker(num, target, tmp, 0, res);
    return res;
}

////////////////////////////////////////////////
//
//
// Permutation II
//
//
////////////////////////////////////////////////
void worker(vector<int>& num, int level, vector<vector<int>>& res){
    if(level == num.size())
    {
        res.push_back(num);
        return;
    }
    else
    {
        unordered_set<int> set;
        for(int i = level;i< num.size();i++)
        {
            if(set.find(num[i]) == set.end())
            {
                swap(num, i, level);
                worker(num, level+1, res);
                swap(num, i, level);
                set.insert(num[i]);
            }
        }
    }
}
vector<vector<int>> permuteUnique(vector<int> &num) {
    vector<vector<int>> res;
    sort(num.begin(), num.end());
    worker(num, 0, res);
    return res;
}

////////////////////////////////////////////////
//
//
// Subset II
//
//
////////////////////////////////////////////////
void worker(vector<int> &S, int level, vector<int>& tmp, vector<vector<int>>& res)
{
    res.push_back(tmp);
    for(int i = level;i< S.size();i++)
    {
        if(i == level || S[i] != S[i-1])
        {
            tmp.push_back(S[i]);
            worker(S, i+1, tmp, res);
            tmp.pop_back();
        }
    }
}
vector<vector<int> > subsetsWithDup(vector<int> &S) {
    vector<vector<int>> res;
    if(S.size() == 0) return res;
    
    sort(S.begin(), S.end());
    vector<int> tmp;
    worker(S, 0, tmp, res);
    return res;
}


int trap(int A[], int n) {
    if(n <=2) return 0;
    int res = 0;
    int left = 0;
    int right = n-1;
    int leftHight= 0;
    int rightHight = 0;
    
    while(left <= right)
    {
        if(A[left] <=  A[right])
        {
            if(leftHight < A[left])
            {
                leftHight = A[left];
            }
            else
            {
                res += leftHight - A[left];
            }
            left++;
        }
        else
        {
            if(rightHight <=A[right])
            {
                rightHight = A[right];
            }
            else
            {
                res += rightHight -A[right];
            }
            right--;
        }
    }
    return res;
}

string multiply(string num1, string num2) {
    if(num1 == "0" || num2 == "0") return "0";
    int m = (int)num1.length();
    int n = (int)num2.length();
    vector<int> res(m+n, 0);
    vector<char> charnum1 (num1.begin(), num1.end());
    vector<char> charnum2 (num2.begin(), num2.end());
    
    reverse(charnum1.begin(), charnum1.end());
    reverse(charnum2.begin(), charnum2.end());
    
    int carry = 0;
    int tmp = 0;
    
    for(int i = 0;i< charnum1.size();i++)
    {
        for(int j = 0;j<charnum2.size();j++)
        {
            tmp = (charnum1[i]-'0')*(charnum2[j]-'0') + carry;
            tmp += res[i+j];
            carry = tmp/10;
            res[i+j] = tmp%10;
        }
        res[i+charnum2.size()]+= carry;
        carry = 0;
    }
    if (res[n+m] == '0')
    {
        res.pop_back();
    }
    
    string s = "";
    for(int i = 0;i< res.size();i++)
    {
        s += res[i]+ '0';
    }
    
    reverse(s.begin(), s.end());
    return s;
}

int jump(int A[], int n) {
    if(n == 1) return 0;
    unordered_set<int> visited;
    queue<int> tracker;
    tracker.push(0);
    tracker.push(-1);
    int layer = 0;
    while(!tracker.empty())
    {
        int tmp = tracker.front();
        if(tmp == -1)
        {
            layer++;
            tracker.pop();
            if(tracker.empty())
            {
                break;
            }
            else
            {
                tracker.push(-1);
            }
            continue;
        }
        
        tracker.pop();
        for(int i = 1;i<= A[tmp];i++)
        {
            if(tmp+i>=n-1)
            {
                return layer+1;
            }
            
            if(visited.find(tmp+i) == visited.end())
            {
                visited.insert(tmp+i);
                tracker.push(tmp+i);
            }
        }
        
    }
    // should not reach.
    return -1;
}

int minDistance(string word1, string word2) {
    int n = (int)word1.length();
    int m = (int)word2.length();
    
    vector<vector<int>> tracker(n+1, vector<int>(m+1));
    
    for(int i = 1;i< n+1;i++)
    {
        tracker[i][0] = i;
    }
    
    for(int i = 1;i<m+1;i++)
    {
        tracker[0][i] = i;
    }
    
    
    for(int i = 1;i< word1.length();i++)
    {
        for(int j = 1;j< word2.length();j++)
        {
            int tmp = word1[i-1] == word2[j-1] ? tracker[i-1][j-1] : tracker[i-1][j-1] + 1;
            tmp = min (tmp, tracker[i-1][j]+1);
            tmp = min (tmp, tracker[i][j-1]+1);
            tracker[i][j] = tmp;
        }
    }
    return tracker[n][m];
}

void rotate(vector<vector<int> > &matrix) {
    int n = (int)matrix.size();
    if(n <2) return;
    
    for(int i = 0;i< n/2;i++)
    {
        for(int j = i; j< n-i-1;j++)
        {
            int tmp = matrix[i][j];
            matrix[i][j] = matrix[n-j-1][i];
            matrix[n-j-1][i] = matrix[n-i-1][n-j-1];
            matrix[n-i-1][n-j-1] = matrix[j][n-i-1];
            matrix[j][n-i-1] = tmp;
        }
    }
    return;
}

vector<int> spiralOrder(vector<vector<int> > &matrix)
{
    int n = (int)matrix.size();
    int m = (int)matrix[0].size();
    vector<int> res;
    if(n == 0 || m == 0) return res;
    
    int top = 0;
    int bottom = n-1;
    int left = 0;
    int right = m-1;
    int i = 0;int j = 0;
    int dir = 0;
    
    while(top<=bottom && left<=right)
    {
        res.push_back(matrix[i][j]);
        if(dir == 0)
        {
            j++;
            if( j >= right)
            {
                j = right;
                dir ++;
                top ++;
            }
        }
        else if(dir == 1)
        {
            i++;
            if(i >= bottom)
            {
                i = bottom;
                dir++;
                right--;
            }
        }
        else if(dir == 2)
        {
            j --;
            if(j <= left)
            {
                j = left;
                dir ++;
                bottom--;
            }
        }
        else
        {
            i--;
            if(i <= top)
            {
                i = top;
                dir = 0;
                left++;
            }
        }
    }
    return res;
}

bool sortIntervalFunctions(Interval i1, Interval i2){
    return (i1.start < i2.start);
}

vector<Interval> merge(vector<Interval> &intervals) {
    sort(intervals.begin(), intervals.end(), sortIntervalFunctions);
    vector<Interval> res;
    Interval tmp = intervals[0];
    for(int i = 1;i< intervals.size();i++)
    {
        if(tmp.end >= intervals[i].start)
        {
            tmp.start = min(tmp.start, intervals[i].start);
            tmp.end = max(tmp.end, intervals[i].end);
        }
        else
        {
            res.push_back(tmp);
            tmp = intervals[i];
        }
    }
    res.push_back(tmp);
    return res;
}



vector<Interval> insert(vector<Interval> &intervals, Interval newInterval) {
    vector<Interval> res;
    
    int i = 0;
    
    while(intervals[i].end < newInterval.start)
    {
        res.push_back(intervals[i]);
        i++;
    }
    
    while(newInterval.end >= intervals[i].start)
    {
        newInterval.start = min(newInterval.start, intervals[i].start);
        newInterval.end = max(newInterval.end, intervals[i].end);
        i++;
        
    }
    res.push_back(newInterval);
    
    while(i< intervals.size())
    {
        res.push_back(intervals[i]);
        i++;
    }
    return res;
}

vector<string> fullJustify(vector<string> &words, int L) {
    vector<string> tracker;
    vector<string> res;
    if(words.size() == 0)
    {
        return res;
    }
    
    int length = 0;
    int wordLength = 0;
    tracker.push_back(words[0]);
    length = (int)words[0].length();
    wordLength = (int)words[0].length();
    for(int i = 1;i< words.size();i++)
    {
        if(length + words[i].size() +1 > L )
        {
            if(tracker.size() == 1)
            {
                string tmp = tracker[0];
                while(tmp.length() < L)
                {
                    tmp+= " ";
                }
                res.push_back(tmp);
            }
            else
            {
                int space = (L - wordLength) / (tracker.size() -1);
                int over = (L - wordLength) % (tracker.size() -1);
                
                string tmp = "";
                for(int j = 0;j< tracker.size();j++)
                {
                    if(tmp != "")
                    {
                        for(int k = 0;k<space;k++)
                        {
                            tmp += " ";
                        }
                        if(over >0)
                        {
                            tmp += " ";
                            over--;
                        }
                    }
                    tmp += tracker[j];
                }
                res.push_back(tmp);
            }
            tracker.clear();
            tracker.push_back(words[i]);
            length = (int)words[i].length();
            wordLength = (int)words[i].length();
        }
        else
        {
            tracker.push_back(words[i]);
            length += (words[i].size() + 1);
            wordLength += words[i].length();
        }
    }
    
    if(tracker.size() > 0)
    {
        string tmp = "";
        for(int j = 0;j< tracker.size();j++)
        {
            if(tmp != "")
            {
                tmp += " ";
            }
            tmp += tracker[j];
        }
        res.push_back(tmp);
    }
    
    return res;
}

void connect(TreeLinkNode *root) {
    if(!root) return;
    TreeLinkNode* prev = NULL;
    TreeLinkNode* nextStart = NULL;
    while(root)
    {
        prev = NULL;
        nextStart = NULL;
        while(root)
        {
        if(!nextStart)
        {
            nextStart = root->left ? root->left : root->right;
        }
        if(root->left)
        {
            if(prev) prev->next = root->left;
            prev = root->left;
                
        }
        if(root ->right)
        {
            if(prev) prev->next = root->right;
            prev = root->right;
        }
        root = root->next;
        }
        root = nextStart;
    }
    
}

void sortColors(int A[], int n) {
    int left = 0;
    int right = n-1;
    int runner = 0;
    while(runner <=right)
    {
        if(A[runner] == 0)
            swap(A, left++, runner++);
        else if(A[runner] == 2)
            swap(A, runner, right--);
        else
            runner ++;
    }
}

vector<vector<int> > zigzagLevelOrder(TreeNode *root) {
    vector<vector<int>> res;
    if(!root) return res;
    stack<TreeNode*> s1;
    stack<TreeNode*> s2;
    bool inOrder = true;
    s1.push(root);
    vector<int>* tracker = new vector<int>();
    while(!s1.empty())
    {
        TreeNode* tmp = s1.top();
        tracker->push_back(tmp->val);
        if(inOrder)
        {
            if(tmp->left)
                s2.push(tmp->left);
            if(tmp->right)
                s2.push(tmp->right);
        }
        else
        {
            if(tmp->right)
                s2.push(tmp->right);
            if(tmp->left)
                s2.push(tmp->left);
        }
        s1.pop();
        if(s1.empty())
        {
            swap(s1, s2);
            res.push_back(*tracker);
            tracker = new vector<int>();
            inOrder = !inOrder;
        }
    }
    return res;
}

TreeNode *sortedListToBST(ListNode *head) {
    if(!head) return NULL;
    ListNode* fast = head->next;
    ListNode* slow = head;
    while(fast && fast->next && fast->next->next)
    {
        slow = slow->next;
        fast = fast->next->next;
    }
    
    ListNode* mid = slow->next;
    if(!mid)
    {
        return new TreeNode(slow->val);
    }
    else
    {
        ListNode* newHead = mid->next;
        slow->next = NULL;
        TreeNode* root = new TreeNode(mid->val);
        root->left = sortedListToBST(head);
        root->right = sortedListToBST(newHead);
        return root;
    }
}

int numDistinct(string S, string T) {
    int n = (int)S.length();
    int m = (int)T.length();
    if(m>n) return 0;
    
    int tracker[m+1][n+1];
    fill_n(&tracker[0][0], m*n, 0);
    for(int i = 0;i<=n;i++)
    {
        tracker[0][i] = 1;
    }
    for(int i = 1;i<=m; i++)
    {
        tracker[i][0] = 0;
    }
    
    for(int i = 1;i<=n;i++)
    {
        for(int j = 1; j<=m; j++)
        {
            tracker[i][j] = tracker[i][j-1] + S[i] == T[j]? tracker[i-1][j-1] : 0;
        }
    }
    
    return tracker[m+1][n+1];
}


int MAX_SUM=INT_MIN;

int worker(TreeNode *root) {
    if(!root) return 0;
    int left = worker(root->left);
    int right = worker(root->right);
    
    int sum = root->val;
    if(left >0) sum += left;
    if(right>0) sum+= right;
    MAX_SUM = max(MAX_SUM, sum);
    
    return max(left, right) > 0? root->val + max(left, right) : root->val;
}

int maxPathSum(TreeNode *root) {
    int val = worker(root);
    return max(MAX_SUM, val);
}

int longestConsecutive(vector<int> &num) {
    unordered_map<int, bool> tracker;
    for(int i = 0;i< num.size();i++)
    {
        if(tracker.find(num[i]) == tracker.end())
        {
            tracker.insert(make_pair(num[i], false));
        }
    }
    
    int maxLength =0;
    for(int i = 0;i< num.size();i++)
    {
        if(!(tracker[num[i]]))
        {
            tracker.find(num[i])->second = true;
            int tmp = 1;
            
            int runner = num[i]+1;
            while(tracker.find(runner) != tracker.end() && !tracker[runner])
            {
                tmp ++;
                tracker[runner]= true;
                runner++;
            }
            runner = num[i] - 1;
            while(tracker.find(runner) != tracker.end() && !tracker[runner])
            {
                tmp ++;
                tracker[runner] = true;
                runner --;
            }
            maxLength = max(maxLength, tmp);
            
        }
    }
    return maxLength;
}

int solve(vector <bool>& m, string& ans, int k) {
    if (k == m.size()) {
        cout << ans << endl;
        return 0;
    }
    
    for (int i=0; i<m.size(); ++i) {
        if (!m[i]) {
            ans = ans+char(i+'1');
            m[i] = true;
            solve(m, ans, k+1);
            m[i] = false;
            ans.erase(ans.size()-1);
        }
    }
    return 0;
}

int permutation(int n) {
    assert(n>0);
    vector <bool> m(n, false);
    string ans;
    solve(m, ans, 0);
    return 0;
}


void fill(vector<vector<char>>& board, int i, int j, char origin, char target)
{
    if(i <0 || i >= board.size() || j < 0 || j >= board[0].size())
    {
        return;
    }
    
    if(board[i][j] == origin)
    {
        board[i][j] = target;
        fill(board, i-1, j, origin, target);
        fill(board, i+1, j, origin, target);
        fill(board, i, j-1, origin, target);
        fill(board, i, j+1, origin, target);
    }
}

void solve(vector<vector<char>> &board) {
    for(int i = 0;i<board.size();i++)
    {
        if(board[i][0] == 'O')
        {
            fill(board, i, 0, 'O', 'A');
        }
        if(board[i][(int)board[0].size()-1])
        {
            fill(board, i, (int)board[0].size()-1, 'O', 'A');
        }
    }
    
    for(int j = 0;j<board[0].size();j++)
    {
        if(board[0][j] == 'O')
        {
            fill(board, 0, j, 'O', 'A');
        }
        if(board[(int)board.size()-1][j] == 'O')
        {
            fill(board, (int)board.size()-1, j, 'O', 'A');
        }
    }
    
    for(int i = 0;i< board.size();i++)
    {
        for(int j = 0;j< (int)board[0].size(); j++)
        {
            if(board[i][j] == 'O')
            {
                fill(board, i, j, 'O', 'X');
            }
        }
    }
    
    for(int i = 0;i< board.size();i++)
    {
        for(int j = 0;j< (int)board[0].size(); j++)
        {
            if(board[i][j] == 'A')
            {
                board[i][j] = 'O';
            }
        }
    }
    
}

void interleave(int a[], int n)
{
    for(int i = (n-1)/2 ;i>0;i--)
    {
        for(int j = i;j< (n-1-i);j=j+2)
        {
            swap(a, j, j+1);
        }
    }
}

int worker(int k)
{
    int count = 0;
    if(k == 0)
        return 0;
    if(k == 1)
        return 1;
    while(k >=2)
    {
        count += pow(2,k--);
    }
    return count;
}

int countBits(int n)
{
    int count = 0;
    while(n>1)
    {
        int i = (int)log2(n+1);
        
        count += worker(i);
        n -= (pow(2,i)-1);
    }
    return count;
}

int SET_BIT=0;

void totalBits(int n)
{
    if(n==1)
    {
        SET_BIT++;
        return;
    }
    if(n/2)
        SET_BIT++;
    
    int Pow=log2(n);
    n-=pow(2,Pow);
    
    if(n)
        totalBits(n);

    return;
}

/*
  Brian Kernighan’s Algorithm:
 */
int countSetBits(int n)
{
    unsigned int count = 0;
    while (n)
    {
        n &= (n-1) ;
        count++;
    }
    return count;
}



//2.1
int removeDuplicates(int A[], int n) {
    if(n<2) return n;
    int index = 1;
    int i = index;
    while(i<n)
    {
        if(A[index-1] ==A[i])
        {
            i++;
        }
        else
        {
            A[index++] = A[i++];
        }
    }
    return index;
}

//2.2
int removeDuplicatesII(int A[], int n) {
    if(n<2) return n;
    int index = 2;
    int i = index;
    while(i<n)
    {
        if(A[index-2] != A[i])
        {
            A[index++] = A[i++];
        }
        else
        {
            i++;
        }
    }
    return index+1;
}

//2.1.3 and 2.1.4
int searchInRotateArray(int A[], int n, int target)
{
    int left = 0;
    int right = n-1;
    while(left<=right)
    {
        int mid = left + (right - left) /2;
        if(A[mid] == target) return mid;
        
        if(A[left] < A[mid])
        {
            if(target >= A[left] && target <= A[mid])
            {
                right = mid-1;
            }
            else
            {
                left = mid+1;
            }
        }
        else if(A[mid] < A[right])
        {
            if(target >= A[mid] && target <= A[right])
            {
                left = mid + 1;
            }
            else
            {
                right = mid - 1;
            }
        }
        else
        {
            if(A[left] != target)
                left++;
            else
                return left;
        }
    }
    return -1;
}

//2.1.5
int findKthElement(int A[], int m, int B[], int n, int k)
{
    if(m > n)
    {
        return findKthElement(B, n, A, m, k);
    }
    
    if(m == 0) return B[k-1];
    if(k ==1) return min(A[0], B[0]);
    
    int ia = min(k/2, m);
    int ib = k -ia;
    
    if(A[ia-1] < B[ib-1])
    {
        return findKthElement(A+ia, m-ia, B, n, k-ia);
    }
    else if(B[ib-1] < A[ia-1])
    {
        return findKthElement(A, m, B+ib, n-ib, k-ib);
    }
    else
    {
        return A[ia-1];
    }
}


double findMedianSortedArraysII(int A[], int m, int B[], int n) {
    if((n+m) & 1)
    {
        int left = findKthElement(A, m, B, n, (n+m)/2);
        int right= findKthElement(A,m,B,n,(n+m)/2 +1);
        return double(left+right)/2.0;
    }
    else
    {
        return findKthElement(A, m, B, n, (n+m)/2+1);
    }
    
}

//2.1.6
int longestConsecutiveII(vector<int> &num)
{
    unordered_map<int, bool> tracker;
    for(int i = 0;i< num.size();i++)
    {
        tracker[num[i]] = false;
    }
    
    int res = 0;
    int count = 0;
    for(int i = 0;i< num.size();i++)
    {
        if(tracker[num[i]] == false)
        {
            count = 0;
            for(int j = num[i] + 1;tracker.find(j) !=tracker.end();j++)
            {
                count ++;
                tracker[j] = true;
            }
            
            for(int j = num[i] -1;tracker.find(j) !=tracker.end();j--)
            {
                count ++;
                tracker[j] = true;
            }
            count ++;
        }
        res = max(res, count);
    }
    return res;
}

//2.1.7
vector<int> twoSum(vector<int> &numbers, int target) {
    vector<int> res;
    unordered_map<int,vector<int>> tracker;
    
    for(int i = 0;i< numbers.size();i++)
    {
        if(tracker.find(numbers[i]) == tracker.end())
        {
            vector<int> tmp = {i};
            tracker[numbers[i]] = tmp;
        }
        else
        {
            tracker[numbers[i]].push_back(i);
        }
    }
    
    for(int i = 0;i< numbers.size();i++)
    {
        if(tracker.find(target - numbers[i]) != tracker.end())
        {
            if(target-numbers[i] == numbers[i])
            {
                if(tracker[numbers[i]].size() > 1)
                {
                    res.push_back(tracker[numbers[i]][0]);
                    res.push_back(tracker[numbers[i]][1]);
                    break;
                }
            }
            else
            {
                res.push_back(tracker[numbers[i]][0]);
                res.push_back(tracker[target-numbers[i]][0]);
                break;
            }
        }
    }
    if(res.size() ==2)
    {
        res[0] ++;
        res[1] ++;
    }
    sort(res.begin(), res.end());
    return res;
}

// this method is shorter, while is go through numbers,
// it add it to the map. save logic on check duplicate and
// same item ....
vector<int> twoSumII(vector<int> &numbers, int target) {
    vector<int> res;
    unordered_map<int, int> tracker;
    for(int i = 0;i<numbers.size();i++)
    {
        int tmp = target - numbers[i];
        if(tracker.find(tmp) != tracker.end())
        {
            res.push_back(tracker[tmp]+1);
            res.push_back(i+1);
            break;
        }
        tracker[numbers[i]] = i;
    }
    return res;
}

//2.1.8
vector<vector<int> > threeSum(vector<int> &num) {
    vector<vector<int>> res;
    if(num.size() <3) return res;
    sort(num.begin(), num.end());
    for(int i = 0;i< num.size()-2;i++)
    {
        if(i == 0 || num[i] != num[i-1])
        {
            int left = i+1;
            int right = (int)num.size()-1;
            while(left <right)
            {
                int tmp = num[i] + num[left] + num[right] ;
                if(tmp == 0)
                {
                    res.push_back({num[i], num[left], num[right]});
                    left ++;
                    right --;
                }
                else if(tmp <0)
                {
                    left++;
                }
                else
                {
                    right --;
                }
            }
        }
    }
    sort(res.begin(), res.end());
    res.erase(unique(res.begin(), res.end()), res.end());
    return res;
}

//2.1.9
int threeSumClosest(vector<int> &num, int target) {
    int distance = INT_MAX;
    int res;
    if(num.size() <3) return res;
    sort(num.begin(), num.end());
    for(auto a = num.begin();a<prev(num.end(),2);a++)
    {
        auto b = next(a);
        auto c = prev(num.end());
        while(b<c)
        {
            int tmp = abs(target - (*a + *b + *c));
            
            if(tmp < distance)
            {
                distance = tmp;
                res = *a + *b + *c;
            }
            if(distance == 0) break;
            if(target > (*a + *b + *c))
            {
                b++;
            }
            else
            {
                c--;
            }
        }
    }
    return res;
}

//2.1.10
vector<vector<int> > fourSum(vector<int> &num, int target) {
    vector<vector<int>> res;
    if(num.size() < 4) return res;
    unordered_map<int, vector<pair<int, int>>> tracker;
    for(int i =0;i< num.size();i++)
    {
        for(int j = i+1;j<num.size();j++)
        {
            int tmp = num[i] + num[j];
            if(tracker.find(tmp) == tracker.end())
            {
                tracker[tmp] = {make_pair(i, j)};
            }
            else
            {
                tracker[tmp].push_back(make_pair(i, j));
            }
        }
    }
    
    for(int i = 0;i<num.size();i++)
    {
        for(int j = i+1;j<num.size();j++)
        {
            int key = num[i] + num[j];
            if(tracker.find(target -key) != tracker.end())
            {
                for(auto a = tracker[target-key].begin();a <tracker[target-key].end();a++)
                {
                    if(i < a->second)
                    {
                        continue;
                    }
                    vector<int> tmp = {num[a->first], num[a->second], num[i], num[j]};
                    sort(tmp.begin(), tmp.end());
                    res.push_back(tmp);
                }
            }
        }
    }
    sort(res.begin(), res.end());
    res.erase(unique(res.begin(), res.end()), res.end());
    return res;
}

//2.1.11
int removeElement(int A[], int n, int elem) {
    int index = 0;
    for(int i = 0;i<n;i++)
    {
        if(A[i] != elem)
        {
            A[index++] = A[i];
        }
    }
    return index+1;
}

//2.1.12
void nextPermutationII(vector<int> &num) {
    int pivot = -1;
    for(int i = (int)num.size()-2;i>=0;i--)
    {
        if(num[i] < num[i+1])
        {
            pivot = i;
            break;
        }
    }
    
    if(pivot ==- 1)
    {
        reverse(num.begin(), num.end());
        return;
    }
    
    int index = -1;
    for(int i = (int)num.size() -1;i>pivot;i--)
    {
        if(num[i] > num[pivot])
        {
            index = i;
            break;
        }
    }
    
    swap(num, pivot, index);
    reverse(num, pivot+1, (int)num.size()-1);
}

//2.1.15
int trapII(int A[], int n) {
    int res = 0;
    if(n<3) return res;
    int left = 0;
    int right = n - 1;
    int leftMax = A[0];
    int rightMax = A[n-1];
    while(left <=right)
    {
        if(leftMax <= rightMax)
        {
            if(A[left] <= leftMax)
            {
                res += leftMax - A[left++];
            }
            else
            {
                leftMax = A[left++];
            }
        }
        else
        {
            if(A[right] <=rightMax)
            {
                res += rightMax - A[right--];
            }
            else
            {
                rightMax = A[right--];
            }
        }
    }
    return res;
}



//2.1.16

//2.1.17
vector<int> plusOne(vector<int> &digits) {
    int carry = 1;
    for(int i = (int)digits.size()-1;i>=0;i--)
    {
        digits[i] += carry;
        carry = digits[i]/10;
        digits[i] = digits[i]%10;
    }
    if(carry >0)
    {
        digits.insert(digits.begin(), carry);
    }
    return digits;
}

//2.1.18
int climbStairs(int n) {
    if(n ==0) return 0;
    if(n ==1) return 1;
    if (n ==2) return 2;
    int prev = 1;
    int cur = 2;
    for( int i = 3;i<=n;i++)
    {
        int tmp = prev+cur;
        prev = cur;
        cur = tmp;
    }
    return cur;
    
}

//2.1.20
void setZeroes(vector<vector<int> > &matrix) {
    int n = (int)matrix.size();
    int m = (int)matrix[0].size();
    
    bool a[n];
    bool b[m];
    for(int i = 0;i<n;i++) a[i] = false;
    for(int j = 0;j<m;j++) b[j] = false;
    for(int i = 0;i< n;i++)
    {
        for(int j = 0;j<m;j++)
        {
            if(matrix[i][j] == 0)
            {
                a[i] = true;
                b[j] = true;
            }
        }
    }
    
    for(int i = 0;i< n;i++)
    {
        if(a[i])
        {
            for(int j = 0;j<m;j++)
            {
                matrix[i][j] =0;
            }
        }
    }
    for(int j = 0;j<m;j++)
    {
        if(b[j])
        {
            for(int i = 0;i<n;i++ )
            {
                matrix[i][j] = 0;
            }
        }
    }
}

//2.1.21
int canCompleteCircuit(vector<int> &gas, vector<int> &cost) {
    int sum = 0;
    int total = 0;
    int start = 0;
    for(int i = 0;i<gas.size();i++)
    {
        sum += (gas[i] - cost[i]);
        total += (gas[i] - cost[i]);
        if(sum <0)
        {
            sum = 0;
            start = i+1;
        }
    }
    if(total>=0) return start>=gas.size()? start-1: start;
    else return -1;
}

//2.1.22
int candy(vector<int> &ratings) {
    int n = (int)ratings.size();
    vector<int> res(n);
    int count = 1;
    for(int i = 1;i< n;i++)
    {
        if(ratings[i] > ratings[i-1])
        {
            res[i] = max(count++, res[i]);
        }
        else
        {
            count = 1;
        }
    }
    count = 1;
    for(int i = n-2;i>=0;i++)
    {
        if(ratings[i] > ratings[i+1])
        {
            res[i] = max(count++, res[i]);
        }
        else
        {
            count = 1;
        }
    }
    return accumulate(&res[0], &res[0]+n, n);
}

//2.1.23
int singleNumber(int A[], int n) {
    int res = 0;
    for(int i = 0;i< n;i++)
    {
        res ^= A[i];
    }
    return res;
}

//2.1.24
int singleNumberII(int A[], int n) {
    int tracker[32];
    fill_n(&tracker[0], 32, 0);
    for(int i = 0;i< n;i++)
    {
        for(int j = 0;j<32;j++)
        {
            tracker[j] += (A[i]>>j) & 1;
        }
    }
    for(int j = 0;j< 32;j++)
    {
        tracker[j] %= 3;
    }
    int res = 0;
    for(int j = 0;j<32;j++)
    {
        if(tracker[j])
        {
            res += 1<<j;
        }
    }
    return res;
}

//2.2.1
ListNode *addTwoNumbers(ListNode *l1, ListNode *l2) {
    if(!l1) return l2;
    if(!l2) return l1;
    
    ListNode* dummy = new ListNode(1);
    ListNode* tmp = dummy;
    int carry = 0;
    while(l1 && l2)
    {
        carry = carry + l1->val + l2->val;
        tmp -> next = new ListNode(carry %10);
        tmp = tmp->next;
        l1 = l1->next;
        l2 = l2->next;
        carry /=10;
    }
    ListNode* remain;
    if(l1) remain = l1;
    else remain = l2;
    while(remain)
    {
        carry += remain->val;
        tmp->next = new ListNode(carry%10);
        tmp = tmp->next;
        remain = remain->next;
        carry /=10;
    }
    if(carry)
    {
        tmp->next = new ListNode(carry);
    }
    return dummy->next;
}

//2.2.2
ListNode *reverseII(ListNode * head)
{
    ListNode *newHead = NULL;
    ListNode *tmp = head;
    while(tmp)
    {
        head = head->next;
        tmp->next = newHead;
        newHead = tmp;
        tmp = head;
    }
    return newHead;
}

ListNode *reverseBetween(ListNode *head, int m, int n) {
    ListNode* dummyHead = new ListNode(-1);
    dummyHead->next = head;
    ListNode* firstEnd = dummyHead;
    ListNode* secondEnd = dummyHead;
    for(int i = 0;i<m-1;i++)
    {
        firstEnd = firstEnd->next;
    }
    for(int j = 0;j<n;j++)
    {
        secondEnd = secondEnd->next;
    }
    
    ListNode *middleHead = firstEnd->next;
    ListNode *lastHead = secondEnd->next;
    secondEnd->next = NULL;
    
    //reverse(middleHead);
    firstEnd->next = reverse(middleHead);
    middleHead->next = lastHead;
    return dummyHead->next;
}

//2.2.3
ListNode *partition(ListNode *head, int x) {
    ListNode* LargHead = new ListNode(-1);
    ListNode* SmallHead = new ListNode(-1);
    ListNode* largeTmp = LargHead;
    ListNode* smallTmp = SmallHead;
    while(head)
    {
        ListNode* tmp = head;
        head = head->next;
        if(tmp->val >= x)
        {
            largeTmp->next = tmp;
            largeTmp = tmp;
            largeTmp->next = NULL;
        }
        else
        {
            smallTmp->next = tmp;
            smallTmp = tmp;
            smallTmp ->next = NULL;
        }
    }
    smallTmp->next = LargHead->next;
    return SmallHead->next;
}

//2.2.4
ListNode *deleteDuplicates(ListNode *head) {

    ListNode* tmp = head;
    while(tmp && tmp->next)
    {
        if(tmp->val == tmp->next->val)
        {
            tmp->next = tmp->next->next;
        }
        else
        {
            tmp = tmp->next;
        }
    }
    return head;
}

//2.2.5
ListNode *deleteDuplicatesII(ListNode *head) {
    ListNode* dummyHead = new ListNode(-1);
    dummyHead->next = head;
    
    ListNode* tmp = dummyHead;
    while(tmp)
    {
        ListNode* runner = tmp->next;
        if(runner)
        {
            int val = runner->val;
            int count = 0;
            while(runner && runner->val == val)
            {
                count ++;
                runner = runner->next;
            }
            if(count == 1 )
            {
                runner = tmp->next;
            }
        }
        
        if(tmp->next == runner)
        {
            tmp = tmp->next;
        }
        else
        {
            tmp->next = runner;
        }
        
    }
    return dummyHead->next;
    
}

//2.2.6
ListNode *rotateRight(ListNode *head, int k) {
    if(k == 0 || !head) return head;
    ListNode* front = head;
    ListNode* back = head;
    int count = -1;
    for(int i =0;i< k;i++)
    {
        if(!front)
        {
            count = i;
            break;
        }
        
        front = front->next;
    }
    if(!front)
    {
        k =k%count;
        if( k == 0) return head;
        front = head;
        for(int i= 0;i<k;i++)
        {
            front = front->next;
        }
    }
    while(front->next)
    {
        back = back->next;
        front = front->next;
    }
    
    ListNode* newHead = back->next;
    back->next = NULL;
    front->next = head;
    return newHead;
}

//2.2.6
ListNode *rotateRightII(ListNode *head, int k) {
    if(k == 0 || !head) return head;
    int n = 1;
    ListNode* p = head;
    while(p->next)
    {
        n++;
        p = p->next;
    }
    k = k%n;
    p->next = head;
    
    ListNode* front = head;
    ListNode* back = head;
    for(int i = 0;i< k;i++)
    {
        front = front->next;
    }
    while(front->next != head)
    {
        front = front->next;
        back = back->next;
    }
    ListNode* newHead = back->next;
    back->next = NULL;
    return newHead;
}

//2.2.7
ListNode *removeNthFromEnd(ListNode *head, int n) {
    //assume n < the length of head.
    ListNode* dummyHead = new ListNode(-1);
    dummyHead->next = head;
    ListNode* front = dummyHead;
    ListNode* back = dummyHead;
    for(int i = 0;i< n;i++)
    {
        front = front->next;
    }
    while(front->next)
    {
        front= front->next;
        back = back->next;
    }
    back->next = back->next->next;
    return dummyHead->next;
}

//2.2.7
int indexFromEndX = 0;
ListNode* swapper;
ListNode* removeNthFromEndIII(ListNode* head, int n)
{
    if(head->next != NULL)
 
    {
        removeNthFromEndIII(head->next, n);
    }
    if(head->next == NULL)
    {
        return NULL;
    }
    else
    {
        indexFromEndX ++;
    }
    
    if(indexFromEndX == n - 1)
    {
        swapper=head;
    }
    
    if(indexFromEndX == n + 1)
    {
        head->next = swapper;
    }
    return head;
}

//2.2.8
ListNode *swapPairs(ListNode *head) {
    ListNode* dummyHead = new ListNode(-1);
    dummyHead->next = head;
    ListNode* tmp = dummyHead;
    while(tmp && tmp->next && tmp->next->next)
    {
        ListNode* nextHead = tmp->next->next->next;
        ListNode* first = tmp->next;
        ListNode* second = tmp->next->next;
        tmp->next = second;
        second->next = first;
        first->next = nextHead;
        tmp = first;
    }
    return dummyHead->next;
}

//2.2.9
ListNode *reverseKGroupII(ListNode *head, int k) {
    ListNode* dummyHead = new ListNode(-1);
    dummyHead->next = head;
    ListNode* firstTail = dummyHead;
    ListNode* secondHead = dummyHead->next;
    ListNode* secondTail = dummyHead;
    ListNode* thirdHead = NULL;
    while(secondTail)
    {
        for(int i = 0;i<k;i++)
        {
            if(!secondTail)
            {
                break;
            }
            secondTail = secondTail->next;
        }
        if(!secondTail)
        {
            break;
        }
        
        thirdHead = secondTail->next;
        secondTail->next = NULL;
        reverse(secondHead);
        firstTail->next = secondTail;
        secondHead->next = thirdHead;
        firstTail = secondHead;
        secondHead = thirdHead;
        secondTail = firstTail;
        thirdHead = NULL;
    }
    return dummyHead->next;
}

//2.2.10
RandomListNode *copyRandomList(RandomListNode *head) {
    unordered_map<RandomListNode*, RandomListNode*> map;
    RandomListNode* tmp = head;
    RandomListNode* prev = NULL;
    while(tmp)
    {
        map[tmp] = new RandomListNode(tmp->label);
        if(prev)
        {
            map[prev]->next = map[tmp];
        }
        prev = tmp;
        tmp = tmp->next;
    }
    tmp = head;
    while(tmp)
    {
        map[tmp]->random = map[tmp->random];
        tmp = tmp->next;
    }
    return map[head];
}

RandomListNode *copyRandomListII(RandomListNode *head) {
    while(!head) return NULL;
    
    RandomListNode* tmp = head;
    while(tmp)
    {
        RandomListNode* newNode = new RandomListNode(tmp->label);
        newNode->next = tmp->next;
        tmp->next = newNode;
        tmp= newNode->next;
    }
    tmp = head;
    while(tmp)
    {
        if(tmp->random)
        {
            tmp->next->random = tmp->random->next;
        }
        tmp = tmp->next->next;
    }
    tmp = head;
    RandomListNode* res = tmp->next;
    tmp = head;
    RandomListNode* tmp2 = res;
    while(tmp)
    {
        tmp->next = tmp->next->next;
        tmp = tmp->next;
        if(tmp)
        {
            tmp2->next = tmp2->next->next;
            tmp2 = tmp2->next;
        }
    }
    return res;
}

//2.2.11
bool hasCycle(ListNode *head) {
    if(!head) return false;
    ListNode* fast = head;
    ListNode* slow = head;
    while(fast->next && fast->next->next)
    {
        fast = fast->next->next;
        slow = slow->next;
        if(fast == slow)
            return true;
    }
    return false;
}

//2.2.12
ListNode *detectCycle(ListNode *head) {
    if(!head) return NULL;
    ListNode* fast = head;
    ListNode* slow = head;
    while(fast->next && fast->next->next)
    {
        fast = fast->next->next;
        slow = slow->next;
        if(slow == fast) break;
    }
    if(!(fast->next && fast->next->next))
    {
        return NULL;
    }
    fast = head;
    while(fast != slow)
    {
        fast = fast->next;
        slow = slow->next;
    }
    return fast;
}

//2.2.13
void reorderList(ListNode *head) {
    if(!head) return;
    ListNode* fast = head;
    ListNode* slow = head;
    while(fast->next && fast->next->next)
    {
        fast = fast->next->next;
        slow = slow->next;
    }
    ListNode* secondHead = slow->next;
    slow->next = NULL;
    secondHead = reverse(secondHead);
    ListNode* dummyHead = new ListNode(-1);
    ListNode* tmp = dummyHead;
    while(head && secondHead)
    {
        tmp->next = head;
        head = head->next;
        tmp = tmp->next;
        tmp->next = secondHead;
        secondHead = secondHead->next;
        tmp = tmp->next;
    }
    if(head) tmp->next = head;
    if(secondHead) tmp->next = secondHead;
    return;
}





//2.2.14

class LRUCache{
public:
    LRUCache(int capacity) {
        cap = capacity;
        dummyHead = new DoubleLinkedListNode(-1, -1);
        dummyTail = new DoubleLinkedListNode(-1, -1);
        dummyHead->next = dummyTail;
        dummyTail->prev = dummyHead;
        count = 0;
        
    }
    
    int get(int key) {
        if(map.find(key) != map.end())
        {
            DoubleLinkedListNode* target = map[key];
            target->prev->next = target->next;
            target->next->prev = target->prev;
            
            target->next = dummyHead->next;
            dummyHead->next->prev = target;
            target->prev = dummyHead;
            dummyHead->next = target;
            return target->val;
        }
        return -1;
        
    }
    
    void set(int key, int value) {
        if(cap == 0) return;
        if(map.find(key) == map.end())
        {
            if(count == cap)
            {
                int keyToRemove = dummyTail->prev->key;
                dummyTail->prev = dummyTail->prev->prev;
                dummyTail->prev->next = dummyTail;
                map.erase(keyToRemove);
                count--;
                //need call delete method on the last item.
            }

            DoubleLinkedListNode* tmp = new DoubleLinkedListNode(key, value);
            map[key] = tmp;
            count++;
            dummyHead->next->prev = tmp;
            tmp->next = dummyHead->next;
            dummyHead->next = tmp;
            tmp->prev = dummyHead;
        }
        else
        {
            DoubleLinkedListNode* target = map[key];
            target->val = value;
            target->prev->next = target->next;
            target->next->prev = target->prev;
            
            target->next = dummyHead->next;
            dummyHead->next->prev = target;
            target->prev = dummyHead;
            dummyHead->next = target;

        }
        
    }
private:
    int cap;
    int count;
    unordered_map<int, DoubleLinkedListNode*> map;
    DoubleLinkedListNode* dummyHead;
    DoubleLinkedListNode* dummyTail;
    
};



//3.1
bool isNumAlph(char c)
{
    if((c>='0' && c<='9') ||
       (c>='A' && c<='Z') ||
       (c>='a' && c<='z'))
        return true;
    else
        return false;
}

bool isPalindrome(string s) {
    int left = 0;
    int right = (int)s.length()-1;
    while(left < right)
    {
        while(left < s.length() && !isNumAlph(s[left])) left++;
        while(right >=0 && !isNumAlph(s[right])) right --;
        if(left<=right)
        {
            
            char a = tolower(s[left]);
            char b = tolower(s[right]);
            if(a == b)
            {
                left++;
                right--;
            }
            else
            {
                return false;
            }
        }
    }
    return true;
}

//3.2
char *strStrII(char *haystack, char *needle) {
    int m = (int)strlen(needle);
    int n = (int)strlen(haystack);
    if(m>n) return NULL;
    if(m==0  && n==0) return haystack;
    if(m==0) return haystack;
    
    int tracker[m];
    tracker[0] = -1;
    int j = tracker[0];
    for(int i = 1;i < strlen(needle);i++)
    {
        while(j > -1 && needle[j+1] != needle[i] )
        {
            j = tracker[j];
        }
        if(needle[j+1] == needle[i])
        {
            j++;
        }
        tracker[i] = j;
    }
    
    j = -1;
    int i = 0;
    while(i<n)
    {
        while(j != -1 && needle[j+1] != haystack[i])
        {
            j = tracker[j];
        }
        if(needle[j+1] == haystack[i])
        {
            j++;
        }
        if(j == m-1)
        {
            return haystack+i-j;
        }
        i++;
    }
    return NULL;
}

//3.3
int atoiII(const char *str) {
    while(*str == ' ') str++;
    int sign = 1;
    int res = 0;
    if(str[0] == '-')
    {
        sign = -1;
        str = str+1;
    }
    else if(str[0] == '+')
    {
        sign = 1;
        str = str+1;
    }
    
    int n = (int)strlen(str);
    for(int i = 0;i< n;i++)
    {
        if(str[i] == ' ') continue;
        int tmp = str[i]-'0';
        if(tmp < 0 || tmp > 9) return 0;
        if(sign ==1)
        {
            if(INT_MAX /10 > res && INT_MAX - res* 10 > tmp)
            {
                res = res * 10 + tmp;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            if(INT_MIN /10 < res && INT_MIN + res*10 < (-tmp))
            {
                res = res *10 + tmp;
            }
            else
            {
                return 0;
            }
        }
    }
    res = res * sign;
    return res;
}


//3.4
string addBinary(string a, string b) {
    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());
    int n = (int)a.length();
    int m = (int)b.length();
    int t = max(n,m);
    int carry = 0;
    string res;
    for(int i = 0;i< t;i++)
    {
        int ai = i<n? a[i] - '0' : 0;
        int bi = i<m? b[i] - '0' : 0;
        int tmp = ai+ bi+ carry;
        carry = tmp / 2;
        tmp = tmp % 2;
        res.insert(res.begin(), tmp+'0');
    }
    if(carry ==1) res.insert(res.begin(), '1');
    return res;
}

//3.6
string longestPalindromeII(string s) {
    if(s.length() <=1 ) return s;
    int maxLength = 1;
    string res;
    int left = -1;
    int right = -1;
    for(int i = 0;i< s.length();i++)
    {
        left = i;
        right = i;
        while(s[left] == s[right] && left >=0 && right< s.length())
        {
            if(right-left+1 > maxLength)
            {
                maxLength = right-left+1;
                res = s.substr(left, right-left+1);
            }
            left--;
            right++;
        }
        left = i;
        right = i+1;
        while(s[left] == s[right] && left >=0 && right< s.length())
        {
            if(right-left+1 > maxLength)
            {
                maxLength = right-left+1;
                res = s.substr(left, right-left+1);
            }
            left--;
            right++;
        }
    }
    return res;
}

/**************************************
 *
 *
 * DP solution for longest palindrome substr
 *
 *
 ***************************************/
string longestPalindromeIII(string s) {
    int n = (int)s.length();
    bool tracker[n][n];
    fill_n(&tracker[0][0], n*n, false);
    int maxLength = 1;
    int start = 0;
    for(int j = 0;j<n;j++)
    {
        tracker[j][j] = true;
        for(int i = j-1;i>=0;i-- )
        {
            tracker[i][j] = (s[i] == s[j] && (j-i <2 || tracker[i+1][j-1]));
            if(tracker[i][j] && maxLength < (j-i+1))
            {
                maxLength = j-i+1;
                start = i;
            }
        }
    }
    return s.substr(start, maxLength);
}

//3.6 regular expression
bool isMatchX(const char *s, const char *p) {
    if(*s == '\0' && *p == '\0') return true;
    if(*p == '\0') return *s == '\0';
    if(*s == '\0')
    {
        if(*(p+1) == '*')
        {
            return isMatchX(s, p+2);
        }
        else
        {
            return false;
        }
    }
    if(*(p+1) == '*')
    {
        if(*p == '.' || *p == *s)
        {
            return isMatchX(s, p+2) || isMatchX(s+1, p);
        }
        else
        {
            return isMatchX(s, p+2);
        }
    }
    else
    {
        if(*p == '.' || *p == *s)
        {
            return isMatchX(s+1, p+1);
        }
        else
        {
            return false;
        }
    }
}

//3.7 wild card matching
bool isMatchW(const char *s, const char *p) {
    if(*s == '\0' && *p == '\0') return true;
    if(*p == '\0') return *s == '\0';
    if(*s == '\0')
    {
        if(*p == '*')
        {
            return isMatchX(s, p+1);
        }
        else
        {
            return false;
        }
    }
    if(*p == '*')
    {
        return isMatchW(s, p+1) || isMatchW(s+1, p);
    }
    else if(*p == '?' || *p ==*s)
    {
        return isMatchW(s+1, p+1);
    }
    else
    {
        return false;
    }
}

//3.8
string longestCommonPrefix(vector<string> &strs) {
    if(strs.size() ==0) return "";
    TrieNode* trie = new TrieNode();
    for(vector<string>::iterator it = strs.begin();it != strs.end();it++)
    {
        int n = (int) (*it).length();
        TrieNode* tmp = trie;
        for(int i = 0;i < n; i++)
        {
            if(tmp->val[(*it)[i]] == NULL)
            {
                tmp->val[(*it)[i]] = new TrieNode();
            }
            tmp->count[(*it)[i]] ++;
            tmp = tmp->val[(*it)[i]];
        }
    }
    
    TrieNode* tmp = trie;
    int n = (int)strs.size();
    string res;
    while(tmp)
    {
        TrieNode* runner = NULL;
        for(int i = 0;i< 256;i++)
        {
            if(tmp->val[i] && tmp->count[i] == n)
            {
                res += i;
                runner = tmp->val[i];
                break;
            }
        }
        tmp = runner;
        
    }
    return res;
}

//3.12

string countAndSay(int n) {
    if(n <1) return "";
    string res = "1";
    for(int i = 1;i<n;i++)
    {
        string tmp;
        int count = 0;
        char c = '\0';
        for(int j = 0;j<res.length();j++)
        {
            if(c != res[j] && count != 0)
            {
                tmp += to_string(count);
                tmp += c;
                count = 1;
                c = res[j];
            }
            else
            {
                c = res[j];
                count++;
            }
            
            if(j == res.length()-1)
            {
                tmp+=to_string(count);
                tmp+=c;
            }
        }
        
        res = tmp;
        tmp = "";
    }
    return res;
}

//3.13
vector<string> anagrams(vector<string> &strs) {
    unordered_map<string, vector<string>> tracker;
    for(vector<string>::iterator it = strs.begin();it!=strs.end();it++)
    {
        string tmp = *it;
        sort(tmp.begin(), tmp.end());
        tracker[tmp].push_back(*it);
    }
    
    vector<string> res;
    for(auto it = tracker.cbegin();it != tracker.cend();it++)
    {
        if(it->second.size() >1)
        {
            res.insert(res.end(), it->second.begin(), it->second.end());
        }
    }
    return res;
}

//3.14
string simplifyPath(string path) {
    vector<string> res;
    if(path[path.size()-1] != '/')
    {
        path = path + '/';
    }
    int n = (int)path.length();
    int start = -1;
    for(int i = 0;i<n;i++)
    {
        if(path[i] == '/')
        {
            string tmp = path.substr(start+1, i-start-1);
            if(tmp == "..")
            {
                if(res.size() >0)
                {
                    res.pop_back();
                }
            }
            else if(tmp != "." && !tmp.empty())
            {
                res.push_back(tmp);
            }
            start = i;
        }
    }
    
    //stringstream s;
    if(res.empty())
    {
        cout<<"/";
    }
    else{
        for(auto it:res)
        {
            cout<<"/"<<it;
        }
    }
    return path;
}

//3.14
int lengthOfLastWord(const char *s) {
    int res = 0;
    int n = (int)strlen(s);
    for(int i = 0;i<n-1;i++)
    {
        if(s[i] !=  ' ')
        {
            res ++;
        }
        else if( i!= n-1 && s[i+1] != ' ')
        {
            res = 0;
        }
    }
    return res;
}

//4.1.1
bool isValid(string s) {
    stack<char> tracker;
    for(int i =0;i<s.length();i++)
    {
        if(s[i] == '(' || s[i] == '{' || s[i] == '[')
        {
            tracker.push(s[i]);
        }
        if(s[i] == ')')
        {
            if(!tracker.empty() && tracker.top() == '(')
            {
                tracker.pop();
            }
            else
            {
                return false;
            }
        }
        else if (s[i] == '}')
        {
            if(!tracker.empty() && tracker.top() == '{')
            {
                tracker.pop();
            }
            else
            {
                return false;
            }
        }
        else if(s[i] == ']')
        {
            if(!tracker.empty() && tracker.top() == '[')
            {
                tracker.pop();
            }
            else
            {
                return false;
            }
            
        }
    }
    return tracker.empty();
}

//4.1.2
int longestValidParenthesesII(string s) {
    stack<int> p;
    int last = -1;
    int res = 0;
    for(int i = 0;i<s.length();i++)
    {
        if(s[i] == '(')
        {
            p.push(i);
        }
        else
        {
            if(s[i] == ')')
            {
                if(p.empty())
                {
                    last = i;
                }
                else
                {
                    p.pop();
                    if(p.empty())
                    {
                        res = max(res, i - last);
                    }
                    else
                    {
                        res = max(res, i - p.top());
                    }
                }
            }
        }
    }
    return res;
}

//4.1.3
int largestRectangleAreaII(vector<int> &height) {
    int res = 0;
    stack<int> p;
    height.push_back(0);
    for(int i = 0;i< height.size();i++)
    {
        if(p.empty() || height[i] > height[p.top()])
        {
            p.push(i);
        }
        else
        {
            while(!p.empty() && height[i]<= height[p.top()])
            {
                int right = p.top();
                p.pop();
                // very tricky here to calculate the length.
                int length = p.empty()? i: i-p.top()-1;
                res = max(res, length*height[right]);
            }
            p.push(i);
        }
    }
    return res;
}

//4.1.4
// assume valid RPN expression
int evalRPN(vector<string> &tokens) {
    stack<int> num;
    for(int i = 0;i<tokens.size();i++)
    {
        if(tokens[i] == "+")
        {
            int second = num.top();
            num.pop();
            int first = num.top();
            num.pop();
            int result = first + second;
            num.push(result);
            
        }
        else if(tokens[i] == "-")
        {
            int second = num.top();
            num.pop();
            int first = num.top();
            num.pop();
            int result = first - second;
            num.push(result);
            
        }
        else if(tokens[i] == "*")
        {
            int second = num.top();
            num.pop();
            int first = num.top();
            num.pop();
            int result = first * second;
            num.push(result);
            
        }
        else if(tokens[i] == "/")
        {
            int second = num.top();
            num.pop();
            int first = num.top();
            num.pop();
            int result = first / second;
            num.push(result);
            
        }
        else
        {
            num.push(atoi(&tokens[i][0]));
        }
    }
    return num.top();
}

//5.1.1
vector<int> preorderTraversal(TreeNode *root) {
    vector<int> res;
    if(!root) return res;
    stack<const TreeNode*> tracker;
    tracker.push(root);
    while(!tracker.empty())
    {
        const TreeNode* tmp = tracker.top();
        tracker.pop();
        if(tmp->right) tracker.push(tmp->right);
        if(tmp->left) tracker.push(tmp->left);
        res.push_back(tmp->val);
    }
    return res;
}

//5.1.2
vector<int> inorderTraversal(TreeNode *root) {
    vector<int> res;
    if(!root) return res;
    const TreeNode* tmp = root;
    stack<const TreeNode*> tracker;
    while(!tracker.empty() || tmp)
    {
        if(tmp)
        {
            tracker.push(tmp);
            tmp= tmp->left;
        }
        else
        {
            tmp = tracker.top();
            tracker.pop();
            res.push_back(tmp->val);
            tmp = tmp->right;
        }
    }
    return res;
}


//5.2 Morris inorder
//vector<int> inorderTraversal(TreeNode *root) {
//}

//5.1.3
vector<int> postorderTraversal(TreeNode *root) {
    vector<int> res;
    if(!root) return res;
    TreeNode* prev = nullptr;
    TreeNode* cur;
    stack<TreeNode*> tracker;
    tracker.push(root);
    while(!tracker.empty())
    {
        cur = tracker.top();
        if((cur->left && cur->left != prev) &&
           (!cur->right || cur->right != prev))
        {
            cur = cur->left;
            tracker.push(cur);
            
        }
        else if(prev == cur->right || !(cur->right))
        {
            res.push_back(cur->val);
            tracker.pop();
            prev = cur;
        }
        else
        {
            tracker.push(cur->right);
        }
    }
    return res;
}

void inorderTraversalTreeNode(TreeNode* root)
{
    stack<TreeNode*> s;
    bool backTrack = false;
    s.push(root);
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        if(tmp->left && !backTrack)
        {
            s.push(tmp->left);
            continue;
        }
        cout<<tmp->val<<endl;
        s.pop();
        backTrack = true;
        if(tmp->right)
        {
            s.push(tmp->right);
            backTrack = false;
        }
    }
}

//5.1.4
vector<vector<int> > levelOrder(TreeNode *root) {
    vector<vector<int>> res;
    vector<int> tmp;
    if(!root) return res;
    queue<TreeNode*> q;
    q.push(root);
    q.push(NULL);
    while(!q.empty())
    {
        if(q.front() == NULL)
        {
            res.push_back(tmp);
            tmp.clear();
            q.pop();
            if(q.empty()) break;
            else q.push(NULL);
        }
        else{
            tmp.push_back(q.front()->val);
            if(q.front()->left) q.push(q.front()->left);
            if(q.front()->right) q.push(q.front()->right);
            q.pop();
        }
    }
    return res;
    
}

//5.1.5
vector<vector<int> > levelOrderBottom(TreeNode *root) {
    vector<vector<int>> res = levelOrder(root);
    reverse(res.begin(), res.end());
    return res;
}

//5.1.6
vector<vector<int> > zigzagLevelOrderII(TreeNode *root) {
    vector<vector<int>> res;
    if(!root) return res;
    vector<int> tmp;
    stack<TreeNode*> first;
    stack<TreeNode*> second;
    bool direction = true;
    first.push(root);
    while(!first.empty())
    {
        tmp.push_back(first.top()->val);
        if(direction)
        {
            if(first.top()->left) second.push(first.top()->left);
            if(first.top()->right) second.push(first.top()->right);
        }
        else
        {
            if(first.top()->right) second.push(first.top()->right);
            if(first.top()->left) second.push(first.top()->left);
        }
        first.pop();
        if(first.empty())
        {
            res.push_back(tmp);
            tmp.clear();
            swap(first, second);
            direction = !direction;
        }
    }
    return res;
}

//5.1.8
bool isSameTree(TreeNode *p, TreeNode *q) {
    if(!p && !q) return true;
    if(!p || !q) return false;
    if(p->val != q->val) return false;
    else return isSameTree(p->left, q->left) &&
        isSameTree(p->right, q->right);
}

//5.1.9
bool isSymmetricWorker(TreeNode* t1, TreeNode* t2)
{
    if(!t1 && !t2) return true;
    if(t1->val != t2->val)
    {
        return false;
    }
    else
    {
        return isSymmetricWorker(t1->left, t2->right) && isSymmetricWorker(t1->right, t2->left);
    }
}

bool isSymmetric(TreeNode *root) {
    if(!root) return true;
    return isSymmetricWorker(root->left, root->right);

}

//5.1.10
bool is_Balanced(TreeNode *root, int& height)
{
    if(!root)
    {
        height = 0;
        return true;
    }
    int leftHeight = 0;
    int rightHeight = 0;
    bool res = is_Balanced(root->left, leftHeight) && is_Balanced(root->right, rightHeight);
    height = max(leftHeight, rightHeight)+1;
    res = res && abs(leftHeight - rightHeight) <= 1;
    return res;
    
}
bool isBalanced(TreeNode *root) {
    if(!root)
    {
        return true;
    }
    int leftHeight = 0;
    int rightHeight = 0;
    bool res = is_Balanced(root->left, leftHeight) && is_Balanced(root->right, rightHeight);
    res = res && abs(leftHeight - rightHeight) <= 1;
    return res;
    
}

//5.1.11
void flatten(TreeNode *root) {
    if(!root) return;
    stack<TreeNode*> s;
    TreeNode* tmp = new TreeNode(-1);
    s.push(root);
    while(!s.empty())
    {
        TreeNode* p = s.top();
        s.pop();
        if(p->right) s.push(p->right);
        if(p->left) s.push(p->left);
        tmp->right = p;
        tmp ->left = NULL;
        tmp = p;
    }
}

//5.1.12
void connectII(TreeLinkNode *root) {
    if(!root) return;
    TreeLinkNode* tmp= NULL;

    TreeLinkNode* p = root->next;
    while(p && (!p->left && !p->right))
    {
        p = p->next;
    }
    
    if(!p)
        tmp = NULL;
    else if(p->left)
        tmp = p->left;
    else if(p->right)
        tmp = p->right;
    
    if(root->left && root->right)
    {
        root->left->next = root->right;
        root->right->next = tmp;
    }
    else if(root->left)
    {
        root->left->next = tmp;
    }
    else if(root->right)
    {
        root->right->next = tmp;
    }
    
    connectII(root->right);
    connectII(root->left);
}

//5.2.1
TreeNode* buildTreeWorker(vector<int> inorder, int start1, int end1, vector<int>& postorder, int start2, int end2)
{
    if(start1>end1) return NULL;
    int pivot = postorder[end2];
    int index = -1;
    for(int i = start1; i<=end1;i++)
    {
        if(inorder[i] == pivot)
        {
            index = i;
            break;
        }
    }
    int firstLength = index - start1;
    TreeNode* left = buildTreeWorker(inorder, start1, index-1, postorder, start2, start2+firstLength-1);
    TreeNode* right = buildTreeWorker(inorder, index+1, end1, postorder, start2+firstLength, end2-1);
    TreeNode* root = new TreeNode(postorder[end2]);
    root->left = left;
    root->right = right;
    return root;
}

TreeNode *buildTree(vector<int> &inorder, vector<int> &postorder) {
    if(inorder.size() ==0 || inorder.size() != postorder.size())
        return NULL;
    return  buildTreeWorker(inorder, 0, (int)inorder.size()-1, postorder, 0 , (int)postorder.size()-1);
}

//5.2.2
TreeNode* buildTreeWorkerII(vector<int> inorder, int start1, int end1, vector<int>& pre, int start2, int end2)
{
    if(start1>end1) return NULL;
    int pivot = pre[start2];
    int index = -1;
    for(int i = start1; i<=end1;i++)
    {
        if(inorder[i] == pivot)
        {
            index = i;
            break;
        }
    }
    int firstLength = index - start1;
    TreeNode* left = buildTreeWorkerII(inorder, start1, index-1, pre, start2+1, start2+firstLength);
    TreeNode* right = buildTreeWorkerII(inorder, index+1, end1, pre, start2+firstLength+1, end2);
    TreeNode* root = new TreeNode(pre[start2]);
    root->left = left;
    root->right = right;
    return root;
}

TreeNode *buildTreeII(vector<int> &preorder, vector<int> &inorder) {
    if(inorder.size() ==0 || inorder.size() != preorder.size())
        return NULL;
    return  buildTreeWorkerII(inorder, 0, (int)inorder.size()-1, preorder, 0 , (int)preorder.size()-1);
}

//5.3.1
int numTrees(int n) {
    if(n<2) return 1;
    int f[n];
    f[0] = 1;
    f[1] = 1;
    for(int i = 2;i<n;i++)
    {
        for(int j = 1;j<i;j++)
        {
            f[i] = f[j-1] * f[i-j];
        }
    }
    return f[n];
}

//5.3.2
vector<TreeNode*> generateTrees(int left, int right)
{
    if(left > right)
    {
        vector<TreeNode*> tmp = {NULL};
        return tmp;
    }
    vector<TreeNode*> res;
    for(int i = left;i<=right;i++)
    {
        vector<TreeNode*> leftNode = generateTrees(left, i-1);
        vector<TreeNode*> rightNode = generateTrees(i+1, right);
        for(auto l: leftNode)
        {
            for(auto r : rightNode)
            {
                TreeNode* root = new TreeNode(i);
                root->left = l;
                root->right = r;
                res.push_back(root);
            }
        }
    }
    return res;
}
vector<TreeNode *> generateTrees(int n) {
    return generateTrees(1, n);
}

//5.3.3

bool isValidBST(TreeNode* root, int left, int right)
{
    if(!root) return true;
    if(root->val < left || root->val > right) return false;
    return isValidBST(root, left, root->val) && isValidBST(root, root->val, right);
}
bool isValidBST(TreeNode *root) {
    return isValidBST(root, INT_MIN, INT_MAX);
}

//5.3.4
TreeNode* sortedArrayToBST(vector<int>& num, int left, int right)
{
    if(left>right) return NULL;
    int mid = left + (right-left)/2;
    TreeNode* root = new TreeNode(num[mid]);
    root->left = sortedArrayToBST(num, left, mid-1);
    root->right = sortedArrayToBST(num, mid+1, right);
    return root;
}

TreeNode *sortedArrayToBST(vector<int> &num) {
    return sortedArrayToBST(num, 0, (int)num.size()-1);
}

//5.3.5
TreeNode *sortedListToBSTI(ListNode *head) {
    if(!head) return NULL;
    ListNode* fast = head->next;
    ListNode* slow = head;
    while(fast && fast->next && fast->next->next)
    {
        slow = slow->next;
        fast = fast->next->next;
    }
    
    ListNode* mid = slow->next;
    if(!mid)
    {
        return new TreeNode(slow->val);
    }
    else
    {
        ListNode* newHead = mid->next;
        slow->next = NULL;
        TreeNode* root = new TreeNode(mid->val);
        root->left = sortedListToBST(head);
        root->right = sortedListToBST(newHead);
        return root;
    }
}

//5.4.1
int minDepth(TreeNode *root) {
    if(!root) return 0;
    if(!root->left && !root->right) return 1;
    int left = minDepth(root->left);
    int right =  minDepth(root->right);
    int res;
    if(left == 0 && right != 0) res = left;
    else if(right == 0 && left != 0) res = left;
    else if(right != 0 && left != 0) res = min(left, right);
    else res = 0;
    
    return res+1;
}


//5.4.2
int maxDepth(TreeNode *root) {
    if(!root) return 0;
    return 1 + max(maxDepth(root->left), maxDepth(root->right));
}

//5.4.3
bool hasPathSum(TreeNode *root, int sum) {
    if(!root) return false;
    stack<TreeNode*> s;
    s.push(root);
    TreeNode* prev = NULL;
    int t = root->val;
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        if(!s.top()->left && !s.top()->right)
        {
            if(t == sum)
                return true;
        }
        if(tmp ->left && prev != tmp -> left && (!(tmp->right) || prev != tmp->right))
        {
            s.push(tmp->left);
            t += tmp->left->val;
        }
        else if(!tmp->right || prev == tmp->right)
        {
            s.pop();
            prev = tmp;
            t -= tmp->val;
        }
        else
        {
            s.push(tmp->right);
            t += tmp->right->val;
        }

    }
    return false;
}

//5.4.4
vector<vector<int> > pathSum(TreeNode *root, int sum) {
    vector<vector<int>> res;
    if(!root) return res;
    stack<TreeNode*> s;
    s.push(root);
    TreeNode* prev = NULL;
    int t = root->val;
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        if(!s.top()->left && !s.top()->right)
        {
            if(t == sum)
            {
                stack<TreeNode*> copy = s;
                vector<int> v;
                while(!copy.empty())
                {
                    v.push_back(copy.top()->val);
                    copy.pop();
                }
                reverse(v.begin(), v.end());
                res.push_back(v);
            }

        }
        if(tmp ->left && prev != tmp -> left && (!(tmp->right) || prev != tmp->right))
        {
            s.push(tmp->left);
            t += tmp->left->val;
        }
        else if(!tmp->right || prev == tmp->right)
        {
            s.pop();
            prev = tmp;
            t -= tmp->val;
        }
        else
        {
            s.push(tmp->right);
            t += tmp->right->val;
        }
    }
    return res;
}

//5.4.5
int maxSum = INT_MIN;
int workerII(TreeNode *root) {
    if(!root) return 0;
    
    int left = worker(root->left);
    int right = worker(root->right);
    
    int sum = root->val;
    if(left>0) sum+= left;
    if(right>0) sum+=right;
    maxSum = max(sum, maxSum);
    
    return max(left,right)>0 ? root->val + max(left,right): root->val;
}

int maxPathSumII(TreeNode *root) {
    return max(workerII(root), maxSum);
}

//5.4.6
int getInteger(vector<int> input)
{
    int res = 0;
    for(int i = 0;i<input.size();i++)
    {
        res = res * 10 + input[i];
    }
    return res;
}

int sumNumbers(TreeNode *root) {
    if(!root) return 0;
    int res;
    stack<TreeNode*> s;
    s.push(root);
    TreeNode* prev = NULL;
    int t = root->val;
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        if(!s.top()->left && !s.top()->right)
        {
            stack<TreeNode*> copy = s;
            vector<int> v;
            while(!copy.empty())
            {
                v.push_back(copy.top()->val);
                copy.pop();
            }
            reverse(v.begin(), v.end());
            res+= getInteger(v);
        }
        if(tmp ->left && prev != tmp -> left && (!(tmp->right) || prev != tmp->right))
        {
            s.push(tmp->left);
            t += tmp->left->val;
        }
        else if(!tmp->right || prev == tmp->right)
        {
            s.pop();
            prev = tmp;
            t -= tmp->val;
        }
        else
        {
            s.push(tmp->right);
            t += tmp->right->val;
        }
    }
    return res;
}

//6.1
void merge(int A[], int m, int B[], int n) {
    int a = m-1;
    int b = n-1;
    while(a>=0 && b>=0)
    {
        if(B[b] > A[a])
        {
            A[a+b+1] = B[b];
            b--;
        }
        else
        {
            A[a+b+1] = A[a];
            a--;
        }
    }
    
    if(b<0)
    {
        return ;
    }
    else
    {
        while(b>=0)
        {
            A[b] = B[b];
            b--;
        }
    }
    return;
}

//6.2
ListNode *mergeTwoLists(ListNode *l1, ListNode *l2) {
    ListNode* dummyHeader = new ListNode(-1);
    ListNode* tmp = dummyHeader;
    while(l1 && l2)
    {
        if(l1->val < l2->val)
        {
            tmp->next = l1;
            l1 = l1->next;
        }
        else
        {
            tmp->next = l2;
            l2 = l2->next;
        }
        tmp = tmp->next;
    }
    if(l1) tmp->next = l1;
    else tmp->next = l2;
    return dummyHeader->next;
}

//6.3
ListNode* mergeLists(vector<ListNode*>& lists, int start, int end)
{
    if(start > end)
    {
        return NULL;
    }
    if( start == end)
    {
        return lists[start];
    }
    if(start+1 == end)
    {
        return mergeTwoLists(lists[start], lists[end]);
    }
    
    int mid = start + (end - start) / 2;
    
    ListNode* first =mergeLists(lists, start, mid);
    ListNode* second = mergeLists(lists, mid+1, end);
    
    return mergeTwoLists(first, second);
}

ListNode *mergeKLists(vector<ListNode *> &lists) {
    return mergeLists(lists, 0, (int)lists.size()-1);
}

//6.4
ListNode *insertionSortList(ListNode *head) {
    if(!head) return NULL;
    ListNode* dummyHead = new ListNode(INT_MIN);
    ListNode* tmp = dummyHead;
    ListNode* cur = head;
    while(cur)
    {
        head = cur->next;
        while( tmp->next && cur->val >= tmp->next->val)
        {
            tmp = tmp->next;
        }
        cur->next = tmp->next;
        tmp->next = cur;
        tmp = dummyHead;
        cur = head;
    }
    return dummyHead->next;
}

//6.5
ListNode *sortList(ListNode *head) {
    if(!head) return NULL;
    if(!head->next) return head;
    ListNode* fast = head;
    ListNode* slow = head;
    while(fast->next && fast->next->next)
    {
        fast = fast->next->next;
        slow = slow->next;
    }
    ListNode* secondHead = slow->next;
    slow->next = NULL;
    head = sortList(head);
    secondHead = sortList(secondHead);
    ListNode* res = mergeTwoLists(head, secondHead);
    return res;
}

//6.6
int firstMissingPositiveII(int A[], int n) {
    int i = 0;
    while(i<n)
    {
        if((A[i] <1 || A[i]>n) || A[i] == A[A[i]-1])
        {
            i++;
        }
        else
        {
            swap(A, i, A[i]-1);
        }
    }
    for(int i = 0;i< n;i++)
    {
        if(A[i]!= i+1)
            return i+1;
    }
    return n+1;
}

//6.7
void sortColorsII(int A[], int n) {
    int left = 0;
    int right = n -1;
    int tmp = left;
    while(tmp <= right)
    {
        if(A[tmp] == 0)
        {
            swap(A, left++, tmp++);
        }
        else if(A[tmp] == 1)
        {
            tmp ++;
        }
        else if( A[tmp] == 2)
        {
            swap(A, tmp, right--);
        }
    }
}

//7.1
vector<int> searchRangeII(int A[], int n, int target) {
    int left = 0;
    int right = n-1;
    int mid = 0;
    while(left<=right)
    {
        mid = left + (right-left) /2;
        if(A[mid] == target && (mid == 0 || A[mid-1] < target))
        {
            break;
        }
        else if(A[mid]>=target)
        {
            right = mid-1;
        }
        else
        {
            left = mid + 1;
        }
    }
    int leftindex = left<=right? mid:-1;
    left = 0;
    right = n-1;
    while(left<=right)
    {
        mid = left + (right-left) /2;
        if(A[mid] == target && (mid == n-1 || A[mid+1] > target ) )
        {
            break;
        }
        else if(A[mid] <= target)
        {
            left = mid +1;
        }
        else
        {
            right = mid-1;
        }
    }
    int rightIndex = left<=right? mid:-1;
    vector<int> res = {leftindex, rightIndex};
    return res;
}

//7.2
int searchInsert(int A[], int n, int target) {
    if(target < A[0]) return 0;
    if(target > A[n-1]) return n;
    int left = 0;
    int right = n - 1;
    int mid = -1;
    while(left<=right)
    {
        mid = left + (right - left) /2;
        if(A[mid] == target)
            break;
        else if(A[mid] < target)
            left = mid+1;
        else
            right = mid-1;
    }
    int res = -1;
    if(left<=right)
    {
        res = mid;
    }
    else if(target <= A[right])
    {
        res = right;
    }
    else
    {
        res = left;
    }
    return res;
}

//7.3
bool searchMatrix(vector<vector<int> > &matrix, int target) {
    int n = (int)matrix.size();
    int m = (int)matrix[0].size();
    int i = 0;
    int j = m-1;
    while(i<n && j >=0)
    {
        if(matrix[i][j] > target)
        {
            j --;
        }
        else if(matrix[i][j] < target)
        {
            i++;
        }
        else
        {
            return true;
        }
    }
    return false;
}

//8.1
void worker(vector<int> &S, int level, vector<vector<int>>& res, vector<int> tracker)
{
    if(level == S.size())
    {
        res.push_back(tracker);
    }
    worker(S, level+1, res, tracker);
    tracker.push_back(S[level]);
    worker(S, level+1, res, tracker);
    tracker.pop_back();
}

vector<vector<int>> subsets(vector<int> &S) {
    vector<vector<int>> res;
    vector<int> tracker;
    worker(S, 0, res, tracker);
    return res;
}

//8.2
void workerx(vector<int> &S, int level, vector<int>& tmp, vector<vector<int>>& res)
{
    res.push_back(tmp);
    
    for(int i = level;i< S.size();i++)
    {
        if(i == level || S[i] != S[i-1])
        {
            tmp.push_back(S[i]);
            workerx(S, i+1, tmp, res);
            tmp.pop_back();
        }
    }
}

vector<vector<int> > subsetsWithDupII(vector<int> &S) {
    vector<vector<int>> res;
    if(S.size() == 0) return res;
    
    sort(S.begin(), S.end());
    vector<int> tmp;
    workerx(S, 0, tmp, res);
    return res;
}

//8.3
void permuteworker(vector<int>& num, int level, vector<vector<int>>& res)
{
    if(level == num.size())
    {
        res.push_back(num);
        return;
    }
    for(int i = level;i< num.size(); i++)
    {
        swap(num, level, i);
        permuteworker(num, level+1, res);
        swap(num, level, i);
    }
}

vector<vector<int> > permute(vector<int> &num) {
    vector<vector<int>> res;
    sort(num.begin(), num.end());
    permuteworker(num, 0, res);
    return res;
}

//8.4
// the reason I have to use an unordered_set to track which number I have used
// instead of just compare the the item before the current item in the vector
// num, is that, after each swap, the item in vector num is no longer sorted..
void permuteworkerUnique(vector<int>& num, int level, vector<vector<int>>& res)
{
    if(level == num.size())
    {
        res.push_back(num);
        return;
    }
    unordered_set<int> t;
    for(int i = level;i< num.size(); i++)
    {
        if(i == level || t.find(num[i]) == t.end())
        {
            t.insert(num[i]);
            swap(num, level, i);
            permuteworkerUnique(num, level+1, res);
            swap(num, level, i);
        }
    }
}

vector<vector<int> > permuteUniqueII(vector<int> &num) {
    vector<vector<int>> res;
    sort(num.begin(), num.end());
    permuteworkerUnique(num, 0, res);
    return res;
}

//8.5
void combineWorker(int n, int k, int level, vector<vector<int>>& res, vector<int>& tracker)
{
    if(level == k)
    {
        if(tracker.size() == k)
        {
            res.push_back(tracker);
        }
        return;
    }
    for(int i = level+1;i<=n;i++)
    {
        if(tracker.size() == 0 || i > tracker[tracker.size()-1])
        {
            tracker.push_back(i);
            combineWorker(n, k, level+1, res, tracker);
            tracker.pop_back();
        }
    }
}

vector<vector<int> > combine(int n, int k) {
    vector<vector<int>> res;
    vector<int> tracker;
    combineWorker(n, k, 0, res, tracker);
    return res;
    
}


//8.6
void letterCombinationWorker(int level, string digits, string& tmp, vector<string>& res, unordered_map<char, string> tracker)
{
    if(level == digits.length())
    {
        res.push_back(tmp);
        return;
    }
    
    string x = tracker[digits[level]];
    for(int i = 0;i< x.length();i++)
    {
        tmp = tmp + x[i];
        letterCombinationWorker(level+1, digits, tmp, res, tracker);
        tmp.pop_back();
    }
}

vector<string> letterCombinations(string digits) {
    vector<string> res;
    if(digits.length() == 0)
    {
        res.push_back("");
        return res;
    }
    unordered_map<char, string> map;
    map['2']= "abc";
    map['3'] = "def";
    map['4'] = "ghi";
    map['5'] = "jkl";
    map['6'] = "mno";
    map['7'] = "pqrs";
    map['8'] = "tuv";
    map['9'] = "wxyz";
    string tmp;
    letterCombinationWorker(0, digits, tmp, res, map);
    
    return res;
}

//9.1
vector<string> getString(string s)
{
    vector<string> res;
    for(int i = 0;i<s.length();i++)
    {
        char oldChar = s[i];
        string tmp(s);
        for(char c = 'a';c <= 'z';c++)
        {
            if(c != oldChar)
            {
                swap(c, tmp[i]);
                res.push_back(tmp);
                swap(c, tmp[i]);
            }
        }
        //swap(oldChar, s[i]);
    }
    return res;
}

int ladderLength(string start, string end, unordered_set<string> &dict) {
    queue<string> tracker;
    int len = 0;
    unordered_set<string> visited;
    tracker.push(start);
    visited.insert(start);
    tracker.push("");
    while(!tracker.empty())
    {
        string tmp = tracker.front();
        tracker.pop();
        if(tmp == "")
        {
            if(tracker.empty())
            {
                break;
            }
            else
            {
                len++;
                tracker.push("");
            }
        }else
        {
            vector<string> next = getString(tmp);
            for(auto i : next)
            {
                if(i == end)
                {
                    return len+2;
                }
                else
                {
                    if(visited.find(i) == visited.end() && dict.find(i) != dict.end())
                    {
                        visited.insert(i);
                        tracker.push(i);
                    }
                }
            }
        }
    }
    return 0;
}

//9.2
//this method still need some work. if a word
//can be reached by 2 word from previous level
//that will be 2 different path then ...

struct qElement{
    string val;
    qElement* parent;
    qElement(string s) : val(s), parent(NULL) {}
};

vector<vector<string>> findLadders(string start, string end, unordered_set<string> &dict) {
    
    queue<qElement*> tracker;
    vector<vector<string>> res;
    int len = 0;
    bool found = false;
    unordered_set<string> visited;
    tracker.push(new qElement(start));
    visited.insert(start);
    tracker.push(new qElement(""));
    while(!tracker.empty())
    {
        qElement* tmp = tracker.front();
        tracker.pop();
        if(tmp->val== "")
        {
            if(tracker.empty() || found)
            {
                break;
            }
            else
            {
                len++;
                tracker.push(new qElement(""));
            }
        }else
        {
            vector<string> next = getString(tmp->val);
            for(auto i : next)
            {
                if(i == end)
                {
                    found = true;
                    vector<string> swapper;
                    qElement* runner = tmp;
                    while(runner)
                    {
                        swapper.insert(swapper.begin(), runner->val);
                        runner = runner->parent;
                    }
                    //swapper.insert(swapper.begin(), start);
                    swapper.push_back(end);
                    res.push_back(swapper);
                }
                else
                {
                    if(visited.find(i) == visited.end() && dict.find(i) != dict.end())
                    {
                        visited.insert(i);
                        qElement* qe = new qElement(i);
                        qe->parent = tmp;
                        tracker.push(qe);
                    }
                }
            }
        }
    }
    return res;
}

//9.3 surrounding region
void findArea(vector<vector<char>>& board, int k, int l, int m, int n, char old, char newchar){
    
    stack<pair<int,int>> tracker;
    if(board[k][l] == old){
        tracker.push(make_pair(k, l));
    }
    
    while(!tracker.empty()){
        pair<int,int> tmp = tracker.top();
        tracker.pop();
        board[tmp.first][tmp.second] = newchar;
        k = tmp.first;
        l = tmp.second;
        
        if(k>0 && board[k-1][l] == old){
            tracker.push(make_pair( k-1, l));
        }
        if(k+1 <m && board[k+1][l] == old){
            tracker.push(make_pair( k+1, l));
        }
        if(l>0 && board[k][l-1] == old){
            tracker.push(make_pair( k, l-1));
        }
        if(l+1<n && board[k][l+1] == old){
            tracker.push(make_pair( k, l+1));
        }
    }
}
void solveSurroundingRegion(vector<vector<char>> &board) {
    if(&board == NULL) return;
    
    const int m = (int)board.size();
    if(m == 0) return;
    const int n = (int)board[0].size();
    if(n == 0) return;
    
    for(int i = 0;i< m;i++){
        for(int j = 0;j < n;j ++){
            if(board[i][j] == 'O'){
                board[i][j] = 'A';
            }
        }
    }
    
    for(int j = 0; j<n;j++){
        if(board[0][j] == 'A'){
            findArea(board, 0, j, m, n, 'A', 'O');
        }
        if(board[m-1][j] == 'A'){
            findArea(board, m-1, j, m, n, 'A', 'O');
        }
    }
    
    for(int i = 0;i<m;i++){
        if(board[i][0] == 'A'){
            findArea(board, i, 0, m, n, 'A', 'O');
        }
        if(board[i][n-1] == 'A'){
            findArea(board, i, n-1, m, n, 'A', 'O');
        }
    }
    
    for(int i = 0;i< m;i++){
        for(int j = 0;j < n;j ++){
            if(board[i][j] == 'A'){
                board[i][j] = 'X';
            }
        }
    }
}

//10.1
bool isPalindrome(string s, int i, int j)
{
    while(i<=j)
    {
        if(s[i] == s[j])
        {
            i++;
            j--;
        }
        else
        {
            return false;
        }
    }
    return true;
}

void workerPartition(string s, int level, vector<string>& tracker, vector<vector<string>>& res)
{
    if(level == s.length())
    {
        res.push_back(tracker);
        return;
    }
    else
    {
        for(int i = level;i< s.size();i++)
        {
            if(isPalindrome(s, level, i))
            {
                tracker.push_back(s.substr(level, i - level+1));
                workerPartition(s, i+1, tracker, res);
                tracker.pop_back();
            }
        }
    }
}

vector<vector<string>> partition(string s)
{
    vector<string> tracker;
    vector<vector<string>> res;
    workerPartition(s, 0, tracker, res);
    return res;
    
}

//10.2
int uniquePaths(int m, int n) {
    if(m == 0 && n ==0) return 0;
    vector<vector<int>> tracker(m, vector<int>(n, 0));
    for(int i =0;i<m;i++)
    {
        tracker[i][0] = 1;
    }
    for (int j = 0;j<n;j++)
    {
        tracker[0][j] = 1;
    }
    for(int i = 1;i<m;i++)
    {
        for(int j = 1;j<n;j++)
        {
            tracker[i][j] = tracker[i-1][j] + tracker[i][j-1];
        }
    }
    return tracker[m-1][n-1];
}

//10.3
int uniquePathsWithObstacles(vector<vector<int> > &obstacleGrid) {
    int m = (int)obstacleGrid.size();
    int n = (int)obstacleGrid[0].size();
    
    vector<vector<int>> tracker(m, vector<int>(n, INT_MAX));
    if(obstacleGrid[0][0] == 0)
    {
        tracker[0][0] = 1;
    }
    for(int i = 1;i<m ;i++)
    {
        if(obstacleGrid[i][0] == 0)
        {
            tracker[i][0] = tracker[i-1][0];
        }
    }
    for(int j = 1;j<n;j++)
    {
        if(obstacleGrid[0][j] == 0)
        {
            tracker[0][j] = tracker[0][j-1];
        }
    }
    
    for(int i = 1;i<m;i++)
    {
        for(int j = 1;j<n;j++)
        {
            if(obstacleGrid[i][j] ==0)
            {
                if((tracker[i-1][j] == INT_MAX) && (tracker[i][j-1] == INT_MAX))
                {
                    tracker[i][j] = INT_MAX;
                }
                else if(tracker[i-1][j] == INT_MAX)
                {
                    tracker[i][j] = tracker[i][j-1];
                }
                else if(tracker[i][j-1] == INT_MAX)
                {
                    tracker[i][j] = tracker[i-1][j];
                }
                else
                {
                    tracker[i][j] = tracker[i-1][j]+tracker[i][j-1];
                }
            }
        }
    }
    for(int i = 0;i<m;i++)
    {
        for(int j = 0;j<n;j++)
        {
            cout<< tracker[i][j] << " ";
        }
        cout<<endl;
    }
    
    return tracker[m-1][n-1] == INT_MAX? 0: tracker[m-1][n-1];
}

//10.4
bool checkQueeens(vector<int> queens, int k)
{
    int n = (int)queens.size();
    
    for(int i = 0;i< n;i++)
    {
        if(k == queens[i])
        {
            return false;
        }
        if(abs(n - i) == abs(queens[i] - k))
        {
            return false;
        }
    }
    return true;
}

vector<string> makeMatrix(vector<int> queens, int n)
{
    vector<string> res;
    for(int i = 0;i<n;i++)
    {
        string tmp;
        for(int j=0;j<n;j++)
        {
            if(queens[i] == j)
            {
                tmp = tmp + 'Q';
            }
            else
            {
                tmp = tmp + '.';
            }
        }
        res.push_back(tmp);
    }
    return res;
}



void worker(int n, vector<int>& queens, vector<vector<string>>& res)
{
    if(queens.size() == n)
    {
        res.push_back(makeMatrix(queens, n));
        return;
    }
    else
    {
        for(int j = 0;j<n;j++)
        {
            if(checkQueeens(queens, j))
            {
                queens.push_back(j);
                worker(n, queens, res);
                queens.pop_back();
            }
        }
    }
}

vector<vector<string> > solveNQueens(int n) {
    vector<vector<string>> res;
    vector<int> queens;
    worker(n, queens, res);
    return res;
}

//10.5
void worker(int n, vector<int>& queens, int& res)
{
    if(queens.size() == n)
    {
        res++;
        return;
    }
    else
    {
        for(int j = 0;j<n;j++)
        {
            if(checkQueeens(queens, j))
            {
                queens.push_back(j);
                worker(n, queens, res);
                queens.pop_back();
            }
        }
    }
}

int totalNQueens(int n) {
    int res = 0;
    vector<int> queens;
    worker(n, queens, res);
    return res;
}

//10.6
int getInt(string s)
{
    int res;
    for(int i = (int)s.length()-1;i>=0;i--)
    {
        res += s[i] - '0';
    }
    return res;
}

void worker(string s, vector<int>& tracker, vector<string>& res )
{
    if(tracker.size() == 3)
    {
        int a1 =getInt(s.substr(0, tracker[0]));
        int a2 = getInt(s.substr(tracker[0], tracker[1]-tracker[0]));
        int a3 = getInt(s.substr(tracker[1], tracker[2]-tracker[1]));
        int a4 = getInt(s.substr(tracker[2], tracker.size() - tracker[2]));
        if((a1<=255 && a2 <=255 && a3 <=255 && a4<= 255) &&
           (a1 != 0 || a2!=0 || a3!= 0 || a4 !=0))
        {
            string tmp = s;
            for(int i = (int)tracker.size()-1;i>=0;i--)
            {
                tmp.insert(tracker[i], ".");
            }
            res.push_back(tmp);
        }
    }
    else
    {
        int level = 0;
        if(tracker.size()>0)
        {
            level = tracker[tracker.size()-1];
        }
        for(int i = level+1;i<level+4;i++)
        {
            int a =getInt(s.substr(level, i-level));
            bool valid = true;
            
            if(a>255) valid = false;
            if (a == 0 && i-level > 1) valid = false;
            if(a>0 && s[level] == '0') valid = false;
                
                
            if(valid)
            {
                tracker.push_back(i);
                worker(s, tracker, res);
                tracker.pop_back();
            }
        }
    }
}
vector<string> restoreIpAddresses(string s) {
    vector<string> res;
    if(s.length()> 12)
        return res;
    vector<int> tracker;
    worker(s, tracker, res);
    return res;
    
}

//10.6
void dfs(string s, size_t start, size_t step, string ip, vector<string> &result)
{
    if (start == s.size() && step == 4)
    {
        // 找到一个合法解
        ip.resize(ip.size() - 1);
        result.push_back(ip);
        return;
    }
    
    if (s.size() - start > (4 - step) * 3)
        return; // 剪枝
    
    if (s.size() - start < (4 - step))
        return; // 剪枝
    
    int num = 0;
    
    for (size_t i = start; i < start + 3; i++)
    {
        num = num * 10 + (s[i] - '0');
        // 当前结点合法，则继续往下递归
        if (num <= 255)
        {
            ip += s[i];
            dfs(s, i + 1, step + 1, ip + '.', result);
        }
        
        // this is a key!!!!!!!!!!!!!!
        if (num == 0)
            break; // 不允许前缀 0，但允许单个 0
    }
}

vector<string> restoreIpAddressesII(string s)
{
    vector<string> result;
    string ip; // 存放中间结果
    dfs(s, 0, 0, ip, result);
    return result;
}


//10.7
void workerSum(vector<int>& candidates, int level, int& sum, int target, vector<int>& tmp, vector<vector<int>>& res)
{
    if(sum == target)
    {
        res.push_back(tmp);
    }
    else if(sum > target)
    {
        return;
    }
    else
    {
        for(int i = level;i<candidates.size();i++)
        {
            sum += candidates[i];
            tmp.push_back(candidates[i]);
            workerSum(candidates, i, sum, target, tmp, res);
            tmp.pop_back();
            sum -= candidates[i];
        }
    }
}

vector<vector<int> > combinationSum(vector<int> &candidates, int target) {
    sort(candidates.begin(), candidates.end());
    vector<int> tmp;
    vector<vector<int>> res;
    int sum = 0;
    workerSum(candidates, 0, sum, target, tmp, res);
    return res;
    
}

//10.8
void workerII(vector<int>& num, int gap, vector<int>& tmp, int level, vector<vector<int>>& res)
{
    if(0 == gap)
    {
        res.push_back(tmp);
        return;
    }
    else
    {
        int prev = -1;
        for(int i = level;i< num.size();i++)
        {
            if(num[i] == prev)
            {
                continue;
            }
            if(gap < num[i])
            {
                return;
            }
            prev = num[i];
            tmp.push_back(num[i]);
            workerII(num, gap - num[i], tmp, i+1, res);
            tmp.pop_back();
        }
    }
}

vector<vector<int> > combinationSum2(vector<int> &num, int target) {
    
    sort(num.begin(), num.end());
    vector<vector<int>> res;
    vector<int> tmp;
    workerII(num, target, tmp, 0, res);
    return res;
}

//10.9
void worker(int left, int right, int n, string tmp, vector<string>& res)
{
    if(left >n || right >n)
    {
        return;
    }
    if(left == n && right == n)
    {
        res.push_back(tmp);
    }
    else
    {
        if(left >right)
        {
            worker(left+1, right, n, tmp + "(", res);
            worker(left, right+1, n, tmp + ")", res);
        }
        else if(right == left)
        {
            worker(left+1, right, n, tmp + "(", res);
        }
        else
        {
            return;
        }
    }
}

vector<string> generateParenthesis(int n) {
    vector<string> res;
    worker(0,0,n,"", res);
    return res;
}

//10.11
bool existWorker(vector<vector<char>> board, string word, string tmp, int i, int j, vector<vector<bool>> tracker)
{
    if(tmp.size() > word.size())
    {
        return false;
    }
    if(word.size() == tmp.size())
    {
        return true;
    }
    else
    {
        int m = (int)board.size();
        int n = (int)board[0].size();
        if(!tracker[i][j] && word[tmp.size()] == board[i][j])
        {
            tracker[i][j] = true;
                
            if(((i >0 && existWorker(board, word, tmp + board[i][j], i-1, j, tracker)) ||
            (j >0 && existWorker(board, word, tmp + board[i][j], i, j-1, tracker)) ||
            (i < m-1 && existWorker(board, word, tmp + board[i][j], i+1, j, tracker)) ||
                    (j < n-1 && existWorker(board, word, tmp + board[i][j], i, j+1, tracker))))
            {
                return true;
            }
            else
            {
                tracker[i][j] = false;
                return false;
            }
        }
        else
        {
            return false;
        }
    }
}

bool exist(vector<vector<char> > &board, string word) {
    int m = (int)board.size();
    int n = (int)board[0].size();
    vector<vector<bool>> tracker(m, vector<bool>(n, false));
    for(int i = 0;i<m;i++)
    {
        for(int j = 0;j<n;j++)
        {
            if(existWorker(board, word, "", i, j, tracker))
            {
                return true;
            }
        }
    }
    return false;
}

static bool dfs(const vector<vector<char> > &board, const string &word,
                int index, int x, int y, vector<vector<bool> > &visited) {
    if (index == word.size())
        return true; // 收敛条件
    if (x < 0 || y < 0 || x >= board.size() || y >= board[0].size())
        return false; // 越界,终止条件
    if (visited[x][y]) return false; // 已经访问过,剪枝
    if (board[x][y] != word[index]) return false; // 不相等,剪枝
    visited[x][y] = true;
    bool ret = dfs(board, word, index + 1, x - 1, y, visited) ||
    dfs(board, word, index + 1, x + 1, y, visited) ||
    dfs(board, word, index + 1, x, y - 1, visited)||
    dfs(board, word, index + 1, x, y + 1, visited);
    visited[x][y] = false;
    return ret;
}

bool existII(vector<vector<char> > &board, string word) {
    int m = (int)board.size();
    int n = (int)board[0].size();
    vector<vector<bool>> tracker(m, vector<bool>(n, false));
    for(int i = 0;i<m;i++)
    {
        for(int j = 0;j<n;j++)
        {
            if(dfs(board, word, 0, i, j, tracker))
            {
                return true;
            }
        }
    }
    return false;
}

//11.1
double pow(double x, int n) {
    int sign = 1;
    if(n<0)
    {
        sign = -1;
        n = -n;
    }
    double tmp;
    if( n == 0)
        return 1;
    else if( n == 1)
        return x;
    else
    {
        tmp = pow(x, n/2);
        if(n & 1)
        {
            tmp = tmp * tmp * x;
        }
        else
        {
            tmp = tmp * tmp;
        }
    }
    if(sign<0)
        return 1/tmp;
    else
        return tmp;
}

//11.1
int sqrtLocal(int x) {
    if(x < 0) return -1;
    if( x ==1) return 1;
    double left = 0;
    double right = x;
    while(right - left > 0.0001)
    {
        double mid = left + (right-left)/2;
        if(mid* mid < x)
        {
            left = mid;
        }
        else
        {
            right = mid;
        }
    }
    return (int)left;
}

//12.1
bool canJump(int A[], int n) {
    int lastIndex = 0;
    for(int i = 0;i<n;i++)
    {
        if(lastIndex >= i)
        {
            lastIndex = max(lastIndex, i+ A[i]);
        }
    }
    return lastIndex>= n-1;
}

//12.2
// this is dynamic programming. O(n^2)
int jumpII(int A[], int n) {
    int tracker[n];
    fill_n(&tracker[0],n,INT_MAX);
    tracker[0] =0;
    for(int i = 1;i<n;i++)
    {
        int tmp = tracker[i];
        for(int j = i-1;j>=0;j--)
        {
            if(i-j >A[j] && A[j] != INT_MAX)
            {
                tmp = min(tmp, A[j] +1);
            }
            
        }
        A[i] = tmp;
    }
    return A[n-1];
}

//this is a O(n)
int jumpIII(int A[], int n)
{
    if(n<2) return 0;
    int left = 0;
    int right = 0;
    int step = 0;
    while(left <=right && right < n)
    {
        int old_right = right;
        step++;
        for(int i = left;i<=old_right;i++)
        {
            right = max(right, i+ A[i]);
            if(right>=n-1) return step;
        }
        left = old_right+1;
    }
    return step;
}

//12.3
int maxProfit(vector<int> &prices) {
    int res = 0;
    int c = 0;
    for(int i = 1;i<(int)prices.size();i++)
    {
        c = max(prices[i] - prices[i-1], c + prices[i] - prices[i-1]);
        res = max ( c, res);
    }
    return res;
}

//12.4
int maxProfitII(vector<int> &prices) {
    int res = 0;
    for(int i = 1;i<(int)prices.size();i++)
    {
        if(prices[i] - prices[i-1] >0)
        {
            res += prices[i]-prices[i-1];
        }
    }
    return res;
}

//12.5
int lengthOfLongestSubstring(string s) {
    int tracker[256];
    fill_n(&tracker[0], 256, -1);
    int startIndex = -1;
    int res = 0;
    for(int i = 0;i<(int)s.size();i++)
    {
        if(tracker[s[i]] > startIndex)
        {
            startIndex = tracker[s[i]];
        }
        res= max(res, i - startIndex);
        tracker[s[i]] = i;
    }
    return res;
}

//12.6
int maxArea(vector<int> &height) {
    int res = 0;
    int left = 0;
    int right = (int)height.size()-1;
    
    while(left<right)
    {
        res = max(res, (right-left) * min(height[left], height[right]));
        if(height[left] < height[right])
        {
            left++;
        }
        else
        {
            right--;
        }
    }
    return res;
}

//13.1
int minimumTotal(vector<vector<int>>& triangle)
{
    vector<vector<int>> tracker;
    int n = (int)triangle.size();
    if(n ==0) return 0;
    tracker.push_back(triangle[0]);
    for(int i = 1;i<n;i++)
    {
        vector<int> tmp;
        for(int j = 0;j< triangle[i].size();j++)
        {
            int x;
            if(j ==0) x = tracker[i-1][0];
            else if(j == triangle[i].size()-1) x = tracker[i-1][j-1];
            else x = min (tracker[i-1][j-1], tracker[i-1][j]);
                
            tmp.push_back(x + triangle[i][j]);
        }
        tracker.push_back(tmp);
    }
    
    int res = INT_MAX;
    for(int i = 0;i< tracker[n-1].size();i++)
    {
        res = min (res, tracker[n-1][i]);
    }
    return res;
}

//13.2
int maxSubArray(int A[], int n) {
    if(n <1) return 0;
    int res = INT_MIN;
    int f = 0;
    int left = -1;
    int right =-1;
    int maxL = -1;
    int maxR = -1;
    for(int i = 0;i< n;i++)
    {
        if(A[i] > f+A[i])
        {
            left = i;
            right = i;
        }
        else
        {
            right = i;
        }
        f = max(f+A[i], A[i]);
        if(res < f)
        {
            maxL = left;
            maxR = right;
        }
        res = max(f, res);
    }
    cout<<"left "<<maxL<<endl;
    cout<<"right "<<maxR<<endl;
    cout<<"max "<< res<<endl;
    return res;
}


//13.3
void Detect(vector<vector<bool>>& tracker, int n, string s)
{
    for(int i = n-1;i>=0;i--)
    {
        for(int j = i;j <n;j++)
        {
            if(i==j)
            {
                tracker[i][j] = true;
            }
            else
            {
                if(s[i] == s[j] && (j-i ==1 || tracker[i+1][j-1]))
                {
                    tracker[i][j] = true;
                }
                else
                {
                    tracker[i][j] = false;
                }
                
            }
        }
    }
}

int minCut(string s) {
    int n = (int)s.length();
    vector<vector<bool> > tracker(n, vector<bool>(n, false));
    Detect(tracker,n,s);
    for(int i = 0;i<n;i++)
    {
        for(int j = 0;j<n;j++)
        {
            cout<< tracker[i][j] << " ";
        }
        cout<<endl;
    }
    int tmp[n];
    fill_n(&tmp[0], n, INT_MAX);
    
    for(int i =0;i<n;i++)
    {
        for(int j = i;j>=0;j--)
        {
            if(tracker[j][i])
            {
                tmp[i] = min(tmp[i], j == 0 ? 0: tmp[j-1]+1);
                cout<< "tmp["<<i<<"]= "<<tmp[i]<<endl;
            }
        }
    }
    return tmp[n-1];
    
}

//13.4
int maxHistogram(vector<int> input)
{
    input.push_back(0);
    int res = INT_MIN;
    stack<int> s;
    for(int i = 0;i<input.size();i++)
    {
        while(!s.empty() && input[s.top()] > input[i])
        {
            int height = input[s.top()];
            s.pop();
            int length = s.empty()? i: i- s.top()-1;
            res = max(res, height* length);
        }
        s.push(i);
    }
    return res;
}

int maximalRectangleII(vector<vector<char> > &matrix)
{
    vector<int> tracker(matrix[0].size(), 0);
    int res = 0;
    if(matrix.size() == 0 || matrix[0].size() == 0)
        return 0;
    for(int i = 0;i< matrix.size();i++)
    {
        for(int j = 0;j< matrix[i].size();j++)
        {
            if(matrix[i][j] == '0')
            {
                tracker[j] = 0;
            }
            else
            {
                tracker[j] ++;
            }
        }
        res = max(res, maxHistogram(tracker));
    }
    return res;
}

//13.5
int maxProfitIII(vector<int> &prices) {
    int n = (int)prices.size();
    if(n ==0) return 0;
    vector<int> f(n,0);
    vector<int> g(n,0);
    int profit = 0;
    for(int i = 1;i<n;i++)
    {
        profit = max(profit+prices[i]-prices[i-1], prices[i]-prices[i-1]);
        f[i] = max(f[i-1], profit);
    }
    
    profit = 0;
    for(int i = n-2;i>0;i--)
    {
        profit = max(profit+prices[i+1]-prices[i], prices[i+1]-prices[i]);
        g[i] = max(i == n-1? 0: g[i+1], profit);
    }
    int res = 0;
    for(int i = 1;i<n;i++)
    {
        res = max(res, f[i]+g[i]);
    }
    return res;
}

//13.6
template<typename InIt>
bool isInterleave(InIt first1, InIt last1, InIt first2, InIt last2, InIt first3, InIt last3)
{
    if(*first3 == *last3)
    {
        return first1 == last1 && first2 == last2;
    }
    else
    {
        bool success = false;
        
        if(*first1 == * first3)
        {
            success = success || isInterleave(next(first1), last1, first2, last2, next(first3), last3);
        }
        if(*first2 == *first3)
        {
            success = success || isInterleave(first1, last1, next(first2), last2, next(first3), last3);
        }
        return success;
    }
}

bool isInterleaveII(string s1, string s2, string s3)
{
    return isInterleave(s1.begin(), s1.end(), s2.begin(), s2.end(), s3.begin(), s3.end());
}

bool isInterleaveDpII(string s1, string s2, string s3)
{
    int m = (int)s1.size();
    int n = (int)s2.size();
    vector<vector<bool>> tracker (m+1, vector<bool>(n+1, false));
    tracker[0][0] = true;
    for(int j = 1;j<n+1;j++)
    {
        tracker[0][j] = tracker[0][j-1] && s2[j-1] == s3[j-1];
    }
    for(int i = 1;i<m+1;i++)
    {
        tracker[i][0] = tracker[i-1][0] && s1[i-1] == s3[i-1];
    }
    
    for(int i = 1;i<m+1;i++)
    {
        for(int j = 1;j<n+1;j++)
        {
            tracker[i][j] = (s1[i-1] == s3[i+j-1] && tracker[i-1][j]) ||
            (s2[j-1] == s3[i+j-1] && tracker[i][j-1]);
        }
    }
    return tracker[m][n];
        
}

//13.7
bool isScramble(string s1, string s2) {
    const int N = (int)s1.size();
    if (N != s2.size()) return false;
    bool f[N + 1][N][N];
    fill_n(&f[0][0][0], (N + 1) * N * N, false);
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    f[1][i][j] = s1[i] == s2[j];
    for (int n = 1; n <= N; ++n) {
        for (int i = 0; i + n <= N; ++i) {
            for (int j = 0; j + n <= N; ++j) {
                for (int k = 1; k < n; ++k) {
                    if ((f[k][i][j] && f[n - k][i + k][j + k]) ||
                        (f[k][i][j + n - k] && f[n - k][i + k][j])) {
                        f[n][i][j] = true;
                        break;
                    }
                }
            }
        }
    }
    return f[N][0][0];
}

//13.8
int minPathSum(vector<vector<int>> &grid) {
    int n = (int)grid.size();
    int m = (int)grid[0].size();
    vector<vector<int>> tracker(n, vector<int>(m, INT_MAX));
    for(int i =0;i<n;i++)
    {
        for(int j = 0;j<m;j++)
        {
            if(i ==0 && j ==0)
            {
                tracker[i][j] = grid[0][0];
            }
            else
            {
                tracker[i][j] = min( i== 0? INT_MAX : tracker[i-1][j] , j== 0? INT_MAX : tracker[i][j-1])
                + grid[i][j];
            }
        }
    }
    return tracker[n-1][m-1];
}

//13.9
int minDistanceII(string word1, string word2) {
    int m = (int)word1.size();
    int n = (int)word2.size();
    vector<vector<int>> tracker(m+1,vector<int>(n+1, 0));
    for(int i = 1; i<=n;i++)
    {
        tracker[0][i] = i;
    }
    
    for(int j = 1; j<=m;j++)
    {
        tracker[j][0] = j;
    }
    
    for(int i = 1;i<=m;i++)
    {
        for(int j = 1;j<=n;j++)
        {
            int tmp = min(tracker[i-1][j], tracker[i][j-1]);
            tmp = min(tmp+1, word1[i-1] == word2[j-1]? tracker[i-1][j-1]: tracker[i-1][j-1] +1);
            tracker[i][j] = tmp;
        }
    }
    return tracker[m][n];
}

//13.10
int numDecodings(string s) {
    if(s.length() == 0) return 0;
    if(s[0] == '0') return 0;
    if(s.length() == 1) return 1;
    if(s[0] > '2')
    {
        return 1 + numDecodings(s.substr(1, s.length()-1));
    }
    else
    {
        if(s[0] == '1' || (s[0] =='2' && s[1] <= '6'))
        {
            return 1+ numDecodings(s.substr(1, s.length()-1)) +
            numDecodings(s.substr(2, s.length()-2));
        }
        else
        {
            return 1 + numDecodings(s.substr(1, s.length()-1));
        }
    }
}

int numDecodingsDP(string s) {
    if(s.length() == 0) return 0;
    int prev = 0;
    int cur = 1;
    for(int i = 0;i<s.length();i++)
    {
        if(s[i] == '0')
            cur = 0;
        if(i< 1 || !(s[i-1] == '1' || (s[i-1] =='2' && s[i] <= '6')))
        {
            prev = 0;
        }
        int tmp = prev;
        prev = cur;
        cur = prev + tmp;

    }
    return cur;
}

//13.11
int numDistinctII(string S, string T) {
    int m = (int)S.size();
    int n = (int)T.size();

    if(m<n) return 0;
    int f[m];
    fill_n(&f[0], m, 0);
    if(S[0] == T[0]) f[0] = 1;
    for(int i = 1;i<m;i++)
    {
        if(S[i] == T[0])
            f[i] = f[i-1] + 1;
        else
            f[i] = f[i-1];
    }
    
    for(int i = 1;i<n;i++)
    {
        for(int j = m;j>i;j--)
        {
            if(T[i] == S[j])
            {
                f[j] = f[j-1] + f[j];
            }
        }
    }

    return f[m-1];
}

int numDistinctOld(string S, string T) {
    int n = (int)S.length();
    int m = (int)T.length();
    if(m>n) return 0;
    vector<vector<int>> tracker(m, vector<int>(n, 0));
    
    if(T[0] == S[0]) tracker[0][0] = 1;
    else tracker[0][0] = 0;
    
    for(int i = 1;i<n;i++)
    {
        if(T[0] == S[i])
        {
            tracker[0][i] = tracker[0][i-1]+1;
        }
        else
        {
            tracker[0][i] = tracker[0][i-1];
        }
    }
    
    for(int i = 1;i<m;i++)
    {
        for(int j = i;j<n;j++)
        {
            tracker[i][j] = tracker[i][j-1] + (T[i] == S[j] ? tracker[i-1][j-1] : 0);
        }
    }
    return tracker[m-1][n-1];
}

//13.12
bool wordBreak(string s, unordered_set<string> &dict) {
    int n = (int) s.size();
    bool f[n];
    fill_n(&f[0], n, false);
    for(int i = 0;i<n;i++)
    {
        for(int j= 0;j<=i;j++)
        {
            bool tmp = j == 0? true: f[j-1];
            string sub = s.substr(j, i-j+1);
            bool found = dict.find(sub) != dict.end();
            if(tmp && found)
            {
                f[i] = true;
                break;
            }
        }
    }
    return f[n-1];
}

//13.13
string constructStr(string s, vector<int> tracker)
{
    int n = (int)tracker.size();
    if(tracker[tracker.size()-1] != s.length()-1)
    {
        return "";
    }
    tracker.pop_back();
    for(int i = n-2;i>=0;i--)
    {
        s.insert(tracker[i]+1, " ");
    }
    cout<<s<<endl;
    return s;
}

void worker(string s, vector<int>& tracker, unordered_set<string>& dict, vector<string>& res)
{
    int n = (int)s.length();
    if(tracker.size() > 0 && tracker[tracker.size()-1] == n-1)
    {
        res.push_back(constructStr(s, tracker));
    }
    else
    {
        int k = tracker.size() == 0? -1 : tracker[tracker.size()-1];
        for(int i = k+1;i<n;i++)
        {
            if(dict.find(s.substr(k+1, i- (k+1)+1)) != dict.end())
            {
                tracker.push_back(i);
                worker(s, tracker, dict, res);
                tracker.pop_back();
            }
        }
    }
    
}

vector<string> wordBreakII(string s, unordered_set<string> &dict) {
    vector<string> res;
    vector<int> tracker;
    worker(s, tracker, dict, res);
    return res;
    
}

//13.13
void gen_path(const string &s, const vector<vector<bool> > &prev,
              int cur, vector<string> &path, vector<string> &result) {
    if (cur == 0) {
        string tmp;
        for (auto iter = path.crbegin(); iter != path.crend(); ++iter)
            tmp += *iter + " ";
        tmp.erase(tmp.end() - 1);
        result.push_back(tmp);
    }
    for (size_t i = 0; i < s.size(); ++i) {
        if (prev[cur][i]) {
            path.push_back(s.substr(i, cur - i));
            gen_path(s, prev, (int)i, path, result);
            path.pop_back();
        }
    }
}

vector<string> wordBreakIII(string s, unordered_set<string> &dict) {
    // 长度为 n 的字符串有 n+1 个隔板
    vector<bool> f(s.length() + 1, false);
    // prev[i][j] 为 true,表示 s[j, i) 是一个合法单词,可以从 j 处切开 // 第一行未用
    vector<vector<bool> > prev(s.length() + 1, vector<bool>(s.length()));
    f[0] = true;
    for (int i = 1; i <= (int)s.length(); ++i) {
        for (int j = i - 1; j >= 0; --j) {
            if (f[j] && dict.find(s.substr(j, i - j)) != dict.end()) {
                f[i] = true;
                prev[i][j] = true;
            }
        } }
    vector<string> result;
    vector<string> path;
    gen_path(s, prev, (int)s.length(), path, result);
    return result;
}

//15.1
int reverse(int x) {
    long long a = x;
    int sign = x>=0? 1: -1;
    long long res = 0;
    while(a>0)
    {
        res = res * 10 + a % 10;
        a = a/10;
    }
    return (int)(sign*res);
    
}

//15.2
bool isPalindrome(int a) {
    
    long long x = abs((long long)a);
    int i = 0;
    while(pow(10, i)<= x)
    {
        i++;
    }
    i --;
    while(x>0)
    {
        int left = x/pow(10,i);
        int right = x%10;
        if(left != right) return false;
        x = x- right;
        x = x- left * pow(10,i);
        x = x/10;
        i = i-2;
    }
    return true;
}

//15.3
vector<Interval> insertII(vector<Interval> &intervals, Interval newInterval) {
    auto it = intervals.begin();
    while(it != intervals.end())
    {
        if(newInterval.end < it->start)
        {
            intervals.insert(it, newInterval);
            return intervals;
        }
        else if(it->end < newInterval.start)
        {
            it ++;
        }
        else
        {
            newInterval.start = min(it->start, newInterval.start);
            newInterval.end = max(it->end, newInterval.end);
            it = intervals.erase(it);
        }
    }
    intervals.push_back(newInterval);
    return intervals;
}

vector<Interval> insertDCII(vector<Interval> &intervals, Interval newInterval)
{
    if(intervals.size() ==0)
    {
        intervals.push_back(newInterval);
        return intervals;
    }
    int start = -1;
    int left = 0;
    int right = (int)intervals.size()-1;
    while(left<=right)
    {
        int mid = left + (right-left)/2;
        if(newInterval.start < intervals[mid].start)
        {
            right = mid - 1;
        }
        else if(newInterval.start > intervals[mid].end)
        {
            left = mid + 1;
        }
        else
        {
            newInterval.start = min(newInterval.start, intervals[mid].start);
            newInterval.end = max(newInterval.end, intervals[mid].end);
            start = mid;
            break;
        }
    }
    start = start==-1? right+1:start;

    int end = -1;
    left = 0;
    right = (int)intervals.size() -1;
    while(left<=right)
    {
        int mid = left +(right-left)/2;
        if(newInterval.end < intervals[mid].start)
        {
            right = mid -1;
        }
        else if(newInterval.end > intervals[mid].end)
        {
            left = mid +1;
        }
        else
        {
            newInterval.start = min(newInterval.start, intervals[mid].start);
            newInterval.end = max(newInterval.end, intervals[mid].end);
            end = mid;
            break;
        }
    }
    
    end = end == -1? left-1: end;
    if(start > end)
    {
        intervals.insert(intervals.begin() + start, newInterval);
    }
    else
    {
        intervals.erase(intervals.begin()+ start, intervals.begin()+ end +1);
        intervals.insert(intervals.begin() + start, newInterval);
    }
    
    return intervals;
}

vector<Interval> insertDC(vector<Interval> &intervals, Interval newInterval) {
    int first = -1;
    int left = 0;
    int right = (int)intervals.size()-1;
    while(left<=right)
    {
        int mid = left + (right-left)/2;
        if(newInterval.start <=intervals[mid].end && newInterval.start >= intervals[mid].end)
        {
            newInterval.start = min(intervals[mid].start, newInterval.start);
            newInterval.end = max(intervals[mid].end, newInterval.end);
            first = mid;
            break;
        }
        else if(newInterval.start < intervals[mid].start)
        {
            right = mid-1;
        }
        else
        {
            left = mid+1;
        }
    }
    if(first == -1) first = left == intervals.size()? left -1: left;
    int last = -1;
    left = 0;
    right = (int)intervals.size()-1;
    while(left<=right)
    {
        int mid = left + (right-left)/2;
        if(newInterval.end >= intervals[mid].start && newInterval.end <= intervals[mid].end)
        {
            newInterval.start = min(intervals[mid].start, newInterval.start);
            newInterval.end = max(intervals[mid].end, newInterval.end);
            last = mid;
            break;
        }
        else if(newInterval.end < intervals[mid].start)
        {
            right = mid -1;
        }
        else
        {
            left = mid+1;
        }
    }
    if(last == -1) last = right <0 ? 0: right;
    intervals.erase(intervals.begin()+ left, intervals.begin()+ right);
    intervals.insert(intervals.begin() + left, newInterval);
    return intervals;
}

//15.4

bool sortIntervalFunctionsII(Interval i1, Interval i2){
    return (i1.start < i2.start);
}

vector<Interval> mergeII(vector<Interval> &intervals) {
    sort(intervals.begin(), intervals.end(), sortIntervalFunctionsII);
    vector<Interval> res;
    Interval cur;
    if(intervals.size() ==0)
    {
        return res;
    }
    else
    {
        cur = intervals[0];
    }
    for(int i = 1;i<intervals.size();i++)
    {
        if(cur.end < intervals[i].start)
        {
            res.push_back(cur);
            cur = intervals[i];
        }
        else
        {
            cur.end = max(cur.end, intervals[i].end);
        }
    }
    res.push_back(cur);
    return res;
}


//15.5
string minWindow(string S, string T) {
    int TAlphaStat[260];
    fill_n(TAlphaStat, 260, 0);
    for (char& c : T)
    {
        TAlphaStat[c]++;
    }
    int currentAlphaStat[260];
    fill_n(currentAlphaStat, 260, 0);
    int head = 0, tail = 0;
    int containNum = 0;
    int minWidth = INT_MAX, left = -1;
    while (tail < S.length())
    {
        char& tc = S[tail];
        if (TAlphaStat[tc] > 0){
            if (currentAlphaStat[tc] < TAlphaStat[tc]){
                containNum++;
            }
            currentAlphaStat[tc]++;
        }
        if (containNum == T.length()){
            while (head < tail && (currentAlphaStat[S[head]] > TAlphaStat[S[head]]
                                   || TAlphaStat[S[head]] == 0))
            {
                currentAlphaStat[S[head]]--;
                head++;
            }
            if (minWidth > tail - head + 1){
                minWidth = tail - head + 1;
                left = head;
            }
        }
        tail++;
    }
    return left == -1 ? "" : S.substr(left, minWidth);
}

//15.6
string multiplyII(string num1, string num2) {
    if(num1 == "0" || num2 == "0") return "0";
    vector<int> n1;
    vector<int> n2;
    for(int i = (int)num1.size()-1;i>=0;i--)
    {
        n1.push_back(num1[i]-'0');
    }
    for(int i = (int)num2.size()-1;i>=0;i--)
    {
        n2.push_back(num2[i]-'0');
    }
    vector<int> res(n1.size() + n2.size() + 1, 0);
    int carry = 0;
    for(int i = 0;i<n1.size();i++)
    {
        carry = 0;
        for(int j = 0;j<n2.size();j++)
        {
            int tmp = carry + n1[i] * n2[j] + res[i+j];
            carry = tmp /10;
            res[i+j] = tmp %10;
        }
        res[i+n2.size()] += carry;
    }
    while (res[res.size()-1] ==0) {
        res.pop_back();
    }
    reverse(res.begin(), res.end());
    string s;
    for(int i = 0;i<res.size();i++)
    {
        s.append(1, '0' + res[i]);
    }
    return s;
}

//15.7
vector<int> findSubstring(string S, vector<string> &L) {
    vector<int> res;
    int length = (int)L[0].size();
    int totalwords = (int)L.size();
    if(S.size() % length != 0) return res;
    if(S.size() < length * L.size()) return res;
    
    unordered_map<string, int> dict;
    
    for(int i = 0;i < L.size() ;i++)
    {
        if(dict.find(L[i]) == dict.end())
        {
            dict[L[i]] = 1;
        }
        else
        {
            dict[L[i]] +=1;
        }
    }
    
    unordered_map<string, int> tracker;
    int startIndex = -1;
    int curCount= 0;
    for(int j = 0;j<= S.size() - length*totalwords;j++)
    {
        for(int i = j;i<= S.size() -length;i += length)
        {
            if(startIndex <0)
            {
                startIndex = i;
            }
            if(dict.find(S.substr(i, length))!= dict.end())
            {
                curCount++;
                if(tracker.find(S.substr(i, length)) == tracker.end())
                {
                    tracker[S.substr(i, length)] = 1;
                }
                else
                {
                    tracker[S.substr(i, length)]++;
                }
                if(tracker[S.substr(i, length)] > dict[S.substr(i, length)])
                {
                    tracker.clear();
                    tracker[S.substr(i, length)] = 1;
                    curCount = 1;
                    startIndex = i;
                }else if (curCount == totalwords) {
                    res.push_back(startIndex);
                    tracker.clear();
                    curCount = 0;
                    startIndex = -1;
                }
            }
            else if (dict.find(S.substr(i, length))== dict.end())
            {
                startIndex = -1;
                curCount = 0;
                tracker.clear();
            }
        }
    }
    return res;
}

vector<int> findSubstringIII(string s, vector<string>& dict) {
    size_t wordLength = dict.front().length();
    size_t catLength = wordLength * dict.size();
    vector<int> result;
    if (s.length() < catLength) return result;
    unordered_map<string, int> wordCount;
    
    for (auto const& word : dict) ++wordCount[word];
    
    for (auto i = begin(s); i <= prev(end(s), catLength); ++i)
    {
        unordered_map<string, int> unused(wordCount);
        for (auto j = i; j != next(i, catLength); j += wordLength)
        {
            auto pos = unused.find(string(j, next(j, wordLength)));
            if (pos == unused.end() || pos->second == 0)
                break;
            if (--pos->second == 0)
                unused.erase(pos);
        }
        if (unused.size() == 0)
            result.push_back((int)distance(begin(s), i));
    }
    return result;
}

//15.8
vector<vector<int> > generate(int numRows) {
    vector<vector<int>> res;
    if(numRows == 0) return res;
    vector<int> tmp = {1};
    res.push_back(tmp);
    for(int i = 1;i<numRows;i++)
    {
        vector<int> cur;
        for(int j = 0;j<=tmp.size();j++)
        {
            int n = ((j-1) >= 0? tmp[j-1]: 0) + (j<tmp.size()? tmp[j]:0);
            cur.push_back(n);
        }
        res.push_back(cur);
        swap(tmp, cur);
        cur.clear();
    }
    return res;
}


//15.9
vector<int> getRow(int rowIndex) {
    vector<int> tmp = {1};
    for(int i = 1;i<=rowIndex;i++)
    {
        vector<int> cur;
        for(int j = 0;j<=tmp.size();j++)
        {
            int n = ((j-1) >= 0? tmp[j-1]: 0) + (j<tmp.size()? tmp[j]:0);
            cur.push_back(n);
        }
        swap(tmp, cur);
        cur.clear();
    }
    return tmp;
}

//15.10
vector<int> spiralOrderII(vector<vector<int> > &matrix) {
    vector<int> res;
    int m = (int)matrix.size();
    if(m == 0) return res;
    int n = (int)matrix[0].size();
    if(n == 0) return res;
    int top = 0;
    int right = n-1;
    int bottom = m-1;
    int left = 0;
    
    int direction = 0;
    while(top<=bottom && left <= right)
    {
        if(direction == 0)
        {
            for(int j = left; j<=right;j++)
            {
                res.push_back(matrix[top][j]);
            }
            top ++;
        }
        else if(direction == 1)
        {
            for(int i = top;i<=bottom;i++)
            {
                res.push_back(matrix[i][right]);
            }
            right--;
        }
        else if(direction == 2)
        {
            for(int j = right;j>=left;j--)
            {
                res.push_back(matrix[bottom][j]);
            }
            bottom --;
        }
        else if(direction == 3)
        {
            for(int i = bottom;i>=top;i--)
            {
                res.push_back(matrix[i][left]);
            }
            left ++;
        }
        direction = (direction + 1) % 4;
    }
    return res;
}

//15.11
vector<vector<int> > generateMatrix(int n) {
    vector<vector<int>> res (n, vector<int>(n,0));
    if(n ==0) return res;
    int top = 0;
    int right = n-1;
    int bottom = n-1;
    int left = 0;
    int direction = 0;
    int tmp = 1;
    while(top <= bottom && left <= right)
    {
        if(direction == 0)
        {
            for(int j = left; j<=right;j++)
            {
                res[top][j] = tmp++;
            }
            top ++;
        }
        else if(direction == 1)
        {
            for(int i = top;i<=bottom;i++)
            {
                res[i][right] = tmp++;
            }
            right--;
        }
        else if(direction == 2)
        {
            for(int j = right;j>=left;j--)
            {
                res[bottom][j] = tmp ++;
            }
            bottom --;
        }
        else if(direction == 3)
        {
            for(int i = bottom;i>=top;i--)
            {
                res[i][left] = tmp++;
            }
            left ++;
        }
        direction = (direction + 1) % 4;
        
    }
    return res;
}

//15.12
string worker(int L, int curLen, vector<string> tracker, bool lastline)
{
    int spaces = L - curLen;
    int averspace = 0;
    int extra;
    if( tracker.size() > 1)
    {
        averspace = spaces / (tracker.size()-1);
        extra = spaces - averspace * ((int)tracker.size() -1);
    }
    else
    {
        averspace = spaces;
        extra = 0;
    }
    string s = "";
    if(lastline)
    {
        for(int i= 0; i< tracker.size();i++)
        {
            s += tracker[i];
            if( i != tracker.size() -1)
            {
                s += " ";
            }
        }
        
        while(s.length() < L)
        {
            s += " ";
        }
        
    }
    else
    {
        if(tracker.size() == 1)
        {
            s += tracker[0];
            for(int j = 0;j< averspace;j++)
            {
                s += " ";
            }
        }
        else
        {
            for(int j = 0;j< tracker.size();j++)
            {
                s += tracker[j];
                if(j < (tracker.size() - 1))
                {
                    for(int j = 0;j< averspace;j++)
                    {
                        s += " ";
                    }
                    if(extra >0)
                    {
                        s += " ";
                        extra --;
                    }
                }
            }
        }
    }
    return s;
}

vector<string> fullJustifyII(vector<string> &words, int L) {
    vector<string> res;
    if(words.size() == 0 ) return res;
    
    vector<string> tracker;
    int curLen = 0;
    for(int i = 0;i< words.size();i++)
    {
        if(( curLen + words[i].length() + (int)tracker.size() > L) && tracker.size() != 0)
        {
            string s = worker(L, curLen, tracker, false);
            res.push_back(s);
            curLen = 0;
            tracker.clear();
        }
        
        curLen = curLen + (int)words[i].length();
        tracker.push_back(words[i]);
        
    }
    if(tracker.size() >0)
    {
        string s = worker(L, curLen, tracker, true);
        res.push_back(s);
        curLen = 0;
        tracker.clear();
    }
    return res;
}


//15.13
int divideIII(int dividend, int divisor) {
    // 当 dividend = INT_MIN 时,-dividend 会溢出,所以用 long long
    long long a = dividend >= 0 ? dividend : -(long long)dividend;
    long long b = divisor >= 0 ? divisor : -(long long)divisor;
    // 当 dividend = INT_MIN 时,divisor = -1 时,结果会溢出,所以用 long long
    long long result = 0;
    while (a >= b) {
        long long c = b;
        for (int i = 0; a >= c; ++i, c <<= 1) {
            a -= c;
            result += 1 << i;
        }
    }
    return ((dividend^divisor) >> 31) ? (int)(-result) : (int)(result);
}

//15.15
int maxPoints(vector<Point> &points) {
    unordered_map<float, int> tracker;
    if(points.size() == 0) return 0;
    int res = 1;
    int samepoint = 0;
    for(int i = 0;i< points.size();i++)
    {
        for(int j = 0;j< points.size();j++)
        {
            if(i != j)
            {
                if(points[i].x == points[j].x && points[i].y == points[j].y)
                {
                    samepoint++;
                }
                else
                {
                    if(points[i].x == points[j].x)
                    {
                        tracker[INT_MAX] ++;
                        res = max(res, tracker[INT_MAX]+1 + samepoint);
                    }
                }
                
                float sl =(float)(points[i].y - points[j].y) / (float)(points[i].x - points[j].x);
                tracker[sl] ++;
                res = max(res, tracker[sl]+1 + samepoint);
            }
        }
        tracker.clear();
        samepoint = 0;
    }
    return res;
}

//16.1
void reverse_wordWorker(string& s, int left, int right)
{
    while(left <= right)
    {
        char tmp = s[left];
        s[left] = s[right];
        s[right] = tmp;
        left ++;
        right --;
    }
    
}


void reverseWords(string &s) {
    int startIndex = -1;
    for(int i = 0;i<s.length();i++)
    {
        if(s[i] != ' ')
        {
            startIndex = i;
            break;
        }
    }
    if(startIndex == -1)
    {
        s = "";
    }
    else
    {
        s = s.substr(startIndex, s.length() - startIndex);
    }
    
    reverse_wordWorker(s, 0, (int)s.length()-1);
    int left = -1;
    int right = -1;
    for(int i = 0;i<s.length();i++)
    {
        if((i == 0 || s[i-1] == ' ') && s[i] != ' ')
        {
            left = i;
        }
        else if((i == s.length() -1 || s[i + 1] == ' ') && s[i] != ' ')
        {
            right = i;
        }
        
        if(left != -1 && right != -1)
        {
            reverse_wordWorker(s, left, right);
        }
    }
    
    startIndex = -1;
    for(int i = 0;i<s.length();i++)
    {
        if(s[i] != ' ')
        {
            startIndex = i;
            break;
        }
    }
    if(startIndex == -1)
    {
        s = "";
    }
    else
    {
        s = s.substr(startIndex, s.length() - startIndex);
    }
    
    
}

//16.2
int findMin(vector<int>& num)
{
    int size = (int)num.size();
    if(size == 0) return 0;
    int left = 0;
    int right =  size-1;
    while(left <=right)
    {
        int mid = left + (right-left)/2;
        
        if((num[mid] < num[mid-1] || mid == 0) &&
           (num[mid] < num[mid+1] || mid == size -1))
        {
            return num[mid];
        }
        else if(num[left] < num[mid] && num[mid] < num[right])
        {
            return num[left];
        }
        else if(num[left] > num[mid])
        {
            right = mid - 1;
        }
        else if(num[mid] > num[right])
        {
            left = mid + 1;
        }
    }
    return num[left-1];
}

// serialize and deserialize tree
vector<int> res;
void serializeTree(TreeNode* root)
{
    if(!root)
    {
        res.push_back(-1);
        return;
    }
    res.push_back(root->val);
    serializeTree(root->left);
    serializeTree(root->right);
}

int i;
TreeNode* deserialize()
{
    if(i >= res.size()) return NULL;
    if(res[i] < 0)
    {
        i++;
        return NULL;
    }
    TreeNode* root = new TreeNode(res[i++]);
    root->left = deserialize();
    root->right = deserialize();
    return root;
        
}

// serialize and deserialize n-ary tree
struct NaryNode{
    char val;
    vector<NaryNode*> child;
    NaryNode(int v): val(v){};
};

// this method use a global variable res to take the out put.
void serializeNary(NaryNode* root)
{
    if(!root) { return;}
    res.push_back(root->val);
    for(auto item: root->child)
    {
        serializeNary(item);
    }
    res.push_back(-1);
}

NaryNode* deserializeNary()
{
    while(i<(int)res.size())
    {
        if(res[i] <0) {i++; return NULL;}
        NaryNode* tmp = new NaryNode(res[i]);
        i++;
        NaryNode* child = deserializeNary();
        while(child)
        {
            tmp->child.push_back(child);
            child = deserializeNary();
        }
        return tmp;
    }
    return NULL;
}

// print all factors
void findfactor(int target,int start, vector<vector<int>>& result, vector<int> &group)
{
    if(target == 1)
    {
        if(group.size()==1)
        {
            group.insert(group.begin(),1);
        }
        result.push_back(group);
        return;
    }
    for(int i= start;i<=target;i++)
    {
        if(target%i==0)
        {
            group.push_back(i);
            findfactor(target/i, i, result, group);
            group.pop_back();
        }
    }
}

vector<vector<int>> printAllFactor(int n)
{
    vector<vector<int>> res;
    vector<int> tracker;
    findfactor(n,2, res, tracker);
    return res;
}

/*
 一个数组，保证前半部递增，后半部递减，求数组最大值。二分查找，没写出来。大家有近期面的，练一下。
 */
int findMax(vector<int> input)
{
    if(input.size() <1) return -1;
    if(input.size() == 1) return input[0];
    int left = 0;
    int right = (int)input.size()-1;
    while(left<=right)
    {
        int mid = left + (right-left)/2;
        if(mid == 0 and input[0] > input[1])
        {
            return input[0];
        }
        else if(mid == input.size()-1 && input[mid] > input[mid-1])
        {
            return input[mid];
        }
        else if (input[mid] > input[mid-1] && input[mid] > input[mid+1])
        {
            return input[mid];
        }
        else if(input[mid] < input[mid+1])
        {
            left = mid+1;
        }
        else
        {
            right = mid -1;
        }
    }
    return -1;
}

/*
 输入是一个自然数T， 输出是(a_1,a_2,...,a_k)使得a_1^2+a_2^2+...+a_k^2=T， 并且k尽可能小
 */
bool isSqaure(int n)
{
    return (sqrt(n) - (int)sqrt(n)) == 0.0;
}
vector<int> findSumOfPerfectSqaure(int t)
{
    vector<vector<int>> res(t+1);
    res[1].push_back(1);
    int minCount = INT_MAX;
    vector<int> tmp;
    for(int i = 2;i<=t;i++)
    {
        minCount = INT_MAX;
        tmp.clear();
        for(int j = 0;j<i;j++)
        {
            if(isSqaure(i-j))
            {
                if(minCount > (res[j].size() + 1))
                {
                    minCount = (int)res[j].size() + 1;
                    tmp = res[j];
                    tmp.push_back(sqrt(i-j));
                }
            }
        }
        res[i] = tmp;
    }
    return res[t];
}

/*
 1.1. Tokenize a string to words. Ignore any space and punctuator
 */
bool isDelimitor(char c)
{
    if(c == ' ') return true;
    return false;
}
vector<string> tokenize(string s)
{
    vector<string> res;
    if(s.length() ==0) return res;
    string tmp = "";
    for(int i = 0;i< s.length();i++)
    {
        if(!isDelimitor(s[i]))
        {
            tmp = tmp + s[i];
        }
        else if(tmp != "")
        {
            res.push_back(tmp);
            tmp = "";
        }
    }
    if(tmp != "")
    {
        res.push_back(tmp);
    }
    return res;
}

/*
 输入是一个 N*N的矩阵，代表地势高度（elevation）。然后如果下雨，水流只能流去比他矮或者一样高的地势。矩阵上边和左边是太平洋，下边和右边是大西洋。求出所有的能同时流到两个大洋的点
 */
struct Elevation{
    int val;
    bool visited;
    Elevation(int n): val(n), visited(false) {}
};
bool up = false;
bool down = false;
bool upset = false;
bool downset = false;

void worker(vector<vector<Elevation>>& input, int i, int j)
{
    if(upset && downset)
    {
        return;
    }
    else
    {
        if(!input[i][j].visited)
        {
            input[i][j].visited = true;
            if (i == 0 || j == 0)
            {
                up = true;
                upset = true;
            }
            if( i == input.size()-1 || j == input.size()-1)
            {
                down = true;
                downset= true;
            }

            if(i - 1 >=0 && input[i][j].val > input[i-1][j].val)
            {
                worker(input, i-1, j);
            }
            if(i + 1 < (int)input.size() && input[i][j].val > input[i+1][j].val)
            {
                worker(input, i+1, j);
            }
            if(j - 1 >=0 && input[i][j].val > input[i][j-1].val)
            {
                worker(input, i, j-1);
            }
            if(j + 1 < (int)input.size() && input[i][j].val > input[i][j+1].val)
            {
                worker(input, i, j+1);
            }
            
        }
    }
}

void cleanVisited(vector<vector<Elevation>>& input)
{
    int n = (int)input.size();
    for(int i = 0;i< n;i++)
        for(int j = 0;j<n;j++)
            input[i][j].visited = false;

}

vector<Elevation> canFlowinto(vector<vector<Elevation>>& input)
{
    vector<Elevation> res;
    int n = (int)input.size();
    for(int i = 0;i< n;i++)
    {
        for(int j = 0;j<n;j++)
        {
            cleanVisited(input);
            up = false;
            down = false;
            upset = false;
            downset = false;
            
            worker(input, i, j);
            if(up && down)
            {
                res.push_back(input[i][j]);
            }
        }
    }
    return res;
}

/*
 求一个unsorted数组的前k大数字，要求O(n)，这题被烙印坑了。给了个O(n)算法非
 说我是O(nlogn)，最后说服了他，不过时间不够了
 */
void findFirstk(vector<int>& input, int start, int end, int k)
{
    if(start>=end) return;
    int index = rand() % (end- start);
    swap(input, start, start+index);
    int runner = start + 1;
    int right = end;
    while(runner<=right)
    {
        if(input[runner]<= input[start]) runner++;
        else swap(input, runner, right--);
    }

    swap(input, start, runner-1);
    if(runner - start == k) return;
    else if(runner-start >k) findFirstk(input, start, runner-1, k);
    else findFirstk(input, runner, end, k-runner+start);
}

/*
 二维Matrix字符里面找word，所有字符必须相邻，比如microsoft 是否出现，相邻字符必须是neighbor关系，貌似career up上别人post 过。
 */
/*
bool worker(vector<vector<char>>& input, unordered_set<pair<int,int>>& visited, int i, int j, string& tracker, string target)
{
    if(tracker == target)
    {
        return true;
    }
    else
    {
        if(visited.find(make_pair(i, j)) == visited.end())
        {
            visited.insert(make_pair(i, j));
            tracker = tracker + input[i][j];
            if(i-1>=0)
            {
                if(worker(input, visited, i-1, j, tracker, target))
                {
                    return true;
                }
            }
            if(i+1 < input.size())
            {
                if(worker(input, visited, i+1, j, tracker, target))
                {
                    return true;
                }
            }
            if(j-1>0)
            {
                if(worker(input, visited, i, j-1, tracker, target))
                {
                    return true;
                }
            }
            if(j+1<input.size())
            {
                if(worker(input, visited, i, j+1, tracker, target))
                {
                    return true;
                }
            }
        }
        tracker.erase(tracker.size()-1);
    }
    return false;
}


bool findWord(vector<vector<char>>& input, string target)
{
    unordered_set<pair<int,int>> visited;
    for(int i = 0;i< input.size();i++)
    {
        for(int j = 0;j<input.size();j++)
        {
            visited.clear();
            string tracker= "";
            if(worker(input, visited, i, j, tracker, target))
            {
                return true;
            }
        }
    }
    return false;
}
*/

void quickSort(int input[], int n)
{
    if(n <=1) return;
    int index = rand() % n;
    swap(input, 0, index);
    int left = 1;
    int right = n-1;
    while(left<=right)
    {
        if(input[left] <= input[0]) left++;
        else swap(input, left, right--);
    }

    swap(input, 0, left-1);
    quickSort(input, left-1);
    quickSort(input + left, n - left);
}

/*
 一个class {int a,  bool c, int b} 里面每个variable所占的空间都不同
 ，比如a,b是int 所以分别占4byte. bool的c只占1byte。还有其他变量，可能占8bytes
 或者16bytes。都是2的次方就是。
 问题是写一个程序让他们可以很好的被放到8byte为单位的block里面去然后空间不会浪
 费。
 比如如果是 就按照a, c, b的话它一共要占12个byte。因为当把a和c放到一个block的
 时候就会浪费一些空间。
 所以最好摆成a，b，c这样的话更合理。占9个byte。剩下的空间还可以放一些小的
 object。
 其实这个就是用排序，然后从大的变量依次放进block。
 有个followup的问题就是：因为我不想过多移动这些变量，所以怎么才能设计一个算法
 所需要移动的object最少。
 比如如果变量的size一次是4, 4, 1, 1, 8, 8, 1, 1最好的排法是4, 4, 8, 8, 1, 1,
 1, 1.而不是8 8 4 4 1 1 1 1因为前一种所需要移动的cost最小。这个没想出来了。。
 应该用divide and conquer？
 */
void quickArrange(int input[], int n)
{
    int left = 0;
    int right = n-1;
    while(left<=right)
    {
        if(input[left] >=4) left++;
        else swap(input, left, right--);
    }
}
/*
 give an array [1, 3, 4, 6, 9] and 10 return [2,5,7-8,10]
 */
vector<string> complementry(vector<int> input, int n)
{
    vector<string> res;
    int startnumber = 1;
    for(int i = 0;i < (int)input.size();i++)
    {
        if(input[i] == startnumber)
        {
            startnumber = input[i] +1;
        }
        else if(input[i] > startnumber && input[i] < n)
        {
            string tmp = "";
            tmp += to_string(startnumber);
            if(input[i] -1 > startnumber)
            {
                tmp += "-";
                tmp += to_string(input[i] - 1);
            }
            res.push_back(tmp);
            startnumber = input[i] +1;
        }
        else
        {
            break;
        }
    }
    
    if(startnumber <= n)
    {
        string tmp = "";
        tmp += to_string(startnumber);
        if(n > startnumber)
        {
            tmp += "-";
            tmp += to_string(n);
        }
        res.push_back(tmp);
    }
    return res;
}
/*
去除string中的空白
*/
string removeSpace(string s)
{
    int index = 0;
    for(int i = 0;i< s.length();i++)
    {
        if(s[i] != ' ')
        {
            s[index++] = s[i];
        }
        if(s[i] == '\0')
        {
            break;
        }
    }
    return s;
}

/*
 把regular expression tree 转换成表达式string最后一题没写完就到只剩下五分钟了。小哥让我停下来跟我介绍了他的组，并且问问我有什么问题之类。然后就结束了。题目都没答完肯 定是没戏了。move on准备下一场。祝各位找工作的都顺利！
 */
struct ExpressionTreeNode
{
    string val;
    ExpressionTreeNode* left;
    ExpressionTreeNode* right;
    ExpressionTreeNode(string input) : val(input), left(NULL), right(NULL) {}
};

string convert(ExpressionTreeNode * root)
{
    string res ="";
    if(!root) return res;
    if(root->val != "*" &&
       root->val != "/" &&
       root->val != "+" &&
       root->val != "-")
    {
        return root->val;
    }
    else
    {
        string res = "";
        if(root->val == "+" ||
           root->val == "-")
        {
            res += "(";
        }
        // assume this is a perfect format expression tree.
        // if a node is *,/,+,- it must has left and right
        // child as numbers.
        res += convert(root->left);
        res += root->val;
        res += convert(root->right);
        if(root->val == "+" ||
           root->val == "-")
        {
            res += ")";
        }
        
    }
    return res;
}

/*
 有很多 字串， 经常要作的操作有 ， 插入， 删除， 清空，查询以某前缀字串开头
 的所有的字串查询以某前缀字串开头的所有的字串的个数查询前缀字串开头的字串的所有
 可能的下一个字母例如
 
 [abc, abd, abe]
 input ab
 
 return [c, d, e]
 
 要求使3个查寻操作的时间上最优, 插入, 删除, 清空, 的性能表现可以牺牲
 */
void worker(vector<string>& res, string& tmp, TrieNode* root, bool isEnd)
{
    if(isEnd)
    {
        res.push_back(tmp);
    }
    
    for(int i = 0;i< 256;i++)
    {
        if(root->val[i])
        {
            tmp += (char)i;
            worker(res, tmp, root->val[i], root->end[i]);
            tmp.pop_back();
        }
    }
}

vector<string> search(string preFix, TrieNode* root)
{
    TrieNode* tmp = root;
    for(int i = 0;i<preFix.length();i++)
    {
        tmp = tmp->val[preFix[i]];
    }
    vector<string> res;
    string s = "";
    worker(res, s, tmp, false);
    return res;
}

/*
 要顺序打印 power(2,x) * power(3,y) * power(5,z).  x, y, z >= 0.
 print {1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16...}
 */
vector<int> printMultiples(int n)
{
    vector<int> x = {2,3,5};
    vector<int> res;
    res.push_back(1);
    for(int i = 1;i< n;i++)
    {
        int minVal = INT_MAX;
        for(int j = 0;j<res.size();j++)
        {
            for(int k = 0;k<x.size();k++)
            {
                if(x[k] * res[j] > res[res.size()-1])
                {
                    minVal = min (minVal, x[k] * res[j]);
                }
            }
        }
        res.push_back(minVal);
    }
    return res;
}

/*
12. binary search tree deletion node
*/
TreeNode* deleteNote(TreeNode*& root, TreeNode* node)
{
    if(!root) return root;
    if(root == node)
    {
        if(!root->left && !root->right) return NULL;
        if(!root->left && root->right) return root->right;
        if(root->left && !root->right) return root->left;
        if(root->left && root->right)
        {
            TreeNode* left = root->left;
            TreeNode* right = root->right;
            while(left->right) left = left->right;
            left->right = right;
            return left;
        }
    }
    else if(node->val < root->val)
    {
        root->left = deleteNote(root->left, node);
    }
    else if(node->val > root->val)
    {
        root->right = deleteNote(root->right, node);
    }
    return root;
}

/*
 第一个是两个单词最短距离，在版上看到很多人都说过这个题目，应该是L家经常面的。
 */
int findMinDis(vector<string> input, string s1, string s2)
{
    unordered_map<string, vector<int>> tracker;
    for(int i = 0;i< (int)input.size();i++)
    {
        tracker[input[i]].push_back(i);
    }
    
    vector<int> first;
    vector<int> second;
    if(tracker.find(s1) != tracker.end())
    {
        first = tracker[s1];
    }
    else
    {
        return INT_MAX;
    }
    
    if(tracker.find(s2) != tracker.end())
    {
        second = tracker[s2];
    }
    else
    {
        return INT_MAX;
    }
    
    int res = INT_MAX;
    int i = 0;
    int j = 0;
    while(i< (int) first.size() && j < second.size())
    {
        res = min(res, abs(first[i] - second[j]));
        if(first[i] < second[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    return res;
}

bool isMatch4(const char *s, const char *p) {
    if( *s == '\0' && *p == '\0') return true;
    
    if(*(p+1) == '*')
    {
        if(*p == *s || (*p == '.' && *s != '\0'))
        {
            return  isMatch4(s+1, p) || isMatch4(s, p+2);
        }
        else
        {
            return isMatch4(s, p+2);
        }
       
    }
    else if(*p == *s || (*p == '.' && *s != '\0'))
    {
        return isMatch4(s+1, p+1);
    }
    
    return false;
}


double sqrt (double value, double tolerance)
{
    double left;
    double right;
    
    if(value >1)
    {
        left = 0;
        right = value;
    }
    else
    {
        left = value;
        right = 1;
    }
    
    while(abs(right-left) > tolerance)
    {
        double mid = (left+right)/2.0;
        if(mid* mid > value)
        {
            right = mid;
        }
        else
        {
            left = mid;
        }
    }
    return left;
}

/*
 9. sort strings like "TADTTTTBDB", fixed length of 10, made up by only four characters: T A D B want linear time.
 */
void swap(string& s, int i , int j)
{
    char c = s[i];
    s[i] = s[j];
    s[j] = c;
}
void sort(string& s, int left, int right)
{

    while(left<=right)
    {
        if(s[left] == 'D')
        {
            left++;
        }
        else
        {
            swap(s, left, right--);
        }
    }
}
string sort(string& s)
{
    int n = (int)s.length();
    int left = 0;
    int mid = 0;
    int right = n-1;
    while(mid <=right)
    {
        if(s[mid] == 'A')
        {
            swap(s, left++, mid++);
        }
        else if(s[mid] == 'B')
        {
            mid++;
        }
        else if(s[mid] == 'T' || s[mid] == 'D')
        {
            swap(s, mid, right--);
        }
    }
    sort(s, mid, n-1);
    return s;
}

/*
 20: largest sum of value at certain time point, with given: [t0, t1, value1], [t2, t3, value2], [t3, t4, value3]…
 */
struct workItem{
    int start;
    int end;
    int val;
    workItem(int s, int e, int v): start(s), end(e), val(v) {}
};

struct pointEntry{
    int endPoint;
    int val;
    pointEntry(int e, int v): endPoint(e), val(v) {}
};

bool sortPointEntry(pointEntry* p1, pointEntry* p2)
{
    if(p1->endPoint != p2->endPoint)
    {
        return p1->endPoint < p2->endPoint;
    }
    else
    {
        return p1->val < p2->val;
    }
}

int findMax(vector<workItem> input)
{
    vector<pointEntry*> tracker;
    for(auto item: input)
    {
        pointEntry* start = new pointEntry(item.start, item.val);
        pointEntry* end = new pointEntry(item.end, -item.val);
        tracker.push_back(start);
        tracker.push_back(end);
    }
    sort(tracker.begin(), tracker.end(), sortPointEntry);
    int maxVal = INT_MIN;
    int cur = 0;
    for(auto item: tracker)
    {
        cur += item->val;
        maxVal = max(maxVal, cur);
    }
    return maxVal;
}

/*
 Check if the given unsigned char * is a valid utf-8 sequence.
 
 Return value :
 If the string is valid utf-8, 0 is returned.
 Else the position, starting from 0, is returned.
 
 Valid utf-8 sequences look like this :
 0xxxxxxx
 110xxxxx 10xxxxxx
 1110xxxx 10xxxxxx 10xxxxxx
 11110xxx 10xxxxxx 10xxxxxx 10xxxxxx
 111110xx 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx
 1111110x 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx
 */
int is_utf8(unsigned char *str, size_t len, size_t* first_invalid_pos)
{
    size_t i = 0;
    size_t j = 0;
    size_t continuation_bytes = 0;
    
    while (i < len)
    {
        j = i;
        if (str[i] <= 0x7F)
            continuation_bytes = 0;
        else if (str[i] >= 0xC0 /*11000000*/ && str[i] <= 0xDF /*11011111*/)
            continuation_bytes = 1;
        else if (str[i] >= 0xE0 /*11100000*/ && str[i] <= 0xEF /*11101111*/)
            continuation_bytes = 2;
        else if (str[i] >= 0xF0 /*11110000*/ && str[i] <= 0xF4 /* Cause of RFC 3629 */)
            continuation_bytes = 3;
        else
        {
            if (first_invalid_pos) *first_invalid_pos = j;
            return (int)i+1;
        }
        i += 1;
        while (i < len && continuation_bytes > 0
               && str[i] >= 0x80
               && str[i] <= 0xBF)
        {
            i += 1;
            continuation_bytes -= 1;
        }
        if (continuation_bytes != 0)
        {
            if (first_invalid_pos) *first_invalid_pos = j;
            return (int)i+1;
        }
    }
    return 0;
}

/*
1. Given a sorted array of floats, find the index of the number closest to x:
Example: {1.2, 2.5, 9.3}  x = 5,    return 1
*/
int findClosest(vector<float> input, int target)
{
    if(input.size() == 0) return 0;
    int left = 0;
    int right = (int)input.size()-1;
    while(left<=right)
    {
        int mid = left + (right-left)/2;
        if(input[mid] < target)
        {
            if(mid == input.size()-1)
            {
                return mid;
            }
            else if(input[mid+1]> target)
            {
                if(abs(target-input[mid]) < abs(target-input[mid+1]))
                {
                    return mid;
                }
                else
                {
                    return mid+1;
                }
            }
            else
            {
                left = mid+1;
            }
        }
        else
        {
            if(mid == 0)
            {
                return 0;
            }
            else if(input[mid-1] < target)
            {
                if(abs(target-input[mid-1]) < abs(target-input[mid]))
                {
                    return mid-1;
                }
                else
                {
                    return mid;
                }
            }
            else
            {
                right = mid -1;
            }
        }
    }
    return -1;
}

struct NestedInteger
{
    bool isInteger;
    int getInteger;
    vector<NestedInteger> getList;
};

int calculateSum(NestedInteger s, int level)
{
    if(s.isInteger)
    {
        return s.getInteger * level;
    }
    else
    {
        int sum = 0;
        for(auto item: s.getList)
        {
            sum += calculateSum(item, level+1);
        }
        return sum;
    }
}

int depSum(vector<NestedInteger> input)
{
    int sum = 0;
    for(auto item: input)
    {
        sum += calculateSum(item, 1);
    }
    return sum;
}

vector<vector<int>> threeSumIII(vector<int>& num) {
    vector<vector<int>> result;
    if (num.size() < 3) return result;
    sort(num.begin(), num.end());
    const int target = 0;
    auto last = num.end();
    for (auto a = num.begin(); a < prev(last, 2); ++a) {
        auto b = next(a);
        auto c = prev(last);
        while (b < c) {
            if (*a + *b + *c < target) {
                ++b;
            } else if (*a + *b + *c > target) {
                --c;
            } else {
                result.push_back({ *a, *b, *c });
                ++b;
                --c;
            }
        }
    }
    sort(result.begin(), result.end());
    result.erase(unique(result.begin(), result.end()), result.end());
    return result;
}


/*
 Heap operations
 */
/*
int getParentIndex(int nodeIndex) {
    return (nodeIndex - 1) / 2;
}

int getLeftChildIndex(int nodeIndex) {
    return 2 * nodeIndex + 1;
}

int getRightChildIndex(int nodeIndex) {
    return 2 * nodeIndex + 2;
}

int heapSize;
int data[100];
int arraySize;

void siftUp(int nodeIndex) {
    int parentIndex, tmp;
    if (nodeIndex != 0) {
        parentIndex = getParentIndex(nodeIndex);
        if (data[parentIndex] > data[nodeIndex]) {
            tmp = data[parentIndex];
            data[parentIndex] = data[nodeIndex];
            data[nodeIndex] = tmp;
            siftUp(parentIndex);
        }
    }
}

void insert(int value) {
    if (heapSize == arraySize)
        throw string("Heap's underlying storage is overflow");
    else {
        heapSize++;
        data[heapSize - 1] = value;
        siftUp(heapSize - 1);
    }
}

void siftDown(int nodeIndex) {
    int leftChildIndex, rightChildIndex, minIndex, tmp;
    leftChildIndex = getLeftChildIndex(nodeIndex);
    rightChildIndex = getRightChildIndex(nodeIndex);
    if (rightChildIndex >= heapSize)
    {
        if (leftChildIndex >= heapSize)
            return;
        else
            minIndex = leftChildIndex;
    } else {
        if (data[leftChildIndex] <= data[rightChildIndex])
            minIndex = leftChildIndex;
        else
            minIndex = rightChildIndex;
    }
    
    if (data[nodeIndex] > data[minIndex])
    {
        tmp = data[minIndex];
        data[minIndex] = data[nodeIndex];
        data[nodeIndex] = tmp;
        siftDown(minIndex);
    }
}

void removeMin() {
    if (heapSize == 0)
        throw string("Heap is empty");
    else {
        data[0] = data[heapSize - 1];
        heapSize--;
        if (heapSize > 0)
            siftDown(0);
    }
}
*/

vector<int> x;
int heapSize = 0;

void siftDown(vector<int>& d, int k)
{
    int parentIndex = (k-1)/2;
    if(parentIndex == k) return;
    if(d[parentIndex] > d[k])
    {
        swap(d, parentIndex, k);
        siftDown(d, parentIndex);
    }
}

void insert(int input)
{
    x.push_back(input);
    heapSize++;
    siftDown(x, (int)x.size()-1);
}

void siftUp(vector<int>& d, int k)
{
    int leftChildIndex = 2*k +1;
    int rightChildIndex = 2*k + 2;
    if(leftChildIndex >= d.size()) return;
    if(rightChildIndex >= d.size())
    {
        if(d[k] > d[leftChildIndex])
        {
            swap(d, k, leftChildIndex);
        }
    }
    else
    {
        if(d[leftChildIndex] > d[rightChildIndex])
        {
            swap(d, leftChildIndex, rightChildIndex);
        }
        if(d[k] > d[leftChildIndex])
        {
            swap(d, k, leftChildIndex);
            siftUp(d, leftChildIndex);
        }
    }
}


int get()
{
    if(x.size() >0)
    {
        int res = x[0];
        x[0] = x[x.size()-1];
        x.pop_back();
        siftUp(x, 0);
        return res;
    }
    return INT_MIN;
}

/*
 1. find longest substring which contains n distinct characters.
 */
string findSub(string s, int n)
{
    int count = 0;
    if(n ==0) return "";
    int longestLength = 0;
    int startIndex = -1;
    int endIndex = -1;
    unordered_map<char, vector<int>> tracker;
    int left = 0;
    int right = 0;
    while(right < s.length())
    {
        if(tracker.find(s[right]) == tracker.end())
        {
            count++;
        }
        tracker[s[right]].push_back(right);
        
        if(count==n)
        {
            if(right - left + 1>longestLength)
            {
                startIndex = left;
                endIndex = right;
                longestLength = right - left;
            }
            
        }
        while(count >n)
        {
            tracker[s[left]].erase(tracker[s[left]].begin());
            if(tracker[s[left]].size() == 0)
            {
                count--;
            }
            left++;
        }
        right ++;
    }
    return s.substr(startIndex, endIndex - startIndex + 1);
}

/*
 2. a) 给你棵二叉树，节点上有权值，问从一个叶子走到另外一个叶子的路里面权值最大的那条是什么。
 */
// this method assme the tree have both n1 and n2,
TreeNode* findLCA(TreeNode* root, TreeNode* n1, TreeNode* n2, int& maxWeight)
{
    if(!root) return NULL;
    if(root == n1)
    {
        maxWeight = max(maxWeight, root->val);
        return root;
    }
    if(root == n2) return root;
    TreeNode* leftNode = findLCA(root->left, n1, n2, maxWeight);
    TreeNode* rightNode = findLCA(root->right, n1, n2, maxWeight);
    if(leftNode && rightNode)
    {
        maxWeight = max( maxWeight, root->val);
        return root;
    }
    else
    {
        if(leftNode)
        {
            maxWeight = max(maxWeight, leftNode->val);
        }
        else if(rightNode)
        {
            maxWeight = max(maxWeight, rightNode->val);
        }
        return leftNode? leftNode: rightNode;
    }
}

int findMax(TreeNode* root, TreeNode* n1, TreeNode* n2)
{
    int res = 0;
    TreeNode* lca = findLCA(root, n1, n2, res);
    cout<<lca->val<<endl;
    return res;
    
}

/*
 b) 给你数组a1,a2,...,an。输出数组a2*a3*...*an, a1*a3*a4*...*an, ..., a1*
 a2*...*an-1.
 */
vector<int> multiresult(vector<int> a)
{
    int n = (int)a.size();
    vector<int> left(n,0);
    vector<int> right(n,0);
    vector<int> res(n, 0);
    left[0] = 1;
    for(int i = 1;i< n;i++)
    {
        left[i] = left[i-1] * a[i-1];
    }
    right[n-1] = 1;
    for(int i = n-2;i>=0;i--)
    {
        right[i] = right[i+1] * a[i+1];
    }
    for(int i =0;i<n;i++)
    {
        res[i] = left[i] * right[i];
    }
    return res;
}


/*
 isOneEditDistance 判断两个string是不是只差一个编辑距离
 */
bool isOneEditDistance(string s1, string s2)
{
    int l1 = (int)s1.length();
    int l2 = (int)s2.length();
    if(abs(l1 - l2) > 1) return false;
    int diff = 0;
    int first = 0;
    int second = 0;

    
    while(first < l1 && second < l2)
    {
        if(s1[first] != s2[second])
        {
            if(diff >0) return false;
            else if (l1 > l2)
            {
                diff = 1;
                first ++;
            }
            else if(l2 > l1)
            {
                diff = 1;
                second ++;
            }
            else
            {
                diff = 1;
                first ++;
                second ++;
            }
            
        }
        else
        {
            first++;
            second++;
        }
    }
    
    if( diff == 1)
    {
        if( first == l1 && second == l2) return true;
        else return false;
    }
    else
    {
        if(l1>l2)
        {
            if(first == l1-1 && second == l2) return true;
            else return false;
        }
        else if(l2 > l1)
        {
            if(first== l1 && second == l2-1) return true;
            else return false;
        }
        else
        {
            return diff == 1;
        }
    }
}

bool isOneEditDistanceII(string s1, string s2)
{
    int l1 = (int)s1.length();
    int l2 = (int)s2.length();
    int diff = 0;
    int first = 0;
    int second = 0;
    if(l1 == l2)
    {
        while(first <l1 && second<l2)
        {
            if(s1[first] != s2[second])
            {
                if( diff == 0) diff = 1;
                else return false;
            }
            first ++;
            second ++;
        }
        return diff == 1;
    }
    else
    {
        if(abs(l1-l2)>1) return false;
        string longstr;
        string shortstr;
        if(l1>l2)
        {
            longstr = s1;
            shortstr = s2;
        }
        else
        {
            longstr = s2;
            shortstr = s1;
            int tmp = l1;
            l1 = l2;
            l2 = tmp;
        }
        
        while(first < l1 && second< l2)
        {
            if(longstr[first] != shortstr[second])
            {
                diff = 1;
                first ++;
            }
            else
            {
                first++;
                second++;
            }
        }
        
        if(diff == 1)
        {
            return first == l1 && second == l2;
        }
        else
        {
            return first == l1-1;
        }
    }
}


/*
 http://www.mitbbs.com/article_t/JobHunting/32748027.html
 有个followup的问题就是：因为我不想过多移动这些变量，所以怎么才能设计一个算法
 所需要移动的object最少。
 比如如果变量的size一次是4, 4, 1, 1, 8, 8, 1, 1最好的排法是4, 4, 8, 8, 1, 1,
 1, 1.而不是8 8 4 4 1 1 1 1因为前一种所需要移动的cost最小。这个没想出来了。。
 应该用divide and conquer？
 
 [Quick sort to bring 4, 8 byte element to the front]
 */

void sort(vector<int>& input)
{
    if(input.size() == 0) return;
    int left = 0;
    int right = (int)input.size()-1;
    while(left<=right)
    {
        if(input[left] == 4 || input[left] == 8)
        {
            left ++;
        }
        else
        {
            swap(input, left, right);
            right--;
        }
    }
}

/*
 1. 设计算法找出平面上点的convex hull 不用写code
 */
vector<Point> convexHull(vector<Point>& input)
{
    // 1. sort input base on the x of each point increasing order.
    // 2. for left most point find the SMALLEST angle can formed by
    //    left most point and any other points which are not
    //    selected.
    // 3. pick the point which form the SMALLEST angle in selected.
    //    this point is one point on the covex hull.
    // 4. continue until find the first point.
    // 5. then the convex hull is found.
    return input;
}

/*
 2. code 插入元素到max heap。
 */
void insert_intoMaxHeap(vector<int>& heap, int n, int value)
{
    heap.push_back(value);
    int i = (int)heap.size()-1;
    while(i>0)
    {
        if(heap[i/2] > heap[i])
        {
            break;
        }
        else
        {
            swap(heap, i, i/2);
            i = i/2;
        }
    }
}
/*
 1. 一个bit的stream， 每次读取6个bit。转化成char。
 */

void readerStream(istream& stream)
{
    char* buf;
    char res[8];
    int charIndex = 0;
    int buffIndex = 0;
    int readedlength = 0;
    while(true)
    {
        while(charIndex<8 && buffIndex<readedlength)
        {
            res[charIndex++] = buf[buffIndex++];
        }
        if(charIndex == 8)
        {
            cout<<res;
            charIndex = 0;
        }
        if(buffIndex == readedlength)
        {
            stream.read(buf, 6);
            //readedlength = stream.getchar();// pretend this line will return the size of the
            buffIndex = 0;
        }
    }
}

/*
 
 while(stream.read(buf,6))
 {
 j = 0;
 for( int i = index; i< 8;i++)
 {
 if(j<6)
 {
 res[i] = buf[j++];
 }
 else
 {
 break;
 }
 index = i;
 
 }
 
 if(index == 7)
 {
 // read this res is done.
 cout<<res[0];
 index = 0;
 }
 
 if(j <6 )
 {
 while(j<6)
 {
 res[index++] = buf[j++];
 }
 }
 }
 }
 */

int reader4096(char* buf)
{
    return 4096;
}

char buf[4096];
int startIndex = 0;
int totalread = 0;

void readRandomStream(int n)
{
    char res[n];
    int i = 0;
    while(i<n && startIndex < totalread)
    {
        res[i++] = buf[startIndex++];
        if(startIndex == totalread)
        {
            totalread = reader4096(buf);
            startIndex = 0;
        }
    }
    
}




/*
 写出长度小于N的所有旋转对称数. 例子 689 顺时针旋转180度还是689递归。也可以dp。
 */
string creatRotate(string s, bool odd)
{
    int length = (int)s.length();
    int i;
    if(odd)
    {
        i = length - 2;
    }
    else
    {
        i = length - 1;
    }
    for(int j = i; j>=0;j-- )
    {
        if(s[j] == '9') s.push_back('6');
        else if(s[j] == '8') s.push_back('8');
        else if(s[j] == '6') s.push_back('9');
    }
    return s;
}

void worker(int n, string& s, vector<string>& res, bool odd)
{
    if(s.length() == n)
    {
        if(odd)
        {
            s.push_back('8');
        }
        res.push_back(creatRotate(s, odd));
        if(odd)
        {
            s.pop_back();
        }
        return;
    }
    else
    {
        s.push_back('6');
        worker(n, s, res, odd);
        s.pop_back();
        s.push_back('8');
        worker(n, s, res, odd);
        s.pop_back();
        s.push_back('9');
        worker(n, s, res, odd);
        s.pop_back();
    }
}

vector<string> findAllRotate(int n)
{
    vector<string> res;
    if(n ==0) return res;
    res.push_back("8");
    if(n ==1) return res;
    string s = "";
    for(int i = 2;i<=n;i++)
    {
        s = "";
        worker(i/2, s, res, i%2 == 0? false:true);
    }
    return res;
}



/*
 
 Write a function which, given two integers (a numerator and a denominator), prints the decimal representation of the rational number "numerator/denominator".
 Since all rational numbers end with a repeating section, print the repeating section of digits inside parentheses; the decimal printout will be/must be
 
 Example:
 1 , 3 = 0.(3)
 2 , 4 = 0.5(0)
 22, 7 = 3.(142857)
 */
string divideII(int a, int b)
{
    vector<int> tracker;
    
    int wholenumber = a/b;
    int remainder = a%b;
    while(true)
    {
        if(tracker.size() == 0 || find(tracker.begin(), tracker.end(), remainder) == tracker.end())
        {
            tracker.push_back(remainder);
            remainder = (remainder*10%b);
        }
        else
        {
            break;
        }
    }
    
    string res = "";
    res += to_string(wholenumber);
    res += ".";
    for(int i = 0;i< tracker.size();i++)
    {
        if(tracker[i] == remainder)
        {
            res+="(";
        }
        res += to_string(tracker[i] * 10 / b);
    }
    //if(remainder != 0)
    //{
    res += ")";
    //}
    
    
    return res;
}

/*
 
 You have a binary tree where each node knows the number of nodes in its sub-tree (including itself).
 
 Given a node n and an int k,
 write a function to return the kth
 node in an in order traversal.
 Can you do this non recursively
 */
// in the method assume the val of node is the number of node in the tree.
TreeNode* findkth(TreeNode* root, int k)
{
    int leftPatch = 0;
    if(!root) return root;
    if(root->val < k) return NULL;
    while(leftPatch < k-1 && root)
    {
        if(root->left)
        {
            if(leftPatch + root->left->val == k - 1)
            {
                return root;
            }
            else if(leftPatch + root->left->val < k-1)
            {
                leftPatch += (root->left->val + 1);
                root = root->right;
            }
            else
            {
                root = root->left;
            }
        }
        else
        {
            leftPatch ++;
            root = root->right;
        }
    }
    return NULL;
}

/*
 Given an array of integer, find the number of un-ordered pairs in that array, say given {1, 3, 2},
 the answer is 1 because {3, 2} is un-ordered, and for array {3, 2, 1}, the answer is 3 because {3, 2}, {3, 1}, {2, 1}.
 
 Obviously, this can be solved by brute force with O(n^2) running time, or permute all possible
 pairs then eliminate those invalid pairs.
 
 My question is does any body have any better solution and how would you do it because it seems
 like a dynamic programming problem. A snippet of code would be helpful
 */
int totalcount;
int* mergeII(int a[], int n, int b[], int m)
{
    int* res = (int*)malloc((n+m)* sizeof(int));
    int i = 0;
    int j = 0;
    while(i<n && j<m)
    {
        if(a[i] <= b[j])
        {
            res[i+j] = a[i];
            i++;
        }
        else
        {
            totalcount += (n - i);
            res[i+j] = b[j];
            j++;
        }
    }
    if(i ==n)
    {
        for(int k = n-1+j;k< n+m;k++)
        {
            res[k] = b[j++];
        }
    }
    if(j ==m)
    {
        for(int k = m-1+j;k<n+m;k++)
        {
            res[k] = a[i++];
        }
    }
    return res;
}

int* mergeSort(int a[] , int n)
{
    if(n <= 1) return a;
    int k = n/2;
    int* left = mergeSort(a, k);
    int* right = mergeSort(a+k, n-k);
    int* res = mergeII(left, k, right, n-k);
    return res;
}



int countInversion(int a[], int n)
{
    totalcount = 0;
    mergeSort(a,n);
    return totalcount;
}

/*
 bolts and nuts
 */

int compare (int bolt , int nut)
{
    return bolt-nut;
}

void match(int bolt[], int nut[] , int n)
{
    if(n<=1) return;
    int left = 0;
    int right = n -1;
    while(left<=right)
    {
        if(compare(bolt[0], nut[left]) <0 )
        {
            swap(nut, left, right--);
        }
        else if(compare(bolt[0], nut[left]) ==0)
        {
            //put the matching nut to 0 for now.
            swap(nut, left++, 0);
        }
        else // nut is smaller than bolt.
        {
            left++;
        }
    }
    swap(nut, left-1, 0);
    int index = left-1;
    left =1;
    right = n -1;
    while(left<=right)
    {
        if(compare(bolt[left], nut[index]) <0)
        {
            left++;
        }
        else if(compare(bolt[left], nut[index]) >0)
        {
            swap(bolt, left, right--);
        }
    }
    swap(bolt, 0, index);
    match(bolt, nut, index);
    match(bolt+index+1, nut+index+1, n -index -1);
    
}

/*
 1 + b + 2 = b + 3
 
 或者 （x ＋ 1）＊ 3 ＋ 2 *（2x + 5） 化简成7x + 13
 */
vector<int> mult(vector<int> v1, vector<int> v2)
{
    int l1 = (int)v1.size();
    int l2 = (int)v2.size();
    vector<int> res(l1+l2,0);
    for(int i =0;i< l1 && i < l2;i++)
    {
        res[i] = v1[i]-v2[i];
    }
    if(l1>l2)
    {
        for(int i = l2;i< l1;i++)
        {
            res[i] = v1[i];
        }
    }
    else if(l2>l1)
    {
        for(int i = l1;i<l2;i++)
        {
            res[i] = v2[i];
        }
    }
    return res;
}

vector<int> sub(vector<int> v1, vector<int> v2)
{
    int l1 = (int)v1.size();
    int l2 = (int)v2.size();
    vector<int> res(l1+l2,0);
    for(int i =0;i< l1 && i < l2;i++)
    {
        res[i] = v1[i]-v2[i];
    }
    if(l1>l2)
    {
        for(int i = l2;i< l1;i++)
        {
            res[i] = v1[i];
        }
    }
    else if(l2>l1)
    {
        for(int i = l1;i<l2;i++)
        {
            res[i] = -v2[i];
        }
    }
    return res;
}

vector<int> multiply(vector<int> v1, vector<int> v2)
{
    int l1 = (int)v1.size()-1;
    int l2 = (int)v2.size()-1;
    vector<int> res(l1+l2+1,0);
    for(int i = 0;i<=l1;i++)
    {
        for(int j = 0;j<=l2;j++)
        {
            res[i+j] += v1[i] * v2[j];
        }
    }
    return res;
}

string simplilfy(string s)
{
    return "";
}

/*
 Linkedin
 写一个Stack的API，包括push, pop和findMiddle功能
 */
DoubleLinkedListNode* start = NULL;
int ListCount = 0;
DoubleLinkedListNode* mid = NULL;

void push(int i)
{
    DoubleLinkedListNode* tmp = new DoubleLinkedListNode(i,0);
    tmp->next = start;
    if(start)
    {
        start->prev = tmp;
    }
    start = tmp;
    ListCount ++;
    if(ListCount == 1) mid = start;
    else if(ListCount%2 ==0)
    {
        mid = mid->prev;
    }
    
}

int pop()
{
    int res = -1;
    if(start)
    {
        
        res = start->key;
        DoubleLinkedListNode* tmp = start;
        start = start->next;
        start->prev = NULL;
        delete tmp;
        ListCount --;
        if(ListCount%2 ==1) mid = mid->next;
    }
    return res;
}

int findMid()
{
    if(mid) return mid->key;
    return -1;
}

/*
 Linkedin
 Given a list of child->parent relationships, build a binary tree out of it.
 All the element Ids inside the tree are unique.
 
 Example:
 
 Given the following relationships:
 
 Child   Parent  IsLeft
 15      20      true
 19      80      true
 17      20      false
 16      80      false
 80      50      false
 50      null    false
 20      50      true
 
 
 You should return the following tree:
 50
 /  \
 20   80
 / \   /\
 15 17 19 16
 */

struct inputItem
{
    int child;
    int parent;
    bool isLeft;
    inputItem(int c, int p, bool b) : child(c), parent(p), isLeft(b) {}
};

TreeNode* constructTree(vector<inputItem> input)
{
    unordered_map<int, vector<pair<int,bool>>> tracker;
    TreeNode* root = NULL;
    for(auto item:input)
    {
        if(item.parent <0)
        {
            root = new TreeNode(item.child);
            continue;
        }
        tracker[item.parent].push_back(make_pair(item.child, item.isLeft));
    }
    
    queue<TreeNode*> q;
    q.push(root);
    while(!q.empty())
    {
        TreeNode*tmp = q.front();
        q.pop();
        if(tracker.find(root->val) != tracker.end())
        {
            for(auto item : tracker[root->val])
            {
                TreeNode* c = new TreeNode(item.first);
                q.push(c);
                if(item.second)
                {
                    tmp->left = c;
                }
                else
                {
                    tmp->right = c;
                }
            }
        }
    }
    
    return root;
}


/*
 Linkedin
 打印一个数组的所有乘数组合，从大到小，不要有重复
 */
void workerabc(vector<int> input, int level, vector<int> tmp, vector<vector<int>>& res )
{
    if(level == input.size())
    {
        res.push_back(tmp);
        return;
    }
    else
    {
        workerabc(input, level+1, tmp, res);
        tmp.push_back(input[level]);
        workerabc(input, level+1, tmp, res);
        tmp.pop_back();
    }
}

vector<vector<int>> allMulti(vector<int> input)
{
    vector<vector<int>>res;
    if(input.size() ==0) return res;
    sort(input.begin(), input.end());
    vector<int> tmp;
    workerabc(input, 0, tmp, res);
    return res;
}

/*
 Linkedin
 打印一个数的所有乘数组合，从大到小，不要有重复
 */
vector<int> allFactors(int n)
{
    stack<int> secondHalf;
    vector<int> res;
    for(int i = 1;i<=sqrt(n);i++)
    {
        if(n%i ==0)
        {
            res.push_back(i);
            secondHalf.push(n/i);
        }
    }
    
    while(!secondHalf.empty())
    {
        if(res[res.size()-1] == secondHalf.top())
        {
            secondHalf.pop();
        }
        if(!secondHalf.empty())
        {
            res.push_back(secondHalf.top());
            secondHalf.pop();
        }
    }
    return res;
}

TreeNode* findLCANoparent(TreeNode* root, TreeNode* n1, TreeNode* n2)
{
    if(!root) return NULL;
    if(root == n1)
    {
        return root;
    }
    if(root == n2) return root;
    TreeNode* leftNode = findLCANoparent(root->left, n1, n2);
    TreeNode* rightNode = findLCANoparent(root->right, n1, n2);
    if(leftNode && rightNode)
    {
        return root;
    }
    else
    {
        return leftNode? leftNode: rightNode;
    }
}

/*
 Given n intervals [si, fi], find the maximum number of overlapping intervals
 */
struct EndPoint{
    int val;
    int count;
    EndPoint(int v, int c) : val(v), count(c) {}
};

bool sortEndPoint(EndPoint e1, EndPoint e2)
{
    if(e1.val < e2.val)
    {
        return true;
    }
    else if(e1.val == e2.val)
    {
        return e1.count < e2.count;
    }
    else
    {
        return false;
    }
}

int maxOverlap(vector<Interval> input)
{
    vector<EndPoint> tracker;
    for(auto item : input)
    {
        EndPoint s(item.start, 1);
        EndPoint e(item.end, -1);
        tracker.push_back(s);
        tracker.push_back(e);
    }
    
    sort(tracker.begin(), tracker.end(), sortEndPoint);
    int res = 0;
    int tmp = 0;
    for(auto item: tracker)
    {
        tmp += item.count;
        res = max(res, tmp);
    }
    return res;
}

/*
 电面1：
 expr ::= int | ‘(‘ op expr… ‘)’;
 op ::= ‘+’ | ‘*’;
 
 “( * 1 ( + 1 2 3 ) )” => 6
 “( * ( + 1 1 ) 17 )” => 34
 “7” => 7
 ( * ( + 1 1 ) 17 )
 ( * 17 ( + 1 1 ) )
 operator: *+
 oprands: (1 (1 2 3)
 
 这题特别要求一个运算符可以对应任意个数。
 */
int calculat(vector<string> input)
{
    stack<string> operators;
    stack<string> operands;
    
    for(auto item : input)
    {
        if(item == "*" || item == "+")
        {
            operators.push(item);
        }else
        {
            if(item != ")")
            {
                operands.push(item);
            }
            else
            {
                vector<int> tmp;
                while(operands.top() != "(")
                {
                    tmp.push_back(atoi(operands.top().c_str()));
                    operands.pop();
                }
                operands.pop();
                int res;
                if(operators.top() == "*")
                {
                    res = 1;
                    for(auto i: tmp)
                    {
                        res *= i;
                    }
                }
                else
                {
                    res = 0;
                    for(auto i: tmp)
                    {
                        res += i;
                    }
                }
                operators.pop();
                operands.push(to_string(res));
            }
        }
        
    }
    if(operands.size() == 1)
    {
        return atoi(operands.top().c_str());
    }
    return INT_MIN;
}

int longestCommonSubString(string s1, string s2)
{
    int n = (int)s1.length();
    int m = (int)s2.length();
    vector<vector<int>> tracker ( n+1, vector<int>(m+1, 0));
    
    for(int i = 0;i<=m;i++)
    {
        tracker[n][i] = 0;
    }
    
    for(int i = 0;i<=n;i++)
    {
        tracker[i][m] = 0;
    }
    
    int res = 0;
    for(int i = n-1;i>=0;i--)
    {
        for(int j = m-1;j>=0;j--)
        {
            if(s1[i] == s2[j])
            {
                tracker[i][j] = tracker[i+1][j+1] + 1;
                res = max(res, tracker[i][j]);
            }
            else
            {
                tracker[i][j] = 0;
            }
            
        }
    }
    return res;
}

/*
 把一个iterator的iterator 转换成 iterator
 */

vector<vector<int>> backend;
vector<vector<int>>::iterator currentParentIterator;
vector<int> ::iterator currentChildIterator;

bool hasNext()
{
    while(true)
    {
        if(currentChildIterator == (*currentParentIterator).end())
        {
            if(currentParentIterator == backend.end())
            {
                return false;
            }
            else
            {
                currentParentIterator++;
                if(currentParentIterator == backend.end())
                {
                    currentChildIterator = (*currentParentIterator).end();
                    return false;
                }
                currentChildIterator = (*currentParentIterator).begin();
            }
        }
        else
        {
            return true;
        }
    }
}

int getNext()
{
    if(hasNext())
    {
        int res;
        res = *currentChildIterator;
        currentChildIterator++;
        return res;
    }
    else
    {
        return INT_MIN;
    }
}



/*
 2.  Given a binary tree where all the right nodes are leaf nodes, flip it upside down and turn it into a tree
 with left leaf nodes.
 */
TreeNode* InverseATree(TreeNode* root)
{
    TreeNode* tmp = NULL;
    TreeNode* newRoot = NULL;
    TreeNode* currRight = NULL;
    TreeNode* prevright = NULL;
    
    while(root)
    {
        currRight = root->right;
        tmp = root;
        root = root->left;
        tmp->left = newRoot;
        tmp->right = prevright;
        

        newRoot = tmp;
        prevright = currRight;
    }
    return newRoot;
}

/*
 search a number in rotated sorted array (leetcode)
 */
int search(vector<int> input, int target)
{
    int left = 0;
    int right = (int)input.size()-1;
    while(left<=right)
    {
        int mid = left + (right-left) / 2;
        if(target == input[mid])
        {
            return mid;
        }
        
        if(input[left] < input[mid])
        {
            if(target >= input[left] && target <=input[mid])
            {
                right = mid - 1;
            }
            else
            {
                left = mid + 1;
            }
        }
        else if (input[left] > input[mid])
        {
            if(target >= input[mid] && target <= input[right])
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
            left++;
        }
        
    }
    return false;
}

/*
 2. Decide whether a target is covered by a list of intervals (类似merge intervals)
 */
struct IntervalTreeNode{
    Interval val;
    IntervalTreeNode* left;
    IntervalTreeNode* right;
    int max;
    IntervalTreeNode(Interval input, int m): val(input), max(m) {}
};

bool isCovered(IntervalTreeNode* root, int target)
{
    if(!root || target > root->max) return false;
    if(target >= root->val.start && target <= root->val.end)
    {
        return true;
    }
    if(root->left && target <= root->left->max)
    {
        return isCovered(root->left, target);
    }
    else
    {
        return isCovered(root->right, target);
    }
}

/*
 2. 洗牌 要求in-place
 */
void shuffle(int input[], int n)
{
    for(int i = 0;i<n;i++)
    {
        int tmp = rand() % (n-i);
        swap(input, i, i +tmp);
    }
}

/*
 * Return the smallest character that is strictly larger than the search
 character,
 * ['c', 'f', 'j', 'p', 'v'], 'a' => ‘c’
 */
char findStrictly(char input[], int n, char target)
{
    int left = 0;
    int right = n-1;
    while(left <= right)
    {
        int mid = left + (right - left) / 2;
        if(input[mid] > target && (mid == 0 || target > input[mid-1]))
        {
           return input[mid];
        }
        else if(target > input[mid] && (mid <n-1 && target < input[mid + 1]))
        {
            return input[mid + 1];
        }
        else if( target < input[mid])
        {
            right = mid - 1;
        }
        else
        {
            left = mid + 1;
        }
    }
    return '1';
}

/*
 2.binary tree, print all paths from root to leaf.
 */
vector<vector<int>> printAllpath(TreeNode* root)
{
    stack<TreeNode*> tracker;
    vector<vector<int>> res;
    TreeNode* prev = NULL;
    tracker.push(root);
    while(!tracker.empty())
    {
        TreeNode* tmp = tracker.top();
        if(tmp->left &&
           (prev == NULL || (prev != tmp->left && prev != tmp->right)))
        {
            tracker.push(tmp->left);
        }
        else if(prev == tmp->right || !(tmp->right))
        {
            tracker.pop();
            prev = tmp;
        }
        else if(tmp->right)
        {
            tracker.push(tmp->right);
        }
        
        tmp = tracker.top();
        if(!(tmp->left) && !(tmp->right))
        {
            //insert the stack to res vector.
        }
    }
    return res;
}

/*
 Clone graph
 */
/*
vector<GraphNode*> cloneGraph(vector<GraphNode*> input)
{
    unordered_map<GraphNode, GraphNode> tracker;
    for(auto node : input)
    {
        tracker[node] = new GraphNode(node->val);
    }
    
    for(auto node: input)
    {
        for(auto neighbor : node->neighbors)
        {
            tracker[node]->neighbors.push_back(tracker[neighbor]);
        }
    }
    
    vector<GraphNode*> res;
    for(auto node : input)
    {
        res.push_back(tracker[node]);
    }
    
    return res;
}
 */

/*
 2. find the maximum number in an integer array. The numbers in the array increase first, then decreases. Maybe there’s only increase or decrease. 先说了直接扫一遍，是O(n), 然后用binary search 就是O(log n).最后时间不够，没写完，应该是挂了。
 */
int findMax(int input[], int n)
{
    if(n <2) return 0;
    int left = 0;
    int right = n-1;
    while(left<=right)
    {
        int mid = left + (right-left)/2;
        if((mid == 0 && input[mid] > input[mid+1]) ||
           (mid == n-1 && input[mid] > input[mid-1]) ||
           (input[mid] > input[mid-1] && input[mid] > input[mid+1]))
        {
           return mid;
        }
        else if(input[mid+1]> input[mid] )
        {
            left = mid + 1;
        }
        else if(input[mid+1] < input[mid])
        {
            right = mid-1;
        }
    }
    return -1;
}

/*
 string 有多少palindrome substring
 */
int countOfPalindrome(string s1)
{
    int res = 0;
    int n = (int)s1.length();
    for(int i = 0;i<n;i++)
    {
        int left = i-1;
        int right = i+1;
        while(left>=0 && right<n && s1[left] == s1[right])
        {
            res++;
            left--;
            right++;
        }
        left = i;
        right = i+1;
        while(left>=0 && right<n && s1[left] == s1[right])
        {
            res++;
            left--;
            right++;
        }
    }
    return res+n;
}

/*
 example tree
 
 a
 / \
 b   c
 / \   \
 d   g   z
 \     /
 e   i
 /
 q
 /
 x
 /
 x1
 /
 x2
 
 
 sample output
 
 x2
 d x1
 b e x
 a g q
 c
 z
*/
struct TreeNodeWithColumn{
    TreeNode* node;
    int column;
    TreeNodeWithColumn(TreeNode* n, int c): node(n), column(c) {}
};

vector<vector<int>> printout(TreeNode* root)
{
    queue<TreeNodeWithColumn> tracker;
    tracker.push(TreeNodeWithColumn(root,0));
    unordered_map<int, vector<int>> container;
    int mincolumn = INT_MAX;
    int maxcolumn = INT_MIN;
    while(!tracker.empty())
    {
        TreeNodeWithColumn tmp = tracker.front();
        tracker.pop();
        container[tmp.column].push_back((tmp.node)->val);
        mincolumn = min (mincolumn, tmp.column);
        maxcolumn = max (maxcolumn, tmp.column);
        if(tmp.node->left)
        {
            tracker.push(TreeNodeWithColumn(tmp.node->left, tmp.column-1));
        }
        if(tmp.node->right)
        {
            tracker.push(TreeNodeWithColumn(tmp.node->right, tmp.column+1));
        }
    }
    vector<vector<int>> res;
    for(int i = mincolumn;i<=maxcolumn;i++)
    {
        res.push_back(container[i]);
    }
    return res;
}
/*
1. is valid palindrome
*/
bool isValidPalindrome(string s)
{
    int left = 0;
    int right = (int)s.length()-1;
    while(left<=right)
    {
        if(s[left] == s[right])
        {
            left++;
            right--;
        }
        else
        {
            return false;
        }
    }
    return true;
}


/*
 先说了说简历上的project，然后做题
 给一个set，里面是一堆pair，每个pair里是两个string，一个first，一个second，假设这堆pair能够构成一个树状结构，按照一定的格式打印这棵树
 first-second关系类似paretnt-child关系
 eg
 set: (a, b) (b, c) (a, d) (d, e) (d, f) (d, g)
 树状结构是root = a, root.left = b, root.right = d blah blah
 打印结果：[space] 就是一个空格
 a
 [space]b
 [space][space]c
 [space]d
 [space][space]e
 [space][space]f-google 1point3acres
 [space][space]g
 */

void printTree(TreeNode* root, int n)
{
    for(int i = 0;i<n;i++)
    {
        cout<<" ";
    }
    cout<<root->val;
    cout<<endl;
    if(root->left) printTree(root->left, n+1);
    if(root->right) printTree(root->right, n+1);
}

void printTree(vector<Interval> input)
{
    unordered_map<int, vector<int>> tracker;
    unordered_set<int> child;
    for(auto item: input)
    {
        tracker[item.start].push_back(item.end);
        child.insert(item.end);
    }
    int rootkey = -1;
    for(auto root: tracker)
    {
        if(child.find(root.first) == child.end())
        {
            rootkey = root.first;
            break;
        }
    }
    
    TreeNode* root = new TreeNode(rootkey);
    queue<TreeNode*> q;
    q.push(root);
    while(!q.empty())
    {
        TreeNode* tmp = q.front();
        q.pop();
        vector<int> children = tracker[tmp->val];
        TreeNode* left = NULL;
        TreeNode* right = NULL;
        if(children.size()>0)
        {
            left = new TreeNode(children[0]);
        }
        if(children.size() > 1)
        {
            right = new TreeNode(children[1]);
        }
        tmp->left = left;
        tmp->right = right;
        if(left) q.push(left);
        if(right) q.push(right);
    }
    
    printTree(root, 0);
}

/*
pow(int x, int y)
1<=x<=9
1<=y<=99
想了半天发现数太大了int或者long都装不下，只能返回string
*/

string multiplyStr(string num1, string num2) {
    if(num1 == "0" || num2 == "0") return "0";
    vector<int> n1;
    vector<int> n2;
    for(int i = (int)num1.size()-1;i>=0;i--)
    {
        n1.push_back(num1[i]-'0');
    }
    for(int i = (int)num2.size()-1;i>=0;i--)
    {
        n2.push_back(num2[i]-'0');
    }
    vector<int> res(n1.size() + n2.size() + 1, 0);
    int carry = 0;
    for(int i = 0;i<n1.size();i++)
    {
        carry = 0;
        for(int j = 0;j<n2.size();j++)
        {
            int tmp = carry + n1[i] * n2[j] + res[i+j];
            carry = tmp /10;
            res[i+j] = tmp %10;
        }
        res[i+n2.size()] += carry;
    }
    while (res[res.size()-1] ==0) {
        res.pop_back();
    }
    reverse(res.begin(), res.end());
    string s;
    for(int i = 0;i<res.size();i++)
    {
        s.append(1, '0' + res[i]);
    }
    return s;
}

string pow(int x,int y)
{
    if(y == 0) return "1";
    if(y == 1) return to_string(x);
    
    string tmp = pow(x,y/2);
    
    if(y&1)
    {
        return multiplyStr(multiplyStr(tmp, tmp), to_string(x));
    }
    else
    {
        return multiplyStr(tmp, tmp);
    }
}

/*
 给string a, string b,判断b里面是否存在子字符串是a的anagram。
 */

bool isEqual(unordered_map<char, vector<int>> map1, unordered_map<char, vector<int>> map2)
{
    if(map1.size() != map2.size()) return false;
    for(auto item: map1)
    {
        vector<int> v1 = item.second;
        vector<int> v2 = map2[item.first];
        if(v1.size() != v2.size()) return false;
        for(int i = 0;i< v1.size();i++)
        {
            if(v1[i] != v2[1]) return false;
        }
    }
    return true;
}

bool isSubAnagram(string a, string b)
{
    if(b.length() < a.length()) return false;
    unordered_map<char, vector<int>> amap;
    unordered_map<char, vector<int>> bmap;
    for(int i = 0;i<(int)a.length();i++)
    {
        amap[a[i]].push_back(i);
        bmap[b[i]].push_back(i);
    }
    
    int k = (int)a.length();
    while(k<(int)b.length())
    {
        if(!isEqual(amap, bmap)) return false;
        bmap[b[k-a.length()]].erase(bmap[b[k-a.length()]].begin());
        bmap[b[k]].push_back(k);
    }
    return false;
}




vector<int> container;
int k;
vector<int> tracker;
void setContainer(vector<int> input, int count)
{
    container = input;
    k = count;

    for(int i = 0;i<k;i++)
    {
        if(i != k-1)
            tracker.push_back(i);
        else
            tracker.push_back(i-1);
    }
}

bool hasNextPermutation()
{
    int n = (int)container.size();
    for(int i = 0;i<k;i++)
    {
        if(tracker[k-1-i] != n-1-i)
        {
            return true;
        }
    }
    return false;
}

void printNext()
{
    int n = (int)container.size();
    if(hasNextPermutation())
    {
        for(int i = 0;i<k;i++)
        {
            int tmp = tracker[k-1-i] + 1;
            if(tmp<(n-i))
            {
                tracker[k-1-i] = tmp;
                for(int j = k-i;j<k;j++)
                {
                    tracker[j] = ++tmp;
                }
                break;
            }
        }
    }
    for(int i = 0;i<tracker.size();i++)
    {
        cout<< tracker[i]<<"  ";
    }
    cout<<endl;
}


/*
 是给一个int[] array, e.g {1,5,0,6}和一个int target，e.g. target = 21;
 问是否存在某种分法把array分成几部分，每部分看成一个int，这几部分加起来等于target。
 e.g. {1,5}{0}{6},三部分加起来是21。{1,5}{0,6}也是21。target=25则false
 should be able to do it in O(n^2).
 tried: can not do O(n^2), but n!
 */
int construct(vector<int> input, int startIndex, int endIndex)
{
    int res = 0;
    for(int i = startIndex;i<= endIndex;i++)
    {
        res = res*10 + input[i];
    }
    return res;
}

bool canSeparate(vector<int> input, int startIndex, int target)
{
    if(startIndex >= input.size()) return false;
    if(construct(input, startIndex, (int)input.size()-1) == target) return true;
    for(int i = startIndex + 1;i<input.size();i++)
    {
        int tmp = construct(input, startIndex, i);
        if(canSeparate(input, i+1, target-tmp))
        {
            return true;
        }
    }
    return false;
}


bool canSeparate(vector<int> input, int target)
{
    return canSeparate(input, 0, target);
}

/*
 写一个shuffle数组的算法
 给一个数组A[]，里面的数是无序的，比如 A[]={1,2,3,7,3,2,5}
 shuffle完之后，使之满足A[0]<=A[1]>=A[2]<=A[3]..
 比如一种可能的结果是 A[]={1,7,3,5,2,2,1}
 
 void shuffle(int A[], int n)
 {
 }
 写出算法，并且证明其正确性。
 */
vector<int> shuffleSort(vector<int> input)
{
    if(input.size() <3) return input;
    sort(input.begin(), input.end());
    for(int i = 2;i< input.size();i = i+2)
    {
        swap(input, i, i-1);
    }
    return input;
}

void shuffleOn(int A[],int n)
{
    int increase=1;
    for(int i=1;i<n;i++)
    {
        int local=A[i]>=A[i-1]?1:0;
        if(increase!=local) swap(A[i],A[i-1]);
        increase^=1;
    }
}

/*
 1.将一个数组right rotate k次。要求O(N),in-place
 */
void mirror(int input[], int start, int end)
{
    while(start <=end)
    {
        swap(input, start, end);
        start++;
        end--;
    }
}

void rightRotate(int input[], int n, int k)
{
    k = k%n;
    mirror(input,0, n-1-k);
    mirror(input,n-k,n-1);
    mirror(input, 0, n-1);
}

/*
 2sum
 数字有重复，比如如果sum是10，{2,2,2,8,8}里面算两个(2,8)pair。求pair总数。
 */
int countPair(vector<int> input, int target)
{
    unordered_map<int, int> tracker;
    for(auto item: input)
    {
        tracker[item]++;
    }
    int res = 0;
    
    for(auto item: input)
    {
        if(tracker[item]>0 &&
           tracker.find(target -item) != tracker.end()
           && tracker[target-item] > 0)
        {
            res ++;
            tracker[item]--;
            tracker[target-item] --;
        }
    }
    return res;
}

/*
 "1 3 4 5 6"
 "3 4 5 6 3"
 "4 5 6 3 3"
 ...
 每行是一个包含数字的string。去除所有数字完全重复的strings.比如这里的第二和第三行数字完全相同，
 可以合并成一个。要求合并所有数字完全重复的strings。最后表示对我的优化结果不满意。
*/
vector<string> Reduce(vector<string> input)
{
    unordered_set<string> tracker;
    vector<string> res;
    for(auto s : input)
    {
        string tmp = s;
        sort(tmp.begin(), tmp.end());
        if(tracker.find(tmp) == tracker.end())
        {
            tracker.insert(tmp);
            res.push_back(s);
        }
    }
    return res;
}

string intToRoman(int num) {
    const int radix[] = {1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1};
    const string symbol[] = {"M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"};
    string roman;
    for (size_t i = 0; num > 0; ++i)
    {
        int count = num / radix[i];
        num %= radix[i];
        for (; count > 0; --count) roman += symbol[i];
    }
    return roman;
}

/*
 2.给你一个很大的字典。一对词如果不share 任何字母，比如dog, cat不share字母，而dog, 
 boy就share一个o，则是interesting pair.找出所以interesting pairs中长度乘积最大的pair.输出这个乘积。
 */
unordered_map<char,int> getChar()
{
    unordered_map<char, int> res;
    for(char c = 'a'; c <= 'z'; c++)
    {
        int tmp = 1<<(c-'a');
        res[c] = tmp;
    }
    return res;
}


int findMax(vector<string> input)
{
    unordered_map<char, int> map = getChar();
    unordered_map<string, int> tracker;
    for(auto item : input)
    {
        int tmp = 0;
        for(auto c: item)
        {
            tmp = tmp | map[c];
        }
        tracker[item] = tmp;
    }
    
    int res = INT_MIN;
    for(int i = 0;i< input.size() -1;i++)
    {
        for(int j = i+1;i<input.size();i++)
        {
            if((tracker[input[i]] & tracker[input[j]]) == 0)
            {
                res = max(res, ((int)input[i].length() * (int)input[j].length()));
            }
        }
    }
    return res;
}

/*
 简单寒暄一下，然后问了一个找硬币的题，应该属于背包问题，首先给了一个贪心算法，
 然后他指出不是总能得到最优解，然后我也发现了，他就让我举例，举了个例子，然后开始想动态规划，
 然后想出来了，在板上写了code，分析了一下复杂度，然后一起讨论一下边边角角的问题，
 然后11:15的时候，第二个面试官到了，然后就冲冲结束了。
 I guess this is the question of give a list of type of coins and a target
 find the least count of coins can combine to target.
 */
int findChange(vector<int> input, int target)
{
    int res;
    unordered_set<int> coins;
    vector<int> tracker(target+1,0);
    for(auto item: input)
    {
        coins.insert(item);
        
    }
    for(int i = 1;i<=target;i++)
    {
        res = INT_MAX;
        for(int j = i-1;j>=0;j--)
        {
            if(coins.find(i-j) != coins.end())
            {
                res = min(res, tracker[j]+1);
            }
        }
        tracker[i] = res;
    }
    return tracker[target];
}

/*
 比如说常见的俄罗斯方块，每一个图案都是由4个block组成，现在给定一个N表示N个block，
 把所有有效的俄罗斯方块组合都输出出来，（有 效的是指block是横着或者竖着连接的，不是直接斜着连接）
 I should start from top left cornor and only do top line and only go down
 or go right, and j has to >= to i. other wise it will duplicate.
 */
void worker(vector<vector<vector<int>>>& res, vector<vector<int>>& tracker,int i, int j, int level, int target)
{
    if(i < j) return;
    if(tracker[i][j]) return;

    tracker[i][j] = 1;
    if(level == target-1)
    {
        res.push_back(tracker);
        tracker[i][j] = 0;
        return;
    }
    if(i+1<target)
    {
        worker(res, tracker, i+1, j, level+1, target);
    }
    if(i-1>=0)
    {
        worker(res, tracker, i-1, j, level+1, target);
    }
    if(j+1<target)
    {
        worker(res, tracker, i, j+1, level+1, target);
    }
    if(j-1>=0)
    {
        worker(res, tracker,i, j-1, level+1, target);
    }
    tracker[i][j] =0;    
}

bool flip(vector<vector<int>> input, vector<vector<int>> output)
{
    vector<vector<int>> res = input;
    // flip the piece and rotate it 180 degress make sure no dip.
    // then return the new piece.
    
    if( res != input)
    {
        output = res;
        return true;
    }
    return false;
}

vector<vector<vector<int>>> FindAllRussiaBlocks(int target)
{
    vector<vector<int>> tracker(target, vector<int>(target,0));
    vector<vector<vector<int>>> res;
    if(target <1) return res;
    if(target == 1)
    {
        vector<vector<int>> tracker(1, vector<int>(1,1));
        res.push_back(tracker);
        return res;
    }
    worker(res, tracker, 0, 0, 0 , target);
    vector<vector<vector<int>>> mirror;
    for(auto item: res)
    {
        
        vector<vector<int>> tmp(target, vector<int>(target, 0));
        if(flip(item, tmp))
        {
            mirror.push_back(tmp);
        }
    }
    for(auto item: mirror)
    {
        res.push_back(item);
    }
    
    return res;
}

int matrixMult(vector<int> p)
{
    int n = (int)p.size();
    vector<vector<int>> tracker(n-1, vector<int>(n-1, 0));
    for(int l = 2;i< n;i++)
    {
        for(int i = 0;i<n-l+1;i++)
        {
            int j = i+l-1;
            int tmp = INT_MAX;
            for(int k = i;k<j;k++)
            {
                tmp = min (tmp, tracker[i][k] + tracker[k+1][j] + p[i]*p[k]*p[j]);
            }
            tracker[i][j] =tmp;
            
        }
    }
    return tracker[0][n-1];
}

/*
 一个array里面找最大的这样的h:有h个数大于等于h
 */
bool func(int i1, int i2)
{
    if(i1<i2)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int findMaxH(vector<int> input)
{
    sort(input.begin(), input.end(), func);
    int res = 0;
    for(int i = 0;i<(int)input.size();i++)
    {
        if(input[i] == i)
        {
            res = i;
        }
        
    }
    return res;
}

/*
 二面：十进制十八进制转换，十八进制加法
 */
char covertto18(int n)
{
    if(n >=18) return ' ';
    if(n<10) return '0' + n;
    if(n>=10 && n< 18) return 'a' + (n-10);
    else return ' ';
}

string convertFrom10To18(int input)
{
    string res ="";
    while(input)
    {
        int tmp = input%18;
        input = input/18;
        res = res+covertto18(tmp);
    }
    reverse(res.begin(), res.end());
    return res;
}

/*
 三面：majority element
 */
int findMajority(vector<int> input)
{
    int res=INT_MAX;
    int count = 0;
    for(int i = 0;i< (int)input.size();i++)
    {
        if(input[i] == res)
        {
            count ++;
        }
        if(input[i] != res)
        {
            
            count--;
            if(count<=0)
            {
                count = 0;
                res = input[i];
            }
        }
    }
    if(count>0)
    {
        return res;
    }
    return INT_MIN;
}

/*
一轮。 给定一个binary search tree， 知道到第k大的数。
*/
TreeNode* findKthTreeNode(TreeNode* root, int k)
{
    if(!root || k > root->val) return NULL;
    
    while(k>0 && root)
    {
        if(!(root->left))
        {
            if(k == 1) return root;
            k--;
            root = root->right;
        }
        else if(root->left->val == k-1)
        {
            return root;
        }
        else if(root->left->val >= k)
        {
            root = root->left;
        }
        else
        {
            k = k - root->left->val - 1;
            root = root->right;
        }
    }
    return NULL;
}

string get_decimal(int num, int den) {
    string ret = to_string(num / den);
    ret.push_back('.');
    num %= den;
    unordered_map<int,int> rems;
    
    while(num != 0 && rems.find(num) == rems.end()) {
        rems[num] = (int)ret.size();
        num *= 10;
        int tmp = num/den;
        ret.push_back(tmp + '0');
        num %= den;
    }
    
    if (num != 0) {
        ret.insert(ret.begin() + rems[num], '(');
        ret += ")";
    } else {
        ret += "(0)";
    }
    return ret;
}

/*
 1. 有这么个游戏，举个例子：给你5个空_ _ _ _ _, 每次猜一个字母，这里出题人想让你猜出来clock，
 假如你猜a，告诉你这里面没有。你又猜c，他把c全写出来，所以你有c _ _ c _。 让你最多猜10次。
 写一个程序去猜。输入是几个空，要考虑每次猜的反馈，尽量把词猜出来。
 */
// guessedWord can be "_c__c_" in this case c has already been guessed.

bool isAmatch(string word, string guessedWord)
{
    if(word.length() != guessedWord.length()) return false;
    for(int i = 0;i<word.length();i++)
    {
        if(guessedWord[i] != '_' &&
           guessedWord[i] != word[i])
        {
            return false;
        }
    }
    return true;
}

char guess(vector<string> dict, string guessedWord, unordered_set<char> guessedChar)
{
    vector<string> res;
    for(auto item :dict)
    {
        if(isAmatch(item, guessedWord))
        {
            res.push_back(item);
        }
    }
    
    vector<vector<int>> count(guessedWord.length(), vector<int>(26, 0));
    for(auto item : res)
    {
        for(int i = 0;i< guessedWord.length();i++)
        {
            if(guessedWord[i] == '_')
            {
                count[i][item[i] - 'a']++;
            }
        }
    }
    
    char c = '0';
    int x = INT_MAX;
    for(auto item: count)
    {
        for(int i = 0;i< 26;i++)
        {
            if(item[i] > x && guessedChar.find('a'+i) == guessedChar.end())
            {
                x = item[i];
                c = 'a' + i;
            }
        }
    }
    
    return c;
}

/*
             1 2 2 3 (5)
             3 2 3 (4) (4)
 pacific     2 4 (5) 3 1              Atlantic
             (6) (7) 1 4 5
             (5) 1 1 2 4
 
 每个数字代表该地区的海拔，然后西边是太平洋，东边是大西洋，让我返回所有path，
 每个path能连通大西洋和太平洋，水只能从高处往低处走。
 */

// O(n^3)?
struct altitude
{
    int val;
    bool p;
    bool a;
    int x;
    int y;
    altitude(int v) : val(v), p(false), a(false) {}
};

/*
void findAllpath(vector<vector<altitude>> input)
{
    queue<altitude> q;
    for(int i = 0;i< input.size();i++)
    {
        queue<altitude>  empty;
        unordered_set<altitude> visited;
        swap(q, empty);
        input[i][0].p = true;
        q.push(input[i][0]);
        while(!q.empty())
        {
            altitude tmp = q.front();
            if(tmp.x -1 >=0 &&
               input[tmp.x-1][tmp.y].p &&
               input[tmp.x][tmp.y].val < input[tmp.x-1][tmp.y].val)
            {
                q.push(input[tmp.x-1][tmp.y]);
                input[tmp.x-1][tmp.y].p = true;
                
            }
            // do the same thing for other 3 directions.
        }
        
    }
    
    // after we done with pacific
    // do the same thing again for atlantic.
    
    //count all the element which has both
    
}
*/

/*
 第二题，给你一串正整数。1，2，3.。。。10， 11，12.。。。 给你一个int n，要你返回哪一位的数。比如 给你10，返回的就是1.给11，返回的就是0
 */
int findNthDigit(int n)
{
    int k = 0;
    while(n> 9*(int)std::pow(10,k))
    {
        n = n-9*(int)std::pow(10,k);
        k ++;
    }
    
    int index = n/(k+1) -1 ;
    int offset = n%(k+1);
    
    if(offset == 0)
    {
        int num = std::pow(10,k) + index;
        return num%10;
    }
    else
    {
        int num = std::pow(10,k) + index + 1;
        int formRight = k+2-offset;
        int res = 0;
        for(int i = 0;i< formRight;i++)
        {
            res = num%10;
            num = num/10;
        }
        return res;
    }
    
}

/*
 拿着我简历进来，说有人跟你谈过你的简历吗，我说没有，他表示万分惊讶，然后在我
 简历上挑了一个research project让我说说，说完后用c++出了一个题，一个cipher类
 ，有一个member function是对输入加密，加密方法为对input的每16个Byte和一个
 increasing counter做xor，这个increasing counter也是有16Byte，从00..01（前
 15Byte都是0，最后1Byte是1）开始，还有一个要求，举例说：
 第一个input 有20个Byte，前16个Byte就和00..01做xor，后4个Byte和00..02的前
 4Byte做xor
 然后之后再对第二个input加密的时候，对这个input的前12Byte用00..02的后12Byte（
 即11个Byte 0，1个Byte 1）
 */

int offset = 0;
char cipher[16];

void initialize()
{
    for(int i = 0;i< 16;i++)
    {
        cipher[i] = 0;
    }
    cipher[15]=1;
}

void getNextCipher()
{
    // add 1 to cipher
}

void encrypt(char input[], int n)
{
    char res[n];
    for(int i = 0;i<n;i++)
    {
        if(offset >= 16)
        {
            getNextCipher();
        }
        res[i] = input[i] ^ cipher[offset++];
    }
}

/*
1：print all path from root to leaf.
*/
void print(stack<TreeNode*> input)
{
    stack<TreeNode*> tmp;
    while(!input.empty())
    {
        tmp.push(input.top());
        input.pop();
    }
    
    while(!tmp.empty())
    {
        cout<<tmp.top()->val<<" ";
        tmp.pop();
    }
    
    cout<<endl;
}

void printAllPath(TreeNode* root)
{
    stack<TreeNode*> s;
    s.push(root);
    TreeNode* prev = NULL;
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        if(tmp->left &&
           prev != tmp->left &&
           prev != tmp->right)
        {
            s.push(tmp->left);
        }
        else if (tmp->right &&
                 tmp->right != prev)
        {
            s.push(tmp->right);
        }
        else
        {
            if(!tmp->left && !tmp->right)
            {
                print(s);
            }
            s.pop();
            prev = tmp;
            
        }
            
    }
}


/*
 2:  power set.
 */

void powerSetWorker(int level, vector<int> input, vector<int>& tracker, vector<vector<int>>& res)
{
    if(level == (int)input.size())
    {
        res.push_back(tracker);
        return;
    }
    powerSetWorker(level+1, input, tracker, res);
    tracker.push_back(input[i]);
    powerSetWorker(level+1, input, tracker, res);
    tracker.pop_back();
}

vector<vector<int>> powerSet(vector<int> input)
{
    vector<vector<int>> res;
    vector<int> tracker;
    powerSetWorker(0, input, tracker, res);
    return res;
}

/*
 3:  given a list of words, find palindrome pairs
 */
bool isPalindromeString(string s)
{
    int n = (int)s.length();
    int left = 0;
    int right = n-1;
    while(left < right)
    {
        if(s[left++] != s[right--]) return false;
    }
    return false;
}

vector<pair<string, string>> findPalindrome(vector<string> input)
{
    vector<pair<string, string>> res;
    for(int i = 0;i< input.size()-1;i++)
    {
        for(int j = i+1;j<input.size();j++)
        {
            string tmp = input[i] + input[j];
            if(isPalindrome(tmp))
            {
                res.push_back(make_pair(input[i], input[j]));
                continue;
            }
            else
            {
                tmp = input[j] + input[i];
                if(isPalindrome(tmp))
                {
                    res.push_back(make_pair(input[i], input[j]));
                    continue;
                }
            }
        }
    }
    return res;
}

/*
 4:  implement BufferedReader
 */
int read4K(char** s)
{
    char tmp[10];
    *s = tmp;
    for(int i = 0;i< 9;i++)
    {
        tmp[i] = i + '0';
        
    }
    tmp[9] = '\0';
    return 10;
}

int totalLength = 0;
int currentIndex = 0;
int readedLength = 0;
char* str;

char* read(int n)
{

    //read4K(&x);
    char* res = (char*)malloc(n);
    readedLength = 0;
    while(readedLength<n)
    {
        if(currentIndex < totalLength)
        {
            while(readedLength < n && currentIndex < totalLength)
            {
                res[readedLength++] = *(str+currentIndex++);
            }
        }
        else
        {
            totalLength = read4K(&str);
            currentIndex = 0;
        }
    }
    return res;
}

string readLine()
{
    string res = ""; // talking in java this should stringbuilder....
    readedLength = 0;
    while(str[currentIndex] != '\n') // this can be changed to str
    {
        if(currentIndex < totalLength)
        {
            while(currentIndex < totalLength)
            {
                res += *(str+currentIndex++);
            }
        }
        else
        {
            if(totalLength < 4096)
            {
                res += "EOF";
                break;
            }
            totalLength = read4K(&str);
            currentIndex = 0;
        }
    }
    return res;
}




/*
 Given two strings containing digits, return the one which 
 represents the largest integer once the digits have been 
 sorted in non-increasing order.
 
 “245” -> 542
 “178” -> 871
 return 178
 */
string findLargest(string s1, string s2)
{
    vector<int> a1(10, 0);
    vector<int> a2(10, 0);
    for(int i = 0;i< s1.length();i++)
    {
        a1[s1[i] - '0'] ++;
    }
    
    for(int i = 0;i<s2.length();i++)
    {
        a2[s2[i] - '0'] ++;
    }
    
    int res = -1;
    for(int i = 0;i<10;i++)
    {
        res = a1[i] > a2[i]? 1 : 2;
    }
    if(res==1) return s1;
    else return s2;
}

/*
 2. We have N numbers 0 ~ N-1. A list contains k of the N numbers, e.g. 
 [1, 3, 4, 6, 7, 9], N = 10, K = 6. The list is sorted. The program returns 
 1 of N-K missing numbers in the list with probability 1/(N-K)
 */
int getRandom(vector<int> input, int n)
{
    int count = n-(int)input.size();
    int offset = rand() % count;
    for(int i = 0;i<(int)input.size();i++)
    {
        if(offset >=input[i]) offset++;
        else break;
    }
    return offset;
}


/*
 写一个小游戏。MxN 的格子上有一条蛇，蛇头可以向前，左，右移动，撞到自己身体任何部位或者撞到边界就算死。
 */





/*
 b) 有若干个盒子，每个盒子有length和width，不考虑高度。只要尺寸fit，大盒子就可以放小盒子，
 但是一层只能套一个，即便还有空余；但可以多层嵌套。求最小的面积放所有的盒子
 比如 7*7  5*5, 4*6, 3*3
 */


/*
 Q1：一个文件里存着代码和注释，注释在/××/中间，要求print所有line除了注释
 */
void printx(string s)
{
    bool isComments = false;
    int level = 0;
    for(int i = 0;i< s.length();i++)
    {
        if(!isComments)
        {
            if(s[i] != '/')
            {
                cout<<s[i];
            }
            else
            {
                if(i == s.length()-1 || s[i+1] != '*')
                {
                    cout<<s[i];
                }
                else
                {
                    i++;
                    isComments = true;
                    level = 1;
                }
            }
        }
        else
        {
            if(s[i] == '/' && i < s.length()-1 && s[i+1] == '*')
            {
                level++;
            }
            else
            {
                if(s[i] == '*' &&
                   (i<s.length() -1 && s[i+1] == '/'))
               {
                   i++;
                   level --;
                   if( level <=0)
                   {
                       isComments = false;
                   }
               }
            }
        }
    }
    
    if(isComments)
    {
        cout<< endl;
        cout<< "Wrong format, throw exception here"<<endl;
    }
}

/*
 2.convert binary tree to double linked list
 */

void conertBSTtoDoubleLinkedList(TreeNode* root)
{
    TreeNode* tmp = root;
    TreeNode* prev = NULL;
    while(tmp)
    {
        if(tmp->right && tmp->left)
        {
            TreeNode* right = tmp->right;
            TreeNode* left = tmp->left;
            while(left->right) left = left->right;
            left->right = right;
            tmp->right = tmp->left;
            tmp->left = NULL;
        }
        else if(tmp->left && !tmp->right)
        {
            tmp->right = tmp->left;
            tmp->left = NULL;
        }
        tmp->left = prev;
        prev = tmp;
        tmp = tmp->right;
    }
}

/*
1.string serialize & deserialize

serialize: 输入两个string，返回serialized stringdeserialize：输入serialized string，返回原来两个string
*/
char separator= '\y';
string serialize2str(string s1, string s2)
{
    string res = s1 + '\y' + s2;
    return res;
}

void deserialize(string s)
{
    string s1 = "";
    string s2 = "";
    int i = 0;
    bool isFirstStr = true;
    while(i<s.length())
    {
        
        if(s[i]!= '\y')
        {
            if(isFirstStr)
            {
                s1 += s[i];
            }
            else
            {
                s2 += s[i];
            }
        }
        else
        {
            isFirstStr = false;
        }
            
        i++;
    }
    cout<<s1<<endl<<s2<<endl;
    
}

/*
 2.integer divide without using divide/
 */
int divideAgain(int n1, int n2)
{
    long long res = 0;
    long long a = n1>=0? n1: -(long long)n1;
    long long b = n2>=0? n2: -(long long)n2;
    while(a >= b)
    {
        long long c = b;
        for(int i = 0;a>=c;i++, c<<=1)
        {
            a -= c;
            res += 1<<i;
        }
    }
    return ((n1^n2)>>31) ?-(int)res: (int)res;
}

int divideAccurate(int dividend, int divisor) {
    // 当 dividend = INT_MIN 时,-dividend 会溢出,所以用 long long
    long long a = dividend >= 0 ? dividend : -(long long)dividend;
    long long b = divisor >= 0 ? divisor : -(long long)divisor;
    // 当 dividend = INT_MIN 时,divisor = -1 时,结果会溢出,所以用 long long
    long long result = 0;
    while (a >= b) {
        long long c = b;
        for (int i = 0; a >= c; ++i, c <<= 1) {
            a -= c;
            result += 1 << i;
        }
    }
    return ((dividend^divisor) >> 31) ? (-(int)result) : ((int)result);
}


/*
 4.给一个string，比如UAXXBAUB，给一个pattern，比如AB，返回包含pattern的最短substring，结果是AUB
 */
string FindShortestPattern(string s, string p)
{
    int minLength = INT_MAX;
    string res;
    for(int i = 0;i<s.length()- p.length()+1;i++)
    {
        int k = i;
        for(int j = 0;j< p.length();j++)
        {
            while(s[k]!= p[j] && k<s.length())
            {
                k++;
            }
        }
        if(k<s.length() && k-i+1 < minLength)
        {
            minLength = k-i+1;
            res = s.substr(i, k-i+1);
        }
    }
    return res;
}


// isOneEditDistance 判断两个string是不是只差一个编辑距离

bool isOneEditDistanceIII(string s1, string s2)
{
    int l1 = (int)s1.length();
    int l2 = (int)s2.length();
    if(abs(l1-l2) >1)return false;
    if(l2> l1)
    {
        swap(s1, s2);
        swap(l1, l2);
    }
    int first = 0;
    int second = 0;
    int diff = 0;
    if( l1 == l2)
    {
        while(first < l1)
        {
            if(s1[first++] != s2[first++])
            {
                diff++;
            }
            if(diff >1) return false;
        }
    }
    else
    {
        while(second < l2)
        {
            if(s1[first] != s2[second])
            {
                diff++;
                first++;
            }
            else
            {
                first++;
                second++;
            }
            if(diff>1) return false;
        }
    }
    return true;
}


 // 1.    Find successor in BST
 // a. assume there is no dup.
TreeNode* findSuccessor(TreeNode* root, int target)
{
    TreeNode* res = NULL;
    while(root)
    {
        if(root->val <= target)
        {
            root = root->right;
        }
        else
        {
            res = root;
            root = root->left;
        }
    }
    return res;
}
/*
 implement memcpy.
 */
void mymemcpy2(void *dest, const void *source, size_t num) {
    int i = 0;
    // casting pointers
    char *dest8 = (char *)dest;
    char *source8 = (char *)source;
    //printf("Copying memory %d byte(s) at a time\n", sizeof(char));
    for (i = 0; i < num; i++) {
        // make sure destination doesnt overwrite source
        if (&dest8[i] == source8) {
            printf("destination array address overwrites source address\n");
            return;
        }
        dest8[i] = source8[i];
    }
    /* this is also a way to copy mem. it is a little bit neatter.
    while(num--)
    {
        *dest8++=*source8++;
    }
    */
}

/*
//write a function f(int x), so that f(x) returns true with x% probability。
//assme 0<=x<=100 and integer
 */
bool f(int x)
{
    return random()%100 <= x;
}

/*
//1. moving all 0s to the beginning of the array
*/
void moveZeroes(int a[], int n)
{
    int left = 0;
    int right = n-1;
    while(left<=right)
    {
        if(a[left] == 0) left++;
        else swap(a, left, right--);
    }
}

/*
 implement strStr()
 */
int strStrIII(char *haystack, char *needle) {
    int l1 = (int)strlen(haystack);
    int l2 = (int)strlen(needle);
    if(l1<l2) return -1;
    int next[l2];
    next[0] = -1;
    int j = next[0];
    for(int i = 1;i< l2;i++)
    {
        while(j!= -1 && needle[i] != needle[j+1])
            j = next[j];
        if(needle[j+1] == needle[i])
            j++;
        next[i] = j;
    }
    j = next[0];
    for(int i = 0;i< l1;i++)
    {
        while(j!= -1 && haystack[i] != needle[j+1])
            j = next[j];
        if(haystack[i] == needle[j+1])
            j++;
        if(j == l2-1)
        {
            return i - j;
        }
    }
    return -1;
}


/*
 1) 给个数组seq， 和一个total，找 if there is a contiguous sequence in seq
 which sums to total.都是正数， 第一次没注意contiguous，给了个back tracking的解法。然后说是
 contiguous， 给了个维护窗口的解法，不过犯了个小错误。时间过去了半小时。。。
 */
bool findContigious(int a[], int n, int target)
{
    // couple assumptions: the sum array -> tracker all values are unique.
    unordered_map<int, int> s;
    int tracker[n];
    int sum = 0;
    for(int i = 0;i< n;i++)
    {
        tracker[i] = sum+ a[i];
        s[tracker[i]]= i;
    }
    
    for(int i=0;i<n;i++)
    {
        if(s.find(tracker[i] + target) != s.end() &&
           s[tracker[i]+target] > i)
        {
            return true;
        }
    }
    return false;
}

/*
//先问怎么求submatrix的和
 */
int findSubMatrix(vector<vector<int>> input, int i1, int j1, int i2, int j2)
{
    int n = (int)input.size();
    int m = (int)input[0].size();
    //check i1, j1, i2, j2 must be smaller than n,m...
    vector<vector<int>> tracker(n, vector<int>(m,0));
    tracker[0][0] = input[0][0];
    for(int i = 1;i<m;i++)
    {
        tracker[0][i] = input[0][i] + tracker[0][i-1];
    }
    for(int i = 1;i<n;i++)
    {
        tracker[i][0] = input[i][0] + tracker[i-1][0];
    }
    for(int i = 1;i<n;i++)
    {
        for(int j = 1;j<m;j++)
        {
            tracker[i][j] = tracker[i-1][j] + tracker[i][j-1] - tracker[i-1][j-1] + input[i][j];
        }
    }
    
    if(i2<i1) swap(i1, i2);
    if(j2<j1) swap(j1, j2);
    
    int res = tracker[i2][j2] -
    (i1 == 0? 0 : tracker[i1-1][j2]) -
    (j1 == 0 ? 0 : tracker[i2][j1-1]) +
    (i1==0 || j1==0 ? 0:tracker[i1-1][j1-1]);

    return res;
}

//单链表k个k个分组，反转奇数组。比如 link = 0->1->2->3->4->5->6->7, k = 3返回 2->1->0->3->4->5->7->6
ListNode* reverseList(ListNode * head)
{
    ListNode* newHead = NULL;
    ListNode* tmp = head;
    while(head)
    {
        tmp = head;
        head = head->next;
        tmp->next = newHead;
        newHead = tmp;
    }
    return newHead;
}

ListNode* reverseLinkedList(ListNode* head, int k)
{
    if(k<2) return head;
    bool shouldReverse = true;
    ListNode* dummyHeader = new ListNode(-1);
    dummyHeader->next = head;
    
    ListNode* thisEnd;
    ListNode* secondStart;
    ListNode* secondEnd;
    ListNode* thirdStart;
    
    thisEnd = dummyHeader;
    while(thisEnd)
    {
        secondStart = thisEnd->next;
        ListNode* tmp = secondStart;
        for(int i = 0;i<k-1;i++)
        {
            if(tmp)
            {
                tmp = tmp->next;
            }else
            {
                break;
            }
        }
        secondEnd = tmp;
        thirdStart = tmp? tmp->next : NULL;

        if(shouldReverse)
        {
            if(secondEnd)
            {
                secondEnd->next = NULL;
            }
            ListNode* newhead = reverseList(secondStart);
            thisEnd->next = newhead;
            
            secondStart->next = thirdStart;
            thisEnd = secondStart;
        }
        else
        {
            thisEnd = secondEnd;
        }
        shouldReverse = !shouldReverse;
    }
    head = dummyHeader->next;
    delete dummyHeader;
    delete thisEnd;
    delete secondEnd;
    delete thirdStart;
    
    return head;
    
}

double mysqrt(double x)
{
    if(x<0) return -1;
    double left = x<1? x: 0;
    double right = x<1? 1.0: x;
    while(abs(right-left) > 0.0001)
    {
        double mid = left + (right-left)/2.0;
        if(mid*mid > x)
        {
            right = mid;
        }
        else
        {
            left = mid;
        }
    }
    return left;
}


/*
 wildcard matching.
 */
bool isMatch4(char* s, char* p)
{
    bool star = false;
    char *str, *ptr;
    for (str = s, ptr = p; *str != '\0'; str++, ptr++) {
        switch (*ptr) {
            case '?':
                break;
            case '*':
                star = true;
                s = str;
                p = ptr;
                while (*p == '*') p++; //skip continuous '*'
                if (*p == '\0') return true;
                str = s - 1;
                ptr = p - 1;
                break;
            default:
                if (*str != *ptr) {
                    if (!star) return false;
                    s++;
                    str = s - 1;
                    ptr = p - 1;
                }
        }
    }
    while (*ptr == '*') ptr++;
    return (*ptr == '\0');
}

bool isMatchRecursive(char* s, char* p)
{
    if(*p =='*')
    {
        while(*p=='*') p++;
        if(*p == '\0') return true;
        while(*s!='\0' && !isMatchRecursive(s, p)) s++;
        return *s != '\0';
    }
    else if(*p =='\0' || *s == '\0') return *p == *s;
    else if(*p =='?' || *p == *s) return isMatchRecursive(++s, ++p);
    else return false;
}

TreeNode* a;
bool started = false;
stack<TreeNode*> s;
bool backTrack =false;
void initializeTree(TreeNode* input)
{
    a = input;
    started= false;
}

bool haveNextTreeNode()
{
    if(started && s.empty() && a)
    {
        return false;
    }
    if(s.empty())
    {
        s.push(a);
    }
    return true;
}

TreeNode* getNextNode()
{
    if(!haveNextTreeNode())
        return NULL;
    TreeNode* res = NULL;
    started= true;
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        if(tmp->left && !backTrack)
        {
            s.push(tmp->left);
            continue;
        }
        res = tmp;
        s.pop();
        backTrack = true;
        if(tmp->right)
        {
            s.push(tmp->right);
            backTrack = false;
        }
        if(res)
        {
            return res;
        }
    }
    return NULL;
}


string longestCommonSubStr(string s1, string s2)
{
    TrieNode* root = new TrieNode();
    for(int i = 0;i<s1.length();i++)
    {
        string tmp = s1.substr(i, s1.length()-i);
        TrieNode* t = root;
        for(int j = 0;j<tmp.length();j++)
        {
            if(!t->val[tmp[j]]) t->val[tmp[j]] = new TrieNode();
            t = t->val[tmp[j]];
        }
    }
    
    int longest = INT_MIN;
    string res = "";
    for(int i = 0;i<s2.length();i++)
    {
        string tmp = s2.substr(i, s2.length()-1);
        TrieNode* t= root;
        int curLength = 0;
        for(int j = 0;j< tmp.length();j++)
        {
            if(t->val[tmp[j]])
            {
                curLength++;
                t = t->val[tmp[j]];
            }
            else
            {
                break;
            }
            
            if(longest < curLength)
            {
                longest = curLength;
                res = tmp.substr(0, longest);
            }

        }
    }
    return res;
}

vector<vector<int>> input;
vector<vector<int>>::iterator parentIterator;
vector<int>::iterator item;
int currenIndex = 0;
void initializerVector(vector<vector<int>> x)
{
    input = x;
    parentIterator = input.begin();
    item = (*parentIterator).begin();
}

bool hasNextIterator()
{
    if(item != (*parentIterator).end())
        return true;
    else
    {
        parentIterator++;
        while(parentIterator!= input.end() &&
              (*parentIterator).begin() == (*parentIterator).end())
        {
            parentIterator++;
        }
        if(parentIterator != input.end())
        {
            item = (*parentIterator).begin();
            return true;
        }
        else
        {
            return false;
        }
    }
}

vector<int>::iterator getNextIterator()
{
    if(hasNextIterator())
    {
        vector<int>::iterator res = item;
        item++;
        return res;
    }
    return (*parentIterator).end();
    /*
    if(item != (*parentIterator).end())
    {
        item++;
        return item;
    }
    else
    {
        parentIterator++;
        if(parentIterator != input.end())
        {
            parentIterator++;
            item = (*parentIterator).begin();
            return item;
        }
        else
        {
            return item;
        }
    }
    */
}

/*
从矩阵一段走对角线多少种走法，然后用组合数写了，然后加入障碍，然后用动态规划写了
 */
int findHowManyWays(int n, int m)
{
    //check n or m
    vector<vector<int>> tracker(n, vector<int>(m,1));
    for(int i = 1;i<n;i++)
    {
        for(int j = 1;j<m;j++)
        {
            tracker[i][j] = tracker[i-1][j] + tracker[i][j-1];
        }
    }
    return tracker[n-1][m-1];
}

int findHowManyWaysOneArray(int n, int m)
{
    vector<int> tracker (m, 1);
    for(int i = 1;i<n;i++)
    {
        for(int j = 1;j<m;j++)
        {
            tracker[j] = tracker[j]+ tracker[j-1];
        }
    }
    return tracker[m-1];
}

/*
 二面也很简单，问了一道一堆数据，从一个开始之后是另外一个类型，然后binary search轻松搞定
 */
int findChange(vector<int> input)
{
    int left = 0;
    int right = (int)input.size()-1;
    if( input[0] == input[input.size()-1]) return 0;
    while(left<=right)
    {
        int mid = left + (right-left)/2;
        if(mid > 0 && input[mid-1] != input[mid])
        {
            return mid;
        }
        else if(mid < input.size()-1 && input[mid] != input[mid+1])
        {
            return mid+1;
        }
        else if(input[left] == input[mid])
        {
            left = mid+1;
        }
        else
        {
            right = mid-1;
        }
    }
    return -1;
}

/*
 一个美国小哥，树的遍历，当然顺序是反向的，然后要求用另外一种实现（空间换时间），还有直接打印树中第n个节点
 */
int printNth(TreeNode* root, int n)
{
    stack<TreeNode*> s;
    s.push(root);
    int count = 0;
    bool backTrack= false;
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        if(tmp->left && !backTrack)
        {
            s.push(tmp->left);
            continue;
        }
        
        s.pop();
        count++;
        if(count == n) return tmp->val;
        backTrack = true;
        if(tmp->right)
        {
            s.push(tmp->right);
            backTrack = false;
        }
            
        
    }
    return -1;
}

/*
还是一个美国小哥，打印所有树的path 使得sum为一个给定的数.
*/


int printNthWithNumber(TreeNode* root, int n)
{
    // that means n is larger than the count of the all nodes in the tree
    if(n>root->val) return -1;
    while(root)
    {
        int leftcount = (root->left)? root->left->val: 0;
        if(n > leftcount+1)
        {
            root = root->right;
            n = n-leftcount-1;
        }
        else if(n<leftcount+1)
        {
            root = root->left;
        }
        else
        {
            return root->val;
        }
    }
    return -1;
}

/*
 打印所有树的path 使得sum为一个给定的数.（假设path是from root 到leaf 的path）
 */
void printPath(stack<TreeNode*> input, int target)
{
    stack<TreeNode*> s;
    int sum = 0;
    while(!input.empty())
    {
        sum += input.top()->val;
        s.push(input.top());
        input.pop();
    }
    
    if(target == sum)
    {
        while(!s.empty())
        {
            cout<< s.top()->val<<"->";
            s.pop();
        }
        cout<<endl;
    }
    
}

void FindAllPath(TreeNode* root, int target)
{
    if(!root) return;
    stack<TreeNode*> s;
    TreeNode* prev = NULL;
    s.push(root);
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        if(tmp->left &&
           tmp->left != prev &&
           (tmp->right && tmp->right != prev))
        {
            s.push(tmp->left);
        }
        else if( tmp->right && tmp->right != prev)
        {
            s.push(tmp->right);
        }
        else
        {
            if(!tmp->left && !tmp->right)
            {
                printPath(s, target);
            }
            s.pop();
            prev = tmp;
        }
    }
}

/*
 roman number to integer
 */

int ConvertRomanToInteger(string roman)
{
    unordered_map<char, int> tracker;
    tracker['I'] = 1;
    tracker['V'] = 5;
    tracker['X'] = 10;
    tracker['L'] = 50;
    tracker['C'] = 100;
    tracker['D'] = 500;
    tracker['M'] = 1000;
    
    int res = 0;
    for(int i = 0;i< roman.length();i++)
    {
        if(i>0 && tracker[roman[i-1]] < tracker[roman[i]])
        {
            res += (tracker[roman[i]] - 2 * tracker[roman[i-1]]);
        }
        else
        {
            res += tracker[roman[i]];
        }
    }
    return res;
}

/*
 integer to roman
 */
void convert(int k, char romanchars[], string& res)
{
    if(k<=0)
        return;
    else if(k<=3)
        res.append(k, romanchars[0]);
    else if(k==4)
    {
        res.append(1, romanchars[0]);
        res.append(1, romanchars[1]);
    }
    else if(k<=8)
    {
        res.append(1, romanchars[1]);
        res.append(k-5, romanchars[0]);
    }
    else if(k==9)
    {
        res.append(1, romanchars[0]);
        res.append(1, romanchars[2]);
    }
    
}

string ConvertIntToRoman(int n)
{
    char roman[] = {'I', 'V', 'X', 'L', 'C', 'D', 'M'};
    string res;
    int factor = 1000;
    int i = 6;
    while(n!= 0)
    {
        convert(n/factor, roman+i, res );
        i -= 2;
        n %= factor;
        factor /=10;
    }
    return res;
}

/*
 就是给你一个string判断是否是palindrome，忽略大小写，忽略非字母符号
 */
bool isNumAlpha(char c)
{
    if(( c>='0' && c <='0') ||
       (c >= 'a' && c<='z'))
    {
        return true;
    }
    return false;
}

bool isPalindromeX(string s)
{
    int left = 0;
    int right = (int)s.length()-1;
    while(left <= right)
    {
        while(!isNumAlph(tolower(s[left]))) left++;
        while(!isNumAlph(tolower(s[right]))) right --;
        if(tolower(s[left]) == tolower(s[right]))
        {
            left++;
            right--;
        }
        else
        {
            return false;
        }
    }
    return true;
}

int strStr4(string haystack, string needle)
{
    int l1 = (int)haystack.length();
    int l2 = (int)needle.length();
    if(l2>l1) return -1;
    int next[l2];
    next[0] =-1;
    int j = -1;
    for(int i = 1;i<l2;i++)
    {
        while(j != -1 && needle[j+1] != needle[i])
            j = next[j];
        if(needle[j+1]== needle[i])
            j++;
        next[i] = j;
    }
    j = -1;
    for(int i = 0;i<l1;i++)
    {
        while(j != -1 && haystack[i] != needle[j+1])
            j = next[j];
        if(haystack[i] == needle[j+1])
            j++;
        if(j == l2-1)
            return i - j;
    }
    return -1;
}

bool threeSum(vector<int> input, int target)
{
    sort(input.begin(), input.end());
    for(int i = 0;i< input.size()-2;i++)
    {
        int left = i+1;
        int right = (int)input.size() -1;
        while(left<right)
        {
            int sum = input[left] + input[right] + input[i];
            if(sum == target)
            {
                return true;
            }
            if(sum < target)
            {
                left++;
            }
            else
            {
                right--;
            }
        }
    }
    return false;
    
}

bool threeSumWithHasTable(vector<int> input,  int target)
{
    unordered_map<int, vector<pair<int, int>>> tracker;
    int n = (int)input.size();
    for(int i = 0;i<n;i++)
    {
        for(int j = i+1;j<n;i++)
        {
            tracker[input[i]+input[j]].push_back(make_pair(i, j));
        }
    }
    
    for(int i = 0;i<n;i++)
    {
        if(tracker.find(target-input[i]) != tracker.end())
        {
            for(auto item : tracker[target-input[i]])
            {
                if(i != item.first && i!= item.second)
                {
                    return true;
                }
            }
        }
    }
    return false;
    
}

/*
 Letter Combinations of a Phone Number，
 */
void letterCombinationWorkeri( string digits, string& tmp, vector<string>& res, unordered_map<char, string> map)
{
    if(tmp.length() == digits.length())
    {
        res.push_back(tmp);
    }
    else
    {
        string x = map[digits[tmp.length()]];
        for(int i = 0;i<(int)x.length();i++)
        {
            tmp.push_back(x[i]);
            letterCombinationWorkeri(digits, tmp, res, map);
            tmp.pop_back();
        }
    }
}

vector<string> letterCombinationsi(string digits) {
    vector<string> res;
    
    unordered_map<char, string> map;
    map['2']= "abc";
    map['3'] = "def";
    map['4'] = "ghi";
    map['5'] = "jkl";
    map['6'] = "mno";
    map['7'] = "pqrs";
    map['8'] = "tuv";
    map['9'] = "wxyz";
    string tmp = "";
    letterCombinationWorkeri(digits, tmp, res, map);
    
    return res;
}

/*
 Anagram groups
 */
unordered_map<string, vector<string>> AnagramGroups(vector<string> input)
{
    unordered_map<string, vector<string>> res;
    for(auto item: input)
    {
        string s = item;
        sort(s.begin(), s.end());
        res[s].push_back(item);
    }
    return res;
    
}

/*
 2. Decode ways
 */
int DecodeWays(string s)
{
    int prev = 0;
    int cur = 1;
    for(int i = 0;i<s.length();i++)
    {
        if(s[i] == '0') cur = 0;
        if(i == 0 || !(s[i-1] == '1' ||(s[i-1]=='2' && s[i] <= '6'))) prev = 0;
        int tmp = cur;
        cur = cur + prev;
        prev = tmp;
    }
    return cur;
}

/*
2.Find minimum number in a rotated sorted array (当时这个题还没在
leetcode里，所以写得代码有些繁琐，估计因为这个要再电面一轮）
*/
int findMinX(vector<int> num)
{
    int size = (int)num.size();
    if(size == 0) return 0;
    int left = 0;
    int right =  size-1;
    while(left <=right)
    {
        int mid = left + (right-left)/2;
        
        if((num[mid] < num[mid-1] || mid == 0) &&
           (num[mid] < num[mid+1] || mid == size -1))
        {
            return num[mid];
        }
        else if(num[left] < num[mid] && num[mid] < num[right])
        {
            return num[left];
        }
        else if(num[left] > num[mid])
        {
            right = mid - 1;
        }
        else if(num[mid] > num[right])
        {
            left = mid + 1;
        }
    }
    return num[left-1];
}

/*
 
 1.    Insert a node into a sorted circular linked list ( all next element is
 larger except for the last one), the given head can point to any node
 
 1 -> 3 -> 5 ->7
 ^             |
 |             |
 |  _  _  _  _ |
 
 如果node的值是2，则插入1和3之间；如果node的值是8或者0，插入7和1之间。
 
 要考虑node值重复的情况，虽然结果一样，但要和面试官讨论新的节点插入的位置，可
 能插入在最开始或最后，我不记得了。
 
 例如插入3, 结果是1->3->3'->5->7或者1->3'->3->5->7
*/
void insertCircularList(ListNode* head, int target)
{
    if(!head) return;
    while(head)
    {
        if(head->val == target)
        {
            ListNode* tmp = new ListNode(target);
            tmp->next = head->next;
            head->next = tmp;
            return;
        }
        else if(head->val < target)
        {
            if(head->next->val >=target)
            {
                ListNode* tmp = new ListNode(target);
                tmp->next = head->next;
                head->next = tmp;
                return;
            }
            else
            {
                if(head->next->val >= head->val)
                {
                    head = head->next;
                }
                else
                {
                    ListNode* tmp = new ListNode(target);
                    tmp->next = head->next;
                    head->next = tmp;
                    return;
                }
            }
        }
        else
        {
            head = head->next;
        }
    }
}

/*
 2.    Clone graph(leetcode)
 */

vector<GraphNode*> clone(vector<GraphNode*> input)
{
    unordered_map<GraphNode*, GraphNode*> tracker;
    for(auto item: input)
    {
        tracker[item] = new GraphNode(item->val);
    }
    
    for(auto item: input)
    {
        for(auto nei: item->neighbors)
        {
            tracker[item]->neighbors.push_back(tracker[nei]);
        }
    }
    
    vector<GraphNode*> res;
    for(auto item: input)
    {
        res.push_back(tracker[item]);
    }
    return res;
    
}

/*
 regular expression DP
 */
bool canMatch(char a, char b)
{
    return (a == b || b == '.');
}

bool isMatchDP(const char *s, const char *p) {
    int lS = (int)strlen(s);
    int lP = (int)strlen(p);
    vector<vector<bool> > F(lS + 1, vector<bool>(lP + 1));
    F[0][0] = true;
    for (int i = 0; i <= lS; i++)
    {
        for (int j = 1; j <= lP; j++)
        {
            if(i>0)
                // matches one character, index of both string and pattern move forward by 1
                if (F[i-1][j-1] && canMatch(s[i-1], p[j-1]))
                {
                    F[i][j] = true;
                    continue;
                }
            if (i > 0 && j > 1)
                // matches the situation when the next char in the string is the same as the char before the '*' in the pattern
                // pre-match: string xxxCyyy, pattern mmmC*nnn, xxx matches mmmC*
                // current match: xxxC matches mmmC*
                if (F[i-1][j] && canMatch(s[i-1], p[j-2]) && p[j-1] == '*')
                {
                    F[i][j] = true;
                    continue;
                }
            if (j > 1)
                // matches the situation when the next two chars in the pattern is in the form of C*
                // pre-match: string xxxyyy, pattern mmmC*nnn, xxx matches mmm
                // current match: xxx matches mmmC*
                if (F[i][j-2] && p[j-1] == '*')
                {
                    F[i][j] = true;
                    continue;
                }
        }
    }
    return F[lS][lP];
}

/*
 告诉面试官有DP 解法, 告知我写递归. 我也写出个递归解法, 并且用hashmap 保存出
 现过的子串来优化
 */
bool wordBreak(string s, unordered_set<string> dict)
{
    int n = (int)s.length();
    vector<bool> tracker(n);
    for(int i=0;i<n;i++)
    {
        for(int j =-1;j<i;j++)
        {
            bool prev = j<0? true: tracker[j];
            if(dict.find(s.substr(j+1, i-j)) != dict.end() && prev)
            {
                tracker[i] = true;
                break;
            }
        }
    }
    return tracker[n-1];
}


bool wordBreakRecursive(string s, unordered_set<string> dict)
{
    int n = (int)s.length();
    for(int i = 1;i<=n;i++)
    {
        bool prev = dict.find(s.substr(0, i)) != dict.end();
        if(prev && wordBreakRecursive(s.substr(i, n-i), dict))
        {
            return true;
        }
    }
    return false;
}

/*// Reverse a Singly Linked List
 // Example Input: A -> B -> C
 // Example Output: C -> B -> A
 他先让说思路，然后问时间和空间复杂度，然后写代码。说思路说了半天，这种list的
 题，就是画图，英语不好说起来真费劲。。。这道题应该是Leetcode上一道的一个小部
 分，所以很快就写完了。
 */

ListNode* reversex(ListNode* head)
{
    ListNode* newHead = NULL;
    ListNode* tmp;
    while(head)
    {
        tmp = head;
        head = head->next;
        tmp->next = newHead;
        newHead = tmp;
    }
    return newHead;
}

/*
 2.第二题直接copy 题目，感觉跟leetcode上面的interval那题很相似，简单一点点。
 // Given a array of pairs where each pair contains the start and end time of a meeting (as in int),
 // Determine if a single person can attend all the meetings
 // Input array(pair(1,4), pair(4, 5), pair(3,4), pair(2,3))
 // Output: False
 同样的思路+复杂度。 同样是变种题，还变简单了，很快写完。（主要是考比较器的
 override吧）
 */

bool compareX(EndPoint ep1, EndPoint ep2)
{
    if(ep1.val < ep2.val) return true;
    if(ep1.val == ep2.val && ep1.count < ep2.count) return true;
    return false;
}

bool canAttendAllmeeting(vector<Interval> input)
{
    vector<EndPoint> tracker;
    for(auto item: input)
    {
        EndPoint ep1(item.start, 1);
        EndPoint ep2(item.end, 1);
        tracker.push_back(ep1);
        tracker.push_back(ep2);
    }
    
    sort(tracker.begin(), tracker.end(), compareX);
    int maxCount = 0;
    int cur = 0;
    for(auto ep: tracker)
    {
        cur+=ep.count;
        maxCount = max(cur, maxCount);
    }
    return maxCount <=1;
    
}

/*
 3. follow up第二题
 // determine the minimum number of meeting rooms needed to hold all the meetings.
 // Input array(pair(1, 4), pair(2,3), pair(3,4), pair(4,5))
 // Output: 2
 */
int FindHowMayMeetingRoomNeeded(vector<Interval> input)
{
    vector<EndPoint> tracker;
    for(auto item: input)
    {
        EndPoint ep1(item.start, 1);
        EndPoint ep2(item.end, 1);
        tracker.push_back(ep1);
        tracker.push_back(ep2);
    }
    
    sort(tracker.begin(), tracker.end(), compareX);
    int maxCount = 0;
    int cur = 0;
    for(auto ep: tracker)
    {
        cur+=ep.count;
        maxCount = max(cur, maxCount);
    }
    return maxCount ;
    
}

/*
 Given an integer array, adjust each integers so that the difference of every adjcent integers are not greater than a given number target.If the array before adjustment is A, the array after adjustment is B, you should minimize the sum of |A[i]-B[i]|
 注意
 You can assume each number in the array is a positive integer and not greater than 100
 样例
 Given [1,4,2,3] and target=1, one of the solutions is [2,3,2,3], the adjustment cost is 2 and it's minimal. Return 2.
 */
int findMinAdjust(vector<int> input, int target)
{
    int n = (int)input.size();
    vector<vector<int>> tracker(n, vector<int>(101, INT_MIN));
    for(int i = 0;i< n;i++)
    {
        for(int j = 1;j<100;j++)
        {
            int minPrevadjust = INT_MAX;
            if(i >0)
            {
                for(int k = 0;k<=target;k++)
                {
                    if(j-k>=0) minPrevadjust = min(minPrevadjust, tracker[i-1][j-k]);
                    if(j+k<=100) minPrevadjust = min(minPrevadjust, tracker[i-1][j+k]);
                }
            }
            else
            {
                minPrevadjust = 0;
            }
            tracker[i][j] = abs(j-input[i]) + minPrevadjust;
        }
    }
    
    int res = INT_MAX;
    for(int i = 0;i< 101;i++)
    {
        res = min(tracker[n-1][i], res);
    }
    return res;
}

/*
 题目是：Given a set of n jobs with [start time, end time, cost] find a
 subset so that no 2 jobs overlap and the cost is maximum.
 
 brute force 肯定是不行的了。我知道一个O(N^2)的解答，请问有更好的答案吗？谢谢
 了！
 */
bool sortInterval(workItem t1, workItem t2)
{
    if(t1.start < t2.start)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int findMaxCost(vector<workItem> input)
{
    sort(input.begin(), input.end(), sortInterval);
    vector<int> cost;
    
    int n = (int)input.size();
    for(int i = 0;i< n;i++)
    {
        int maxCost = INT_MIN;
        for(int j= -1;j<i;j++)
        {
            if(j < 0) maxCost = input[i].val;
            else
            {
                if(input[i].start > input[j].end)
                {
                    maxCost = max(maxCost, cost[j]+ input[i].val);
                }
            }
        }
        cost.push_back(maxCost);
    }
    return cost[n-1];
    
}

// the difference between this method and prev method
// is that in this method I use one more item in the cost
// vector, so j's start indext can be 1.
// but I have to shift the cost index. it seem like
// there is no easy way to handle this.
int findMaxCostII(vector<workItem> input)
{
    sort(input.begin(), input.end(), sortInterval);
    
    int n = (int)input.size();
    vector<int> cost(n+1, 0);
    for(int i = 0;i< n;i++)
    {
        int maxCost = INT_MIN;
        for(int j= 0;j<=i;j++)
        {
            if(j < 0) maxCost = input[i].val;
            maxCost = max(maxCost, cost[j]+ input[i].val);
        }
        cost[i+1] = maxCost;
    }
    return cost[n];
    
}



/*
 minimum window substring
 */


string FindMinwindowsSubStr(string s, string p)
{
    int ptracker[256];
    int stracker[256];
    
    fill_n(&ptracker[0], 256, 0);
    fill_n(&stracker[0], 256, 0);
    
    for(int i = 0;i<p.length();i++)
    {
        ptracker[p[i]]++;
    }
    
    int head = 0;
    int start = -1;
    int length = INT_MAX;
    int currentNum = 0;
    
    for(int i = 0;i< s.length();i++)
    {
        if(ptracker[s[i]] >0)
        {
            if(stracker[s[i]] < ptracker[s[i]])
            {
                currentNum++;
            }
            stracker[s[i]]++;
        }
        
        if(currentNum == p.length())
        {
            while(head < i && (stracker[s[head]] > ptracker[s[head]] ||
                               ptracker[s[head]] == 0))
            {
                stracker[s[head]]--;
                head++;
            }
            if((i-head+1) < length)
            {
                start = head;
                length = i-head+1;
            }
        }
        
    }
    
    return start ==-1? "": s.substr(start, length);
}


/*{ "face", "ball", "apple", "art", "ah" }
 "htarfbp..."
 根据下面的string去给上面list words排序。
 就是平常我们按abcd。。。排，这次按string里的letter顺序排
 */
unordered_map<char, int> trackerx;
bool sortString(string s1, string s2)
{
    int l1 = (int)s1.length();
    int l2 = (int)s2.length();
    for(int i =0;i<l1 && i<l2;i++)
    {
        if(s1[i] == s2[i])
        {
            continue;
        }
        else
        {
            if(trackerx[s1[i]]<trackerx[s2[i]])
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    if( l1<l2) return true;
    else return false;
}

vector<string> mySort(string sequence, vector<string> strings)
{

    for(int i = 0;i< sequence.length();i++)
    {
        trackerx[sequence[i]] = i;
    }
    sort(strings.begin(), strings.end(), sortString);
    return strings;
}

/*
 第二题，将1->2->3->4->5->6->7 变成 1->7->2->6->3->5->4，不能用额外空间
 */
ListNode* shuffleList(ListNode* head)
{
    ListNode* fast = head;
    ListNode* slow = head;
    while(fast &&
          fast->next)
    {
        fast = fast->next->next;
        slow = slow->next;
    }
    
    ListNode* secondHead = slow->next;
    slow->next = NULL;
    
    ListNode* newSecondHead = NULL;
    ListNode* tmp= NULL;
    while(secondHead)
    {
        tmp = secondHead;
        secondHead = secondHead->next;
        tmp->next = newSecondHead;
        newSecondHead = tmp;
    }
    
    ListNode* dummyHead = new ListNode(-1);
    ListNode* tail = dummyHead;
    while(head && newSecondHead)
    {
        ListNode* tmp;
        tmp = head;
        head = head->next;
        
        tail->next = tmp;
        tail = tmp;
        
        tmp = newSecondHead;
        newSecondHead->next = newSecondHead;
        
        tail->next = tmp;
        tail = tmp;
    }
    
    if(head) tail->next = head;
    if(secondHead) tail->next = head;
    return dummyHead->next;
    
}

/*
 2.    Design a data structure supporting two operations
 1）    void addWord(string)
 2）    bool search(string)
 
 search(string) can search word and regular expression ( only consider “.”,
 which means any one character)
 
 例如
 addWord("rat")
 addWord("cat")
 addWord("bat")
 search("dat") -> false
 search("bat") -> true
 search(".at") -> true
 search("r.t") -> true
 */
unordered_set<string> containerx;
void addWord(string s)
{
    containerx.insert(s);
}


bool searchWordWorker(string tmp, vector<int> index, int i)
{
    if(i == index.size())
    {
        return containerx.find(tmp) != containerx.end();
    }
    else
    {
        for(char c = 'a';c<='z';c++ )
        {
            tmp[index[i]] = c;
            
            if(searchWordWorker(tmp, index, i+1))
            {
                return true;
            }
        }
    }
    return false;
}

bool searchWord(string s)
{
    vector<int> index;
    for(int i = 0;i< s.length();i++)
    {
        if(s[i] == '.')
        {
            index.push_back(i);
        }
    }
    return searchWordWorker(s, index, 0);
}

// use Trie to do this code.
TrieNode* host;
int SIZE = 256;
void addWordTrie(string s)
{
    if(!host) host = new TrieNode();
    TrieNode* tmp = host;
    for(int i = 0;i<s.length();i++)
    {
        if(!tmp->val[s[i]])
        {
            tmp->val[s[i]] = new TrieNode();
        }
        if(i == s.length()-1)
        {
            tmp->end[s[i]] = true;
        }
        tmp = tmp->val[s[i]];
    }
}


bool searchWordTrie(string s)
{
    queue<TrieNode*> q;
    TrieNode* tmp = host;
    q.push(tmp);
    q.push(NULL);
    int i = 0;
    while(!q.empty())
    {
        tmp = q.front();
        q.pop();
        // taking care of the end of a level
        if(!tmp)
        {
            i++;
            q.push(NULL);
            continue;
        }
        
        if(s[i] != '.')
        {
            if(tmp->val[s[i]])
            {
                if(tmp->end[s[i]] && i == s.length()-1) return true;
                else q.push(tmp->val[s[i]]);
            }
            else
            {
                return false;
            }
        }
        else
        {
            for(int j = 0;j< SIZE;j++)
            {
                if(tmp->val[j])
                {
                    if((i == s.length()-1) && tmp->end[j]) return true;
                    else q.push(tmp->val[j]);
                }
            }
        }
    }
    return false;
}

bool searchWordTrieDFSWorker(int level, string s, TrieNode* root)
{
    if(level == s.length()-1)
    {
        if(s[level] == '.')
        {
            for(int i = 0;i<SIZE;i++)
            {
                if(root->val[i] && root->end[i]) return true;
            }
            return false;
        }
        else
        {
            if(root->val[s[level]] && root->end[s[level]])
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        if(s[level] == '.')
        {
            for(int i = 0;i<SIZE;i++)
            {
                if(root->val[i])
                {
                    if(searchWordTrieDFSWorker(level+1, s, root->val[i]))
                    {
                        return true;
                    }
                }
            }
            return false;
        }
        else
        {
            if(root->val[s[level]])
            {
                return searchWordTrieDFSWorker(level+1, s, root->val[s[level]]);
            }
            else
            {
                return false;
            }
        }
    }
}

bool searchWordTrieDFS(string s)
{
    
    return searchWordTrieDFSWorker(0, s, host);
}

/*
 phone keyboard letter combination iterative
 */
vector<string> letterCombinationsx (const string &digits) {
    const vector<string> keyboard { " ", "", "abc", "def", // '0','1','2',...
        "ghi", "jkl", "mno", "pqrs", "tuv", "wxyz" };
    vector<string> res = {""};

    for(auto d : digits)
    {
        int n = (int)res.size();
        int m = (int)keyboard[d-'0'].length();
        if(m>0)
        {
            res.resize(n*m);
        }
        
        for(int i = 1;i< m;i++)
        {
            copy(res.begin(), res.begin()+n, res.begin()+n*i);
        }
        
        for(int i=0;i<m;i++)
        {
            for(int j = i*n;j<i*n+n;j++)
            {
                res[j] += keyboard[d-'0'][i];
            }
        }
    }
    return res;
}

/*
 1， 美国人， 给一个词，判断是不是Palindrome, 然后扩展问，给一个字典，找出所有对 单词，
 这两个单词可以组成一个palindom, 然后有问，可以组合任意个单词，怎么找到最长的可能的palindom
 */
bool isPalindrome(string s, int offset)
{
    if( offset >= s.length()-1) return true;
    int left = offset;
    int right = (int)s.length()-1;
    while(left<=right)
    {
        if(s[left]== s[right])
        {
            left--;right++;
        }
        else
        {
            return false;
        }
    }
    return true;
}
vector<pair<string, string>> FindAllPalindromePairs(vector<string> dict)
{
    for(auto s : dict)
    {
        addWordTrie(s);
    }
    vector<pair<string, string>> res;
    for(auto s: dict)
    {
        TrieNode* tmp = host;
        reverse(s.begin(), s.end());
        int index = 0;
        while(tmp &&
              tmp->val[s[index]] &&
              index < s.length())
        {
            if(tmp->end[s[index]])
            {
                if(isPalindrome(s, index+1))
                {
                    string tmp = s.substr(0, index+1);
                    string tmp2 = s;
                    reverse(tmp2.begin(), tmp2.end());
                    res.push_back(make_pair(tmp, tmp2));
                }
            }
            tmp = tmp->val[s[index]];
            index++;
        }
    }
    
    delete host;
    host = new TrieNode();
    for(auto s : dict)
    {
        reverse(s.begin(), s.end());
        addWordTrie(s);
    }
    
    for(auto s: dict)
    {
        TrieNode* tmp = host;
        int index = 0;
        while(tmp &&
              tmp->val[s[index]] &&
              index < s.length())
        {
            if(tmp->end[s[index]])
            {
                if(isPalindrome(s, index+1))
                {
                    string tmp1 = s.substr(0, index+1);
                    reverse(tmp1.begin(), tmp1.end());
                    res.push_back(make_pair(tmp1, s));
                }
            }
            tmp = tmp->val[s[index]];
            index++;
        }
    }
    return res;
}

/*
 第一题，给一个字符数组，要求将其中的'a'加倍，'b'删除，其他字符保持不变。要求
 inplace，线性复杂度。这一题做的很顺利。面试官说good enough
 */
vector<char> modify(vector<char> input)
{
    int a = 0;
    int b = 0;
    for(int i = 0;i<input.size();i++)
    {
        if(input[i] == 'a') a++;
        else if(input[i] == 'b') b++;
    }
    int originalsize = (int)input.size();
    int n = (int)input.size() + a -b;

    int left = 0;
    int runner = 0;
    while(runner < originalsize)
    {
        if(input[runner] != 'b')
        {
            input[left++] = input[runner];
        }
        runner++;
    }
    
    input.resize(n);
    int right = n-1;
    runner = left-1;
    while(runner>=0)
    {
        if(input[runner] == 'a')
        {
            input[right--] = 'a';
            input[right--] = 'a';
        }
        else if(input[runner] != 'b')
        {
            input[right--] = input[runner];
        }
        runner --;
    }
    return input;
}

/*
 1.几个数字array，像这样的：
 1
 11
 21
 1211
 111221
 给n，返回第n行的结果。第二行返回前面一行每个number的count。我用的recursive方
 法。不知道是不是最优的。
 */
string countAndSayII(int n)
{
    string s = "1";
    if(n == 1) return s;
    for(int i = 2;i<=n;i++)
    {
        string tmp = "";
        int count=0;
        char c = 'a';
        for(int j = 0;j<s.length();j++)
        {
            if(s[j] != c)
            {
                if(count >0)
                {
                    tmp+=to_string(count);
                    tmp+=c;
                }
                count = 1;
                c = s[j];
            }
            else
            {
                count++;
            }
        }
        tmp+=to_string(count);
        tmp+=c;
        s= tmp;
    }
    return s;
}

queue<int> thequeue;
mutex mt;
int size = 0;
int maxSize = 10;
condition_variable fullCondition;
condition_variable emptyCondition;

void producer(int n)
{
    unique_lock<mutex> locker(mt);
    while(size >= maxSize)
    {
        fullCondition.wait(locker);
    }
    
    thequeue.push(n);
    cout<<"produced "<<n<<endl;
    size++;
    emptyCondition.notify_all();
    locker.unlock();
    this_thread::sleep_for(chrono::milliseconds(100));

    
}

void consume()
{
    
    unique_lock<mutex> locker(mt);
    while(size<=0)
    {
        emptyCondition.wait(locker);
    }
    int tmp = thequeue.front();
    thequeue.pop();
    size--;
    fullCondition.notify_all();
    cout<<"consume "<<tmp<<endl;
    locker.unlock();
    this_thread::sleep_for(chrono::milliseconds(105));

}

void func1()
{
    for(int i = 0;i<100;i++)
    {
        producer(i);
    }
}

void func2()
{
    for(int i = 0;i<100;i++)
    {
        consume();
    }
}

/*
fabanacci，期待o(lgn)解法,但O(n)也行
*/

vector<vector<int>> multi(vector<vector<int>> m1, vector<vector<int>> m2)
{
    vector<vector<int>> res(2, vector<int>(2,0));

    res[0][0] = m1[0][0]* m2[0][0] + m1[0][1]*m2[1][0];
    res[0][1] = m1[0][0]* m2[0][1] + m1[0][1]*m2[1][1];
    res[1][0] = m1[1][0]* m2[0][0] + m1[1][1]*m2[1][0];
    res[1][1] = m1[1][0]* m2[0][1] + m1[1][1]*m2[1][1];
    
    return res;
}

vector<int> multi(vector<vector<int>> m1, vector<int> m2)
{
    vector<int> res(2);
    res[0] = m1[0][0] * m2[0] + m1[0][1] * m2[1];
    res[1] = m1[1][0] * m2[0] * m1[1][1] * m2[1];
    return res;
}

vector<vector<int>> power(vector<vector<int>> matrix, int n)
{
    if(n ==1) return matrix;
    vector<vector<int>> res = power(matrix, n/2);
    res = multi(res, res);
    if((n%2) ==1)
    {
        res = multi(res, matrix);
    }
    return res;
}

int findFibnacci(int n)
{
    if(n<0) return -1;
    if(n<2) return 1;
    vector<int> fn ={1, 1};
    vector<vector<int>> matrix(2, vector<int>(2, 1));
    matrix[1][1] = 0;
    
    vector<vector<int>> res = power(matrix, n-1);
    vector<int> fb = multi(res, fn);
    return fb[0];
}
/*
 generate all possible paretheses,
 */
//assum n mean pair.
void parenthesesPair(vector<string>& res, string tmp, int left, int right, int n )
{
    if(left == n && right ==n)
    {
        res.push_back(tmp);
    }
    else
    {
        if(left>right)
        {
            parenthesesPair(res, tmp+')', left, right+1, n);
        }
        if(left < n)
        {
            parenthesesPair(res, tmp+'(', left+1, right, n);
        }
    }
}



vector<string> generateAllParentheses(int n)
{
    vector<string> res;
    parenthesesPair(res, "", 0, 0, n);
    return res;
}

/*
 3）divide and mod，但不能用/或者%，基本也是leetcode原题了
 */

int divideIIII(int n1, int n2)
{
    long long a = abs((long long)n1);
    long long b = abs((long long)n2);
    int res = 0;
    while(a>=b)
    {
        long long c = b;
        for(int i = 0;a>=c;i++, c<<=1)
        {
            a -= c;
            res += 1<<i;
        }
    }
    return (n1^n2)>>31? -res: res;
}

/*
 Given a sequence of distinct integers, your program must remove as few
 elements as possible in order for the elements which are not removed to
 appear in ascending order.  If there is more than one way to do this, your
 program must print one solution, then print the number of all solutions.
 
 Example.
 
 Given   1 2 3 8 10 5 6 7 12 9 11 4 0
 Remove        8 10       12      4 0
 Remain  1 2 3      5 6 7    9 11       (ascending)
 */

int MakeAscending(vector<int> input)
{
    int n = (int)input.size();
    vector<int> tracker(n, 0);
    tracker[0] = 1;
    
    for(int i= 1;i<n;i++)
    {
        int maxCount = 0;
        
        for(int j = i-1;j>=0;j--)
        {
            if(input[j] < input[i])
            {
                maxCount = max(maxCount, tracker[j]+1);
            }
        }
        tracker[i] = maxCount;

    }
    return n-tracker[n-1];
    
}

/*
 跳河问题。给一个0/1数组R代表一条河，0代表水，1代表石头。起始位置R[0]等于1，
 初速度为1. 每一步可以选择以当前速度移动，或者当前速度加1再移动。只能停留在石
 头上。问最少几步可以跳完整条河流。
 
 给定数组为R=[1,1,1,0,1,1,0,0]，最少3步能过河：
 第一步先提速到2，再跳到R[2]；
 第二步先提速到3，再跳到R[5]；
 第三步保持速度3，跳出数组范围，成功过河。
 */
int FrogCrossRiver(string river) {
    if (river.empty()) {
        return 0;
    }
    
    vector<vector<pair<size_t, int>>> vp(river.size());
    vp[0].emplace_back(1, 1);
    int res = INT_MAX;
    for (size_t i = 0; i < vp.size(); ++i) {
        if (river[i] == '0')
        {
            continue;
        }
        for (auto pr : vp[i])
        {
            if (i + pr.first >= vp.size())
            {
                res = min(pr.second, res);
            } else if (river[i + pr.first] == '1')
            {
                vp[i + pr.first].emplace_back(pr.first, pr.second + 1);
            }
            if (i + pr.first + 1 >= vp.size())
            {
                res = min(pr.second, res);
            } else if (river[i + pr.first + 1] == '1')
            {
                vp[i + pr.first + 1].emplace_back(pr.first + 1, pr.second + 1);
            }
        }
    }
    return res;
}
/*
 给A，B 2个array，里面都是integer，已经排好序了，由大到小，他们的长度都是N
 现在从A和B里各选出一个数，总成一个sum，请返回前N个最大的sum
 */
/*
struct Item
{
    int val;
    int index1;
    int index2;
    Item(int v, int i1, int i2) : val(v), index1(i1), index2(i2) {}
};

bool compareItem(Item item1, Item item2)
{
    return item1.val < item2.val;
}

class mycomparison
{
public:
    bool operator() (const Item& lhs, const Item&rhs) const
    {
        return (lhs.val<rhs.val);
    }
};

vector<int> FindMaxSum(vector<int> v1, vector<int> v2, int k)
{
    int n = (int)v1.size();
    int n1 = n-1;
    int n2 = n-1;
    priority_queue<Item, mycomparison> heapx;
    heapx.push(Item(v1[n1]+v2[n2], n1, n2));
    vector<int> res;
    while(res.size() < k)
    {
        Item tmp = heapx.pop();
        int l1 = tmp.index1;
        int l2 = tmp.index2;
        if(l1>0 && l2>0)
        {
            heapx.push(Item(v1[l1-1]+v2[l2], l1-1, l2));
            heapx.push(Item(v1[l1]+v2[l2-1], l1, l2-1));
        }
        res.push_back(tmp.val);
    }
    return res;
}
*/

/*
Longest consecutive sequence
*/

/*
 a.    Flattern this multilevel data structure
 b.    Restore the original structure from the flatterned structure
 
 e.g.
 
 L1 --> L2 --> L3 --> L7 --> L8
                 |
                 v
                 L4 --> L5-->L6
 
 
 WIll be flattened to
 L1 --> L2 --> L3 -->L4 -->L5-->L6-->L7-->L8
 */

struct ListNodeWithChild
{
    int val;
    ListNodeWithChild* next;
    ListNodeWithChild* child;
    ListNodeWithChild(int a): val(a) {}
};

ListNodeWithChild* flatten(ListNodeWithChild* head)
{
    ListNodeWithChild* tmp = head;
    while(tmp)
    {
        if(tmp->child)
        {
            ListNodeWithChild* nextPoint = tmp->next;
            tmp->next = tmp->child;
            while(tmp->next)
            {
                tmp = tmp->next;
            }
            tmp->next= nextPoint;
            tmp->child = nextPoint;
        }
        tmp = tmp->next;
    }
    return head;
}

ListNodeWithChild* unflatten(ListNodeWithChild* head)
{
    ListNodeWithChild* tmp = head;
    ListNodeWithChild* prev = NULL;
    while(tmp)
    {
        if(tmp->next && tmp->child)
        {
            prev = tmp;
            tmp->next = NULL;
            tmp = tmp->child;
            while(!(tmp->next && tmp->child)&& tmp->next)
            {
                tmp = tmp->next;
            }
            prev->next = tmp->next;
            tmp->next = NULL;
            tmp->child = NULL;
            tmp = prev->next;
        }
        else
        {
            tmp = tmp->next;
        }
    }
    return head;
}

bool isMatchDP2(const char *s, const char *p) {
    int l1 = (int)strlen(s);
    int l2 = (int)strlen(p);
    vector<vector<bool>> tracker(l1+1, vector<bool>(l2+1, false));
    tracker[0][0] = true;
    for(int i = 0;i< l1+1;i++)
    {
        for(int j = 1;j<l2+1;j++)
        {
            if(i>0)
            {
                if(tracker[i-1][j-1] && canMatch(s[i-1], p[j-1]))
                {
                    tracker[i][j] = true;
                    continue;
                }
            }
            if(i>0 && j>1)
            {
                if(tracker[i-1][j] && canMatch(s[i-1], p[j-2]) && p[j-1] == '*')
                {
                    tracker[i][j] = true;
                    continue;
                }
            }
            if(j>1)
            {
                if(tracker[i-1][j-2] && p[j-1] == '*')
                {
                    tracker[i][j] = true;
                    continue;
                }
            }
        }
    }
    return tracker[l1][l2];
}

void workerx_1(vector<int> &S, int level, vector<int>& tmp, vector<vector<int>>& res)
{
    res.push_back(tmp);
    
    for(int i = level;i< S.size();i++)
    {
        if(i == level || S[i] != S[i-1])
        {
            tmp.push_back(S[i]);
            workerx(S, i+1, tmp, res);
            tmp.pop_back();
        }
    }
}

vector<vector<int> > subsetsWithDupx(vector<int> &S) {
    vector<vector<int>> res;
    if(S.size() == 0) return res;
    
    sort(S.begin(), S.end());
    vector<int> tmp;
    workerx_1(S, 0, tmp, res);
    return res;
}

/*
 1.anagram, 输出一个句子, 里面的单词是空格隔开, 输出list of anagram in this sentence. 就是List<List<String>>.
 */
vector<vector<string>> groupStr(vector<string> input)
{
    unordered_map<string, vector<string>> map;
    for(auto item : input)
    {
        string s = item;
        sort(s.begin(), s.end());
        map[s].push_back(item);
    }
    
    vector<vector<string>> res;
    for(auto item: map)
    {
        res.push_back(item.second);
    }
    return res;
}

/*
 2.sort colors, 三色旗问题, 用swap, O(n)时间, O(1)空间.
 */
void sortColor(vector<int>& input)
{
    int right = (int)input.size()-1;
    int left = 0;
    int runner = 0;
    while(runner<=right)
    {
        if(input[runner] == 2) swap(input, runner, right);
        else if(input[runner] == 0) swap(input, left++, runner++);
        else runner ++;
    }
}


/*
 print tree in vertical order
 后面这个要先遍历, 边遍历边给每个节点一个index, 比如root为0, 做left减1, right加1.
 然后建立一个HashMap, key是index, value是list<TreeNode>
 */


void findVerticalOrder(TreeNode* root, unordered_map<int, vector<int>>& map, int index, int& minIndex, int& maxIndex)
{
    if(!root) return;
    
    map[index].push_back(root->val);
    minIndex = min(index, minIndex);
    maxIndex = max(index, maxIndex);
    findVerticalOrder(root->left, map, index-1, minIndex, maxIndex);
    findVerticalOrder(root->right, map, index+1, minIndex, maxIndex);
    
}

void printVertical(TreeNode* root)
{
    unordered_map<int, vector<int>> map;
    int minIndex = INT_MAX;
    int maxIndex = INT_MAX;
    findVerticalOrder(root, map, 0, minIndex, maxIndex);
    for(int i = minIndex;i<=maxIndex;i++)
    {
        for(auto item : map[i])
        {
            cout<<item<<" ";
        }
        cout<<endl;
    }
}

/*
 14.jump river的题目, 给一个数组[1,0,1,0,1], 1代表可以站, 0不可以站. 从速度为1开始往前跳, 每次跳的时候, 
 可以跳当前速度那么多格, 也可以跳当前速度+1那么多格. 问最少跳几次可以跳过河(即跳出数组), 或者跳不过河. 
 解法直接递归+cache就可以. 上面的例子跳2次就能跳过河了, 第一次从index=0, 速度为1跳到2, 然后速度为2刚好跳出去
 */
bool canJump(vector<int> input)
{
    int n = (int) input.size();
    vector<unordered_set<int>> tracker(n);
    tracker[0].insert(1);
    for(int i = 0;i< n;i++)
    {
        for(auto item: tracker[i])
        {
            int index = i + item; // current location + speed;
            if(index >=n) return true;
            if(input[index] != 0)
            {
                if(tracker[index].find(item) == tracker[index].end()) tracker[index].insert(item);
            }
            if(input[index+1]!=0)
            {
                if(tracker[index+1].find(item+1) == tracker[index+1].end()) tracker[index+1].insert(item+1);
            }
        }
    }
    return false;
}

/*
 There are a row of houses, each house can be painted with three colors red, blue and green. 
 The cost of painting each house with a certain color is different. You have to paint all 
 the houses such that no two adjacent houses have the same color. You have to paint the houses 
 with minimum cost. How would you do it?
 */

int paintHouse(vector<vector<int>> input)
{
    int n = (int)input.size();
    int m = (int)input[0].size();
    
    vector<vector<int>> tracker(n, vector<int>(m, INT_MAX));
    for(int j = 0;j<m;j++)
    {
        tracker[0][j] = input[0][j];
    }
    for(int i = 1;i<n;i++)
    {
        for(int j = 0;j<m;j++)
        {
            int minVal = INT_MAX;
            for(int k = 0;k<m;k++)
            {
                if(j != k)
                {
                    minVal = min(input[i][j] + tracker[i-1][k], minVal);
                }
            }
            tracker[i][j] = minVal;
        }
    }
    int res = INT_MAX;
    for(int j = 0;j<m;j++)
    {
        res = min(res, tracker[n-1][m]);
    }
    return res;
}

/*
 2. a) 给你棵二叉树，节点上有权值，问从一个叶子走到另外一个叶子的路里面权值最大的那条是什么
 */
int findMaxA(TreeNode* root, TreeNode* n1, TreeNode* n2, int& count)
{
    if(!root) return INT_MIN;
    if(n1 == n2) return n1->val;
    int left = findMaxA(root->left, n1, n2, count);
    int right = findMaxA(root->right, n1, n2, count);
    int res;
    
    if((left > INT_MIN && right > INT_MIN) ||
       (root == n1 || root == n2))
    {
        res = max(left, right);
        res = max(res, root->val);
    }
    else
    {
        res = max(left, right);
        if(count == 1)
            res = max(res, root->val);
    }
    if(root == n1 || root== n2) count++;
    
    return res;
}

void findMaxB(TreeNode* root, TreeNode* n1, TreeNode* n2, int& maxVal)
{
    if(!root)
    {
        maxVal = INT_MIN;
        return;
    }
    int leftMax;
    int rightMax;
    findMaxB(root->left, n1, n2, leftMax);
    findMaxB(root->right, n1, n2, rightMax);
    if(leftMax > INT_MIN && rightMax> INT_MIN)
    {
        maxVal = max(max(leftMax, rightMax), root->val);
    }
    else if( leftMax != INT_MIN)
    {
        maxVal = max(leftMax, root->val);
    }else if(rightMax != INT_MIN)
    {
        maxVal = max(rightMax, root->val);
    }
    else if(root != n1 && root != n2)
    {
        maxVal = INT_MIN;
    }
    else
    {
        maxVal = root->val;
    }
    return;
    
}

int findMaxC(TreeNode* root, TreeNode* n1, TreeNode* n2, int & maxVal)
{
    if(!root)
    {
        maxVal = INT_MIN;
        return 0;
    }
    int leftMax;
    int rightMax;
    int leftCount;
    int rightCount;
    leftCount = findMaxC(root->left, n1, n2, leftMax);
    rightCount = findMaxC(root->right, n1, n2, rightMax);
    
    if(leftCount == 2 || rightCount == 2)
    {
        maxVal = max(leftMax, rightMax);
        return 2;
    }
    else if(root == n1 || root == n2)
    {
        maxVal = max(max(leftMax, rightMax), root->val);
        return 1 + leftCount + rightCount;
    }
    else if(leftCount == 0 && rightCount == 0)
    {
        maxVal = INT_MIN;
        return 0;
    }
    else
    {
        maxVal = max(max(leftMax, rightMax), root->val);
        return leftCount + rightCount;
    }
}


TreeNode *LCA(TreeNode *root, TreeNode *p, TreeNode *q) {
    if (!root) return NULL;
    if (root == p || root == q) return root;
    TreeNode *L = LCA(root->left, p, q);
    TreeNode *R = LCA(root->right, p, q);
    if (L && R) return root;  // if p and q are on both sides
    return L ? L : R;  // either one of p,q is on one side OR p,q is not in L&R subtrees
}

/*
 2.给出一个无重复数字的array,输出其subset。 例如[1, 2]，输出[] [1] [2] [1, 2]
 我一开始用了dfs，面试官check过无误后让我写成iterative. 稍微想了一下，就用位运算做了，
 每个独立subset都可以用长度为array大小的二进制串表示，0代表不存在，1代表存在，只要从0 枚举到2^n-1就 ok了。
 */
void PushVector(string s, vector<vector<int>>& res)
{
    vector<int> tmp;
    for(int i = 0;i< s.length();i++)
    {
        if(s[i] == '1')
        {
            tmp.push_back(i+1);
        }
    }
    res.push_back(tmp);
}

vector<vector<int>> allSubset(int n)
{
    string s1;
    for(int i = 0;i<n;i++)
    {
        s1 += '0';
    }
    string s2 = "1";
    int base = 2;
    vector<vector<int>> res;
    for(int i = 0;i < pow((double)base,n);i++)
    {
        PushVector(s1, res);
        s1 = addBinary(s1, s2);
    }
    return res;
}

/*
 给二叉树增加next节点。我开始用了递归，时间复杂度o(n)空间O(h)。
 */
void AddNext(TreeLinkNode* root)
{
    if(!root ||(!root->left && !root->right)) return;
    TreeLinkNode* tmp;
    tmp = root->next;
    while(tmp && !tmp->left && !tmp->right)
    {
        tmp = tmp->next;
    }

    TreeLinkNode* next;
    if(root->left && root->right)
    {
        root->left->next = root->right;
    }
    next = root->right? root->right : root->left;
    
    if(tmp)
    {
        next->next = tmp->left? tmp->left: tmp->right;
    }
    AddNext(root->right);
    AddNext(root->left);
    
}

/*
4. (1) 一维度向量相乘。每个向量很长，billion个数字。
*/
int dotMulti(vector<pair<int, int>> v1, vector<pair<int,int>> v2)
{
    int n = (int)v1.size();
    int m = (int)v2.size();
    int res = 0;
    int first = 0;
    int second = 0;
    while(first<n && second < m)
    {
        if(v1[first].first == v2[second].first)
        {
            res += v1[first].second * v2[second].second;
            first++;
            second++;
        }
        else if(v1[first].first < v2[second].first)
        {
            first++;
        }
        else
        {
            second++;
        }
    }
    return res;
}

/*
 店面，第一轮。上来先问了很多简历的东西。然后就做题，minimum window substring
 , 但字典没有重复。直接给出了最优解，太紧张有两个bug，很快改正。跑了一个test
 case, 通过了。然后就问了一下想去什么样的组？然后又随便聊了一下就挂了。
*/
int minSubStrNoDup(string s1, string p)
{
    bool tracker[256];
    fill_n(&tracker[0], 256, false);
    for(int i = 0;i< p.length();i++)
    {
        tracker[p[i]] = true;
    }
    bool stracker[256];
    fill_n(&stracker[256], 256, false);
    int left = 0;
    while(!tracker[s1[left]]) left++;
    int count = 0;
    int minLength = INT_MAX;
    for(int i = left;i< s1.length();i++)
    {
        if(tracker[p[i]])
        {
            if(!stracker[p[i]])
            {
                stracker[p[i]] = true;
                count++;
            }
            else
            {
                while(left <i)
                {
                    if(s1[left] != s1[i]) break;
                    
                    if(tracker[s1[left]])
                    {
                        count--;
                        stracker[s1[left]] = false;
                    }
                }
            }
        }
        if(count == p.length())
        {
            if(i-left+1 < minLength) minLength = i-left+1;
        }
    }
    return minLength;
}

/*
 1) 给个数组seq， 和一个total，找 if there is a contiguous sequence in seq which sums to total.
 都是正数， 第一次没注意contiguous，给了个back tracking的解法。然后说是contiguous， 给了
 个维护窗口的解法，不过犯了个小错误。时间过去了半小时。。。
 */
bool findTarget(vector<int> input, int target)
{
    int n = (int)input.size();
    unordered_map<int, vector<int>> tracker;
    int sum = 0;
    for(int i = 0;i< n;i++)
    {
        sum += input[i];
        tracker[sum].push_back(i);
    }
    
    for(auto item: tracker)
    {
        if(tracker.find(item.first+target) != tracker.end())
        {
            int minLeft = item.second[0];
            int maxRight = tracker[item.first + target][tracker[item.first + target].size()-1];
            if(maxRight>=minLeft)
                return true;
        }
    }
    return false;
}

int findMinXI(vector<int> input)
{
    int n = (int)input.size();
    int l = 0;
    int r = n-1;
    int minVal = INT_MAX;
    while(l<n-1)
    {
        int mid = l + (r-l)/2;
        if(input[l] < input[mid])
        {
            minVal = min(input[l], minVal);
            l = mid+1;
        }
        else if(input[r] < input[mid])
        {
            minVal = min(input[r], minVal);
            r = mid -1;
        }
        else
        {
            l++;
        }
    }
    minVal = min(input[l], minVal);
    minVal = min(input[r], minVal);
    return minVal;
}


int main(int argc, const char * argv[])
{
    /*
    int A[] = {1};
    int B[] = {2,3,4,5,6,7};
    
    double res = findMedianSortedArrays(A, 1, B, 6);
    cout<<res<<endl;
*/
/*    string a = "aa";
    string b = "aa";
 
    bool res = isMatch(a.c_str(),b.c_str());
    cout << res<<endl;
*/
    //vector<int> a = {1,2};
    //nextPermutation(a);
    
/*    int A[] = {1,2,0};
    int res = firstMissingPositive(A, 3);
    cout<<res<<endl;
*/
    /*
    bool res = isMatchII("", "*a*");
    cout<<res<<endl;
     */
/*    vector<int> A = {2,4};
    int res = largestRectangleArea(A);
    cout<<res<<endl;
 */
/*    string s1 = "aabcc";
    string s2 = "dbbca";
    string s3 = "aadbbcbcac";
    bool res = isInterleave(s1, s2, s3);
    cout<<res<<endl;
    
    s3 = "aadbbbaccc";
    res = isInterleave(s1, s2, s3);
    cout<<res<<endl;
 */
    //int res = divide(5, 2);
    //cout<<res<<endl;
    //int A[] = {1,2,3};
    //int res = jump(A, 3);
    //cout<<res<<endl;
    
    //vector<int> num = {1,1,1,2,2};
    //vector<vector<int>> res = permuteUnique(num);
    
    //ListNode* input = new ListNode(1);
    //TreeNode* res = sortedListToBST(input);
    //cout<<res->val<<endl;
    
    //permutation(3);
    /*
    vector<char> r1 = {'O','X','O'};
    vector<char> r2 = {'X','O','X'};
    vector<char> r3 = {'O','X','O'};
    //vector<char> r4 = {'X','O','X','X'};
    
    vector<vector<char>> board = {r1,r2,r3};
    solve(board);
    
    vector<char> rr1 = {'O','O','O'};
    vector<char> rr2 = {'O','O','O'};
    vector<char> rr3 = {'O','O','O'};
    
    board = {rr1,rr2,rr3};
    solve(board);
     */
    
    /*
    int a[] = { 1,3,5,7,9,2,4,6,8,10};
    interleave(a,10);
    for(int i= 0;i< 10;i++)
    {
        cout<<a[i]<<endl;
    }
     */
    /*
    for(int i = 0;i< INT_MAX;i++)
    {
        //SET_BIT = 0;
        //totalBits(i);
        //cout<< i << " : " << SET_BIT<<endl;
        //cout<< i << " : " <<countSetBits(i)<<endl;
        countSetBits(i);

    }
    cout<<"done"<<endl;
     */
    /*int A[] = {1,2,2};
    int res = removeDuplicatesII(A, 3);
    cout<<res;
     */
    /*
    TreeNode* root = new TreeNode(1);
    TreeNode* left = new TreeNode(2);
    root->left = left;
    vector<int> res = postorderTraversal(root);
    cout<<res[0]<<endl;
     */
    //int A[] = {1,3};
    //int res = searchInRotateArray(A,2,3);
    //cout<<res<<endl;
    //vector<int> input = {0,4,3,0};
    //vector<int> res = twoSum(input, 0);
    //cout<<res[0]<< "   "<<res[1]<<endl;
    //vector<int> input = {0,0,0,0};
    //vector<vector<int>> res = fourSum(input, 0);
    //int a[] = {5,5,1,7,1,1,5,2,7,6};
    //int res = trapII(a, 10);
    //cout<<res;
    //res = trap(a, 10);
    //cout<<res;
    /*
    ListNode* input = new ListNode(5);
    ListNode* res = reverseBetween(input, 1, 1);
    cout<<res->val<<endl;
    */
    //ListNode* input = new ListNode(1);
    //ListNode* res = removeNthFromEnd(input, 1);
    //cout<<res->val  <<endl;
    //RandomListNode* input = new RandomListNode(1);
    //RandomListNode* res = copyRandomListII(input);
    //cout<<res->label<<endl;
    //2,[set(2,1),set(1,1),set(2,3),set(4,1),get(1),get(2)]
    
    /*LRUCache* cache = new LRUCache(2);
    int res;
    
    cache->set(2,1);
    cache->set(1,1);
    cache->set(2,3);
    cache->set(4,1);
    res = cache->get(1);
    cout<<res<<endl;
    res = cache->get(2);
    cout<<res<<endl;
     */
    //bool res = isPalindrome(".,");
    //cout<<res<<endl;
    
    //char* res = strStrII("aaa", "aa");
    //cout<<res<<endl;
    
    //int res = atoi("    -010");
    //cout<<res;
    
    //bool res = isMatchX("a", ".*..a*");
    //cout << res <<endl;
    
    //vector<string> input = {"aa", "aa"};
    //string res = longestCommonPrefix(input);
    //cout<< res <<endl;
    
    //simplifyPath("/../");
    
    //string res = countAndSay(4);
    //cout<<res<<endl;
    //int res = longestValidParenthesesII("(()");
    //cout<<res<<endl;
    
    //vector<int> input = {2,1,2};
    //int res = largestRectangleAreaII(input);
    //cout<<res<<endl;
    
    //TreeNode* root = new TreeNode(1);
    //TreeNode* right = new TreeNode(2);
    ///root->right = right;
    //flatten(root);

    /*
    TreeLinkNode* root = new TreeLinkNode(1);
    TreeLinkNode* r1 = new TreeLinkNode(2);
    TreeLinkNode* r2 = new TreeLinkNode(3);
    TreeLinkNode* r3 = new TreeLinkNode(4);
    TreeLinkNode* r4 = new TreeLinkNode(5);
        TreeLinkNode* r5 = new TreeLinkNode(6);
        TreeLinkNode* r6 = new TreeLinkNode(7);
        TreeLinkNode* r7 = new TreeLinkNode(8);
    root->left = r1;
    root->right = r2;
    r1->left = r3;
    r1->right = r4;
    r3->left = r6;
    r2->right = r5;
    r5->right = r7;
    
    connectII(root);
     */
    /*
    vector<int> inorder = {2,1};
    vector<int> pre = {1,2};
    TreeNode* res = buildTreeII(pre, inorder);
    cout<<res->val;
    */
    /*
    ListNode* root = new ListNode(0);
    ListNode* r1 = new ListNode(2);
    ListNode* r2 = new ListNode(5);
    
    root->next = r1;
    r1->next = r2;
    vector<ListNode*> input = {root};
    ListNode* res = mergeKLists(input);
    */
    /*
    ListNode* root = new ListNode(1);
    ListNode* r1 = new ListNode(-1);
    root->next = r1;
    
    //ListNode* res = insertionSortList(root);
    ListNode* res = sortList(root);
    */
    /*
    int A[] = {7,2,5,4,1,6};
    int res = firstMissingPositiveII(A, 6);
    cout<< res<<endl;
     */
    /*int A[] = {1,0};
    sortColorsII(A, 2);
    int B[] ={2,0};
    sortColorsII(B, 2);
     */
    /*
    vector<vector<int>> input(1,vector<int>(1,1));
    bool res = searchMatrix(input, 0);
    cout<<res<<endl;
    */
    
    /**
     vector<int> input = {1,2,2,4,4};
    vector<vector<int>> res = permuteUniqueII(input);
    for(auto i: res)
    {
        for(auto j : i)
        {
            cout<<j<< " ";
        }
        cout<<endl;
    }
    cout<<res.size()<<endl;
     */
    
    /*
     vector<string> res = letterCombinations("2");
    cout<<res.size()<<endl;
     */
    /*
    vector<string> res = getString("a");
    cout<<res.size()<<endl;
    for(auto i : res)
    {
        //if(i == "miss")
            cout<<i<<endl;
    }
    */
    
   /*
    unordered_set<string> dict;
    dict.insert("most");
    dict.insert("mist");
    dict.insert("miss");
    dict.insert("lost");
    dict.insert("fist");
    dict.insert("fish");
    int res = ladderLength("lost", "miss", dict);
    cout<< res<<endl;
     */
    /*
    vector<vector<string>> res = partition("aab");
    cout<<res.size()<<endl;
     */
    /*
    int res = minCut("aa");
    cout<<res<<endl;
     */
    /*
    vector<vector<int>> input (2, vector<int>(2,0));
    int res = uniquePathsWithObstacles(input);
    cout<<res;
     */
    /**
    vector<string> res = restoreIpAddresses("0000");
    cout<<res.size()<<endl;
     */
    
    /*vector<int> input = {1,1,1,3,3, 5};
    vector<vector<int>> res = combinationSum2(input, 8);
    cout<<res.size()<<endl;
    */
    //int res = divideII(-1010369383, -2147483648);
    //cout<<res<<endl;
    //int res = sqrt(1);
    //cout<<res<<endl;
    //int A[] = {-2,1,-3,4,-1,2,1,-5,4};
    //int res = maxSubArray(A, 9);
    //cout<<res<<endl;
    
    //vector<int> input = {2,4};
    //int res = maxHistogram(input);
    //cout<<res<<endl;
    
    //int res = lengthOfLongestSubstring("wlrbbmqbhcdarzowkkyhiddqscdxrjmowfrxsjybldbefsarcbynecdyggxxpklorellnmpapqfwkhopkmco");
    //cout<<res<<endl;
    //vector<int> input = {1,2,4};
    //int res = maxProfitIII(input);
    //cout<<res<<endl;
    /*
    bool res;
    string s1 = "aabcc";
    string s2 = "dbbca";
    string s3 = "aadbbcbcac";
    string s4 = "aadbbbaccc";
    res = isInterleaveII(s1, s2, s3);
    cout<<res<<endl;
    res = isInterleaveII(s1, s2, s4);
    cout<<res<<endl;
    res = isInterleaveDpII(s1, s2, s3);
    cout<<res<<endl;
    res = isInterleaveDpII(s1, s2, s4);
    cout<<res<<endl;
    */
    /*vector<int> firstRow = {1,2};
    vector<int> secondRow = {1,1};
    vector<vector<int>> input = {firstRow, secondRow};
    int res = minPathSum(input);
    cout<<res;
     */
    //int res = minDistanceII("a", "b");
    //cout<<res<<endl;
    //int res = numDecodings("10");
    //cout<<res<<endl;
    //int res = numDistinctII("ccc", "cc");
    //cout<<res<<endl;
    //res = numDistinctOld("ccc", "cc");
    //cout<<res<<endl;
    /*
    vector<string> res;
    unordered_set<string> dict = {"cat", "cats", "and", "sand", "dog"};
    res = wordBreakII("catsanddog",dict);
    cout<< res.size()<<endl;
    for(int i = 0;i< res.size();i++)
    {
        cout<<res[i]<<endl;
    }
    */
    //bool res = isPalindrome(1);
    //cout<<res<<endl;
    
    /*string s = "aaa";
    vector<string> L = {"a", "a"};
    vector<int> res = findSubstring(s, L);
    for(int i = 0;i<res.size();i++)
    {
        cout<< res[i]<<endl;
    }
    */
    
    //vector<vector<int>> res = generate(3);
    //cout<<res.size()<<endl;
    /*
    Interval newInterval(6,8);
    Interval input1(1,5);
    vector<Interval> input {input1};
    vector<Interval> res = insertDCII(input, newInterval);
    cout<< res.size()<<endl;
    */
    
    /*
     vector<int> a1 = {1,2};
    vector<int> a2 = {3,4};
    vector<vector<int>> input = { a1,a2};
    vector<int> res = spiralOrderII(input);
    for(auto a : res)
    {
        cout<<a<<endl;
    }
     */
    //vector<string> input = {"This", "is", "an", "example", "of", "text", "justification."};

    /*vector<string> input = {""};
    vector<string> res = fullJustifyII(input,2);
    for(auto s : res)
    {
        cout << s<<endl;
    }
    */
    //string s = " ";
    //reverseWords(s);
    //cout<<s<<endl;
    //vector<int> input = {2,3,1};
    //int res = findMin(input);
    //cout<<res<<endl;
    
    //vector<vector<int>> res = combine(3,2);
    /*
    TreeNode* r1 = new TreeNode(1);
    TreeNode* r2 = new TreeNode(2);
    TreeNode* r3 = new TreeNode(3);
    TreeNode* r4 = new TreeNode(4);
    r1->left = r2;
    r1->right = r3;
    r2->right = r4;
    
    res.clear();
    serializeTree(r1);
    i = 0;
    TreeNode* res = deserialize();
    cout<<res->val<<endl;
    */
    //vector<vector<int>> res = printAllFactor(8);
    //vector<int> input = {1,2,3,4,5,6,7,8,5,4,2};
    //int res = findMax(input);
    //cout<<res<<endl;
    
    /*
    vector<int> res = findSumOfPerfectSqaure(41);
    for(auto item:res)
    {
        cout<<item<<endl;
    }
    */

    /*
    vector<string> res = tokenize("abc     bbb ccc");
    for(auto item : res)
    {
        cout<<item<<endl;
    }
    */
    
    /*
    vector<Elevation> p1 = {Elevation(3), Elevation(2), Elevation(1)};
    vector<Elevation> p2 = {Elevation(4), Elevation(5), Elevation(1)};
    vector<Elevation> p3 = {Elevation(5), Elevation(6), Elevation(1)};
    vector<vector<Elevation>> input = {p1,p2,p3};
    
    vector<Elevation> res = canFlowinto(input);
    for(auto item : res)
    {
        cout<<item.val<<endl;
    }
    */
    
//    int input[] = {1,8,10,4,5,6,7,3,9,2};
    //int input[] = {10,4,8,2,1,6,5,3,7,9};
    //quickSort(input, 10);
    /*
    vector<int> input = {1,3,50,75};
    vector<string> res = complementry(input, 99);
    for(auto item : res)
    {
        cout<<item<<", ";
    }
    */
    /*
    string input = "ab cd ef  g   hij  ";
    string res = removeSpace(input);
    cout<<res<<endl;
     */
    /*
    vector<int> res = printMultiples(20);
    for(auto item : res)
    {
        cout<<item<<endl;
    }
     */
    
    // insert code here...
    //std::cout << "Hello, World!\n";
/*
    vector<int> input = {1,8,10,4,5,6,7,3,9,2,-1,-4,-3};
    findFirstk(input, 0, (int)input.size()-1, 7);

    for(int i = 0;i<input.size();i++)
    {
        cout<<input[i] <<endl;
    }
*/
  
    /*
    int res = numDecodings("10");
    cout<<numDecodings("10")<<endl;
    cout<<numDecodings("90")<<endl;
    cout<<numDecodingsDP("10")<<endl;
    cout<<numDecodingsDP("90")<<endl;
    */
    
    //cout<<sqrt(2.0, 0.001) <<endl;
    //cout<<sqrt(0.5, 0.001) <<endl;
    //string s = "TADTTTTBDB";
    //s = sort(s);
    //cout<<s<<endl;
    /*
    workItem w1(1,10,2);
    workItem w2(9,20,3);
    workItem w3(10,25,3);
    workItem w4(26,30,2);
    
    vector<workItem> input = {w1, w2, w3,w4};
    int res = findMax(input);
    cout<<res<<endl;
    */

    /*vector<float> input = {1.2, 2.5, 9.3};
    int res = findClosest(input, 5);
    cout<<res<<endl;
    res = findClosest(input, 1);
    cout<<res<<endl;
    res = findClosest(input, 10);
    cout<<res<<endl;
    input = {1.2};
    res = findClosest(input, 1);
    cout<<res<<endl;
    res = findClosest(input, 10);
    cout<<res<<endl;
    */
    /*
    vector<int> input = {1,8,10,4,5,6,7,3,9,2};
    findFirstk(input, 0, (int)input.size()-1, 8);
    for(auto item : input)
    {
        cout<<item<<endl;
    }
     */
    /*
    vector<int> input = {1,2,3,4};
    vector<vector<int>> res = permute(input);
    int sum = 0;
    for(auto item: res)
    {
        int tmp = 0;
        for(int j = 0;j<item.size();j++)
        {
            tmp = tmp* 10 + item[j];
        }
        sum += tmp;
    }
    
    cout<<sum<<endl;
    */
    /*
    insert(1);
    insert(6);
    insert(4);
    insert(2);
    insert(5);
    insert(3);

    cout<<get()<<endl;
    cout<<get()<<endl;
    cout<<get()<<endl;
    cout<<get()<<endl;
    cout<<get()<<endl;
    cout<<get()<<endl;
    cout<<get()<<endl;
    */
    /*
    int A[] = {1,8,10,4,5,6,7,3,9,2,-1,-4,-3};
    quickSort(A, 13);
    for(int i = 0;i<13; i++)
    {
        cout<<A[i] <<endl;
    }
     */
    
    /*
    TreeNode* n1 = new TreeNode(1);
    TreeNode* n2 = new TreeNode(2);
    TreeNode* n3 = new TreeNode(3);
    TreeNode* n4 = new TreeNode(4);
    TreeNode* n8 = new TreeNode(8);
    TreeNode* n5 = new TreeNode(5);
    TreeNode* n9 = new TreeNode(9);
    TreeNode* n6 = new TreeNode(6);
    TreeNode* n7 = new TreeNode(7);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n3->left = n6;
    n3->right = n7;
    n4->left = n8;
    n4->right = n9;
    
    
    int res = findMax(n1, n2, n9);
    cout<<res<<endl;
    */
    
    /*
    vector<int> a = {1,2,3,4,5};
    vector<int> res = multiresult(a);
    for(int i = 0;i< 5;i++)
    {
        cout<< res[i]<<endl;
    }
    */
    
    /*
    bool res = isOneEditDistance("abcdefg", "abcdefgh");
    cout<<res<<endl;
    res = isOneEditDistance("abcdxfg", "abcdefg");
    cout<<res<<endl;
    res = isOneEditDistance("abcdefg", "abchdefg");
    cout<<res<<endl;
    res = isOneEditDistance("abcdefg", "axcoefg");
    cout<<res<<endl;
    res = isOneEditDistance("abcdefg", "abcdefg");
    cout<<res<<endl;
    
    
    res = isOneEditDistanceII("abcdefg", "abcdefgh");
    cout<<res<<endl;
    res = isOneEditDistanceII("abcdxfg", "abcdefg");
    cout<<res<<endl;
    res = isOneEditDistanceII("abcdefg", "abchdefg");
    cout<<res<<endl;
    res = isOneEditDistanceII("abcdefg", "axcoefg");
    cout<<res<<endl;
    res = isOneEditDistanceII("abcdefg", "abcdefg");
    cout<<res<<endl;
    */
    
    /*
    Interval p1(1, 4);
    Interval p2(2, 5);
    Interval p3(3, 6);
    Interval p4(4, 7);
    vector<Interval> input = {p1,p2,p3,p4};
    
    int res = maxOverlap(input);
    cout<<res<<endl;
    */
    
    /*
     “( * 1 ( + 1 2 3 ) )” => 6
     “( * ( + 1 1 ) 17 )” => 34
     “7” => 7
     ( * ( + 1 1 ) 17 )
     ( * 17 ( + 1 1 ) )
     */
    /*
    int res;
    vector<string> input = {"(",  "*",  "1",  "(",  "+", "1",  "2", "3", ")", ")"};
    res = calculat(input);
    cout<<res<<endl;
    
    input = {"(", "*", "(", "+", "1", "1", ")", "17", ")"};
    res = calculat(input);
    cout<<res<<endl;
    
    input = {"7"};
    res = calculat(input);
    cout<<res<<endl;
    */
    /*
    string s1 = "OldSite:LeetCode.org";
    string s2 = "NewSite:LeetCodeOj.com";
    int res = longestCommonSubString(s1, s2);
    cout<<res<<endl;
    */
    
    /*
    vector<int> n1 = {1,2,3,4};
    vector<int> n2 = {};
    vector<int> n3 = {5};
    vector<int> n4 = {6,7,8};
    
    backend = {n1, n2, n3, n4};
    currentParentIterator = backend.begin();
    currentChildIterator = (*currentParentIterator).begin();
    while(hasNext())
    {
        cout<<getNext()<<endl;
    }
    */
    
    
    /*
    TreeNode* n1 = new TreeNode(1);
    TreeNode* n2 = new TreeNode(2);
    TreeNode* n3 = new TreeNode(3);
    TreeNode* n4 = new TreeNode(4);
    TreeNode* n8 = new TreeNode(8);
    TreeNode* n5 = new TreeNode(5);
    TreeNode* n9 = new TreeNode(9);
    TreeNode* n6 = new TreeNode(6);
    TreeNode* n7 = new TreeNode(7);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n4->left = n6;
    n4->right = n7;
    n6->left = n8;
    n6->right = n9;
    
    TreeNode* res = InverseATree(n1);
    cout<< res->val<<endl;
    */
    
    /*
    Interval n1(16,21);
    Interval n2(8,9);
    Interval n3(25,30);
    Interval n4(5,8);
    Interval n5(15,23);
    Interval n6(0,3);
    Interval n7(6,10);
    Interval n8(17,19);
    Interval n9(26,26);
    Interval n10(19,20);
    
    IntervalTreeNode* it1 = new IntervalTreeNode(n1, 30);
    IntervalTreeNode* it2 = new IntervalTreeNode(n2, 23);
    IntervalTreeNode* it3 = new IntervalTreeNode(n3, 30);
    IntervalTreeNode* it4 = new IntervalTreeNode(n4, 10);
    IntervalTreeNode* it5 = new IntervalTreeNode(n5, 23);
    IntervalTreeNode* it6 = new IntervalTreeNode(n6, 3);
    IntervalTreeNode* it7 = new IntervalTreeNode(n7, 10);
    IntervalTreeNode* it8 = new IntervalTreeNode(n8, 20);
    IntervalTreeNode* it9 = new IntervalTreeNode(n9, 26);
    IntervalTreeNode* it10 = new IntervalTreeNode(n10,20);
    it1->left = it2;
    it1->right = it3;
    
    it2->left = it4;
    it2->right = it5;
    
    it4->left = it6;
    it4->right = it7;
    
    it3->left = it8;
    it3->right = it9;
    
    it8->right = it10;
    
    

    //bool res;
    cout<<isCovered(it1, 23)<<endl;
    cout<<isCovered(it1, 24)<<endl;
    */
    
    /*
    char input[] = {'c', 'f', 'j', 'p', 'v'};
    cout<<findStrictly(input, 5, 'a')<<endl;
    cout<<findStrictly(input, 5, 'g')<<endl;
    */
    
    /*
    int input[] = {1,2};
    int res;
    res = findMax(input, 2);
    cout<<res<<endl;
    
    int input2[] = {2,1};
    res = findMax(input2, 2);
    cout<<res<<endl;
    
    int input3[] = { 1, 2, 3};
    res = findMax(input3, 3);
    cout<<res<<endl;
    
    int input4[] = { 1, 2, 3, 7, 6, 0};
    res = findMax(input4, 6);
    cout<<res<<endl;
    
    int input5[] = { 1, 2, 3, 6, 7, 9};
    res = findMax(input5, 6);
    cout<<res<<endl;
    */
    
    /*
    int res = 0;
    res = countOfPalindrome("aba");
    cout<<res<<endl;
    
    res = countOfPalindrome("abba");
    cout<<res<<endl;
    */
    
    /*
    TreeNode* n1 = new TreeNode(1);
    TreeNode* n2 = new TreeNode(2);
    TreeNode* n3 = new TreeNode(3);
    TreeNode* n4 = new TreeNode(4);
    TreeNode* n8 = new TreeNode(8);
    TreeNode* n5 = new TreeNode(5);
    TreeNode* n9 = new TreeNode(9);
    TreeNode* n6 = new TreeNode(6);
    TreeNode* n7 = new TreeNode(7);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n4->left = n6;
    n4->right = n7;
    n6->left = n8;
    n6->right = n9;
    
    vector<vector<int>> res = printout(n1);
    for(auto item : res)
    {
        for(auto num: item)
        {
            cout<< num << " ";
        }
        cout<<endl;
    }
     */
    
    /*
    cout<<isValidPalindrome("aba")<<endl;
    cout<<isValidPalindrome("abaccc")<<endl;
    */
    
    /*
    Interval v1(1,2);
    Interval v2(1,3);
    Interval v3(2,4);
    Interval v4(3,5);
    Interval v5(3,6);
    Interval v6(4,7);
    Interval v7(5,8);
    Interval v8(5,9);
    vector<Interval> input = {v1,v2,v3,v4,v5,v6,v7,v8};
    printTree(input);
    
    */
    
    /*
    cout<< pow(9,99)<<endl;
    */
    
    /*
    vector<int> input = {1,2,3,4,5};
    setContainer(input, 3);
    while(hasNextPermutation())
    {
        printNext();
    }
    */
    /*
    vector<int> input = {1,5,0,6};
    cout<<canSeparate(input, 21)<<endl;
    cout<<canSeparate(input, 25)<<endl;
    */
    /*
    vector<int> input = {1,2,3,4};
    vector<int> res = shuffleSort(input);
    for(int i = 0;i< (int)res.size();i++)
    {
        cout<<res[i]<<" ";
    }
    cout<<endl;
    */
    /*
    int a[] = {7,2,7,3,2,5};
    shuffleOn(a, 7);
    for(int i = 0;i<7;i++)
    {
        cout<<a[i]<< " ";
    }
    cout<<endl;
     */
    
    /*
    int a[] = {7,2,8,7,3,2,5};
    rightRotate(a, 7, 3);
    for(int i = 0;i<7;i++)
    {
        cout<<a[i]<< " ";
    }
    cout<<endl;
    */
    
    /*
    vector<int> input = {2,2,2,2,2,2,2,2,2,8,8,8};
    int res = countPair(input, 10);
    cout<<res<<endl;
    */
    
    /*
    vector<int> input = { 1,5,20,25};
    int res = findChange(input, 40);
    cout<<res<<endl;
    */
    
    /*
    int x,y,z;
    vector<vector<vector<int>>> res;
    
    
    res = FindAllRussiaBlocks(0);
    cout<<res.size()<<endl;
    
    res = FindAllRussiaBlocks(1);
    cout<<res.size()<<endl;
    res = FindAllRussiaBlocks(2);
    cout<<res.size()<<endl;
    
    res = FindAllRussiaBlocks(3);
    cout<<res.size()<<endl;
    z = (int)res.size();
    y = (int)res[0].size();
    x = (int)res[0][0].size();
    for(int k = 0;k<z;k++)
    {
        for(int j = 0;j<y;j++)
        {
            for(int i = 0;i<x;i++)
            {
                cout<< res[k][j][i] << " ";
            }
            cout<<endl;
        }
        cout<<"-----------------------"<<endl;
    }
    cout<<endl;
    
    
    
    res = FindAllRussiaBlocks(4);
    cout<<res.size()<<endl;
    z = (int)res.size();
    y = (int)res[0].size();
    x = (int)res[0][0].size();
    for(int k = 0;k<z;k++)
    {
        for(int j = 0;j<y;j++)
        {
            for(int i = 0;i<x;i++)
            {
                cout<< res[k][j][i] << " ";
            }
            cout<<endl;
        }
        cout<<"-----------------------"<<endl;
    }
    cout<<endl;
    */
    
    /*
    string res;
    res = convertFrom10To18(29);
    cout<<res<<endl;
    */
    
    
    /*
    int res;
    vector<int> input = {1,2,3,2,3,2,2};
    res = findMajority(input);
    cout<<res<<endl;
    
    vector<int> n2 = {1,2,3,2,3,2};
    res = findMajority(n2);
    cout<<res<<endl;;
    */
    
    /*
     TreeNode* n1 = new TreeNode(9);
     TreeNode* n2 = new TreeNode(5);
     TreeNode* n3 = new TreeNode(3);
     TreeNode* n4 = new TreeNode(3);
     TreeNode* n8 = new TreeNode(1);
     TreeNode* n5 = new TreeNode(1);
     TreeNode* n9 = new TreeNode(1);
     TreeNode* n6 = new TreeNode(1);
     TreeNode* n7 = new TreeNode(1);
     
     n1->left = n2;
     n1->right = n3;
     n2->left = n4;
     n2->right = n5;
     n3->left = n6;
     n3->right = n7;
     n4->left = n8;
     n4->right = n9;
    
    
    TreeNode* res = findKthTreeNode(n1, 10);
    cout<<res->val<<endl;
     */
    
    /*
    string res;
    //res = get_decimal(1,19);
    //cout<<res<<endl;

    res = divideII(20,200);
    cout<<res<<endl;
    */
    /*
    int res = 0;
    res = findNthDigit(10);
    cout<<res<<endl;
    
    res = findNthDigit(11);
    cout<<res<<endl;
    
    res = findNthDigit(12);
    cout<<res<<endl;
    
    res = findNthDigit(13);
    cout<<res<<endl;
    
    res = findNthDigit(100);
    cout<<res<<endl;
    
    res = findNthDigit(101);
    cout<<res<<endl;
    
    res = findNthDigit(102);
    cout<<res<<endl;
    
    res = findNthDigit(105);
    cout<<res<<endl;
    
    res = findNthDigit(106);
    cout<<res<<endl;
    
    res = findNthDigit(108);
    cout<<res<<endl;
    */
    /*
    vector<int> input = {1, 3, 4, 6, 7, 9};
    for(int i = 0;i<100;i++)
    {
        int res = getRandom(input, 10);
        if(find(input.begin(), input.end(), res) != input.end())
        {
            cout<<res<<endl;
        }
    }
    cout<<"*************** done ***********************"<<endl;
     */
    /*
    TreeNode* n1 = new TreeNode(9);
    TreeNode* n2 = new TreeNode(6);
    TreeNode* n3 = new TreeNode(12);
    TreeNode* n4 = new TreeNode(3);
    TreeNode* n8 = new TreeNode(2);
    TreeNode* n5 = new TreeNode(7);
    TreeNode* n9 = new TreeNode(5);
    TreeNode* n6 = new TreeNode(10);
    TreeNode* n7 = new TreeNode(14);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n3->left = n6;
    n3->right = n7;
    n4->left = n8;
    n4->right = n9;
    */
    /*
    TreeNode* res = deleteNote(n1, n4);
    */
    /*
    TreeNode* n1 = new TreeNode(9);
    TreeNode* n2 = new TreeNode(6);
    TreeNode* n3 = new TreeNode(12);
    TreeNode* n4 = new TreeNode(3);
    TreeNode* n8 = new TreeNode(2);
    TreeNode* n5 = new TreeNode(7);
    TreeNode* n9 = new TreeNode(5);
    TreeNode* n6 = new TreeNode(10);
    TreeNode* n7 = new TreeNode(14);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n3->left = n6;
    n3->right = n7;
    n4->left = n8;
    n4->right = n9;
    printAllPath(n1);
    */
    
    /*
    char* a = read(5);
    char* b = read(7);
    for(int i = 0;i< 5;i++)
    {
        cout<<*(a++);
    }
    cout<<endl;
    for(int i = 0;i< 7;i++)
    {
        cout<<*(b++);
    }
    cout<<endl;
    */
    

    //printx(" char* a = read(5); char* b = read(7); for(int i = 0;i< 5;i++){ cout<<*(a++);} cout<<endl; /*for(int i = 0;i< 7;i++){ */
    //cout<<*(b++); /*}cout<</*e*/*/ndl;*/");
    
    
    TreeNode* n1 = new TreeNode(1);
    TreeNode* n2 = new TreeNode(2);
    TreeNode* n3 = new TreeNode(3);
    TreeNode* n4 = new TreeNode(4);
    TreeNode* n8 = new TreeNode(8);
    TreeNode* n5 = new TreeNode(5);
    TreeNode* n9 = new TreeNode(9);
    TreeNode* n6 = new TreeNode(6);
    TreeNode* n7 = new TreeNode(7);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n4->left = n6;
    n4->right = n7;
    n6->left = n8;
    n6->right = n9;
    
    conertBSTtoDoubleLinkedList(n1);
    
    /*
    string res = serialize2str("abc", "defg");
    deserialize(res);
    */
    /*
    int res;
    res = divideAgain(-2147483648,-1);
    cout<< res <<endl;
    res = divideAccurate(-2147483648,-1);
    cout<< res <<endl;
    res =(-2147483648)/(-1);
    cout<<res<<endl;
     */
    /*
    string res;
    res = FindShortestPattern("UAXXBAUBA", "AB");
    cout<<res<<endl;
     */
    
    /*
    bool res;
    res = isOneEditDistanceIII("a", "");
    cout<<res;
    */
    
    /*
    TreeNode* n1 = new TreeNode(9);
    TreeNode* n2 = new TreeNode(6);
    TreeNode* n3 = new TreeNode(12);
    TreeNode* n4 = new TreeNode(3);
    TreeNode* n8 = new TreeNode(2);
    TreeNode* n5 = new TreeNode(7);
    TreeNode* n9 = new TreeNode(5);
    TreeNode* n6 = new TreeNode(10);
    TreeNode* n7 = new TreeNode(14);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n3->left = n6;
    n3->right = n7;
    n4->left = n8;
    n4->right = n9;
    TreeNode* res;
    res = findSuccessor(n1, 13);
    cout<<res->val<<endl;
     */
    /*
    char res[100];
    mymemcpy2(res, "0123456789", 5);
    cout<<res<<endl;
     */
    /*
    int a[] = {8,4,0,3,0,2,0};
    moveZeroes(a, 7);
    for(int i = 0;i<7;i++) cout<<a[i]<<endl;
     */
    /*
    char* a = "aaa";
    char* b = "aaa";
    int res;
    res = strStrIII(a,b);
    cout<<res<<endl;
     */
    
    /*
    vector<vector<int>> input(4, vector<int>(4,1));
    cout<<findSubMatrix(input, 0,0,2,2)<<endl;
     */
    
    
    /*
    //0->1->2->3->4->5->6->7
    //k = 3
    //返回 2->1->0->3->4->5->7->6
    ListNode* n0 = new ListNode(0);
    ListNode* n1 = new ListNode(1);
    ListNode* n2 = new ListNode(2);
    ListNode* n3 = new ListNode(3);
    ListNode* n4 = new ListNode(4);
    ListNode* n5 = new ListNode(5);
    ListNode* n6 = new ListNode(6);
    ListNode* n7 = new ListNode(7);
    n0->next = n1;
    n1->next = n2;
    n2->next = n3;
    n3->next = n4;
    n4->next = n5;
    n5->next = n6;
    n6->next = n7;
    
    ListNode* res;
    res = reverseLinkedList(n0, 3);
    while(res)
    {
        cout<<res->val<<"->";
        res = res->next;
    }
     
     */
    
    /*
    double res;
    res = mysqrt(0.4);
    cout<<res<<endl;
    */
    
    /*
    TreeNode* n1 = new TreeNode(9);
    TreeNode* n2 = new TreeNode(6);
    TreeNode* n3 = new TreeNode(12);
    TreeNode* n4 = new TreeNode(3);
    TreeNode* n8 = new TreeNode(2);
    TreeNode* n5 = new TreeNode(7);
    TreeNode* n9 = new TreeNode(5);
    TreeNode* n6 = new TreeNode(10);
    TreeNode* n7 = new TreeNode(14);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n3->left = n6;
    n3->right = n7;
    n4->left = n8;
    n4->right = n9;
    
    //inorderTraversalTreeNode(n1);
    
    initializeTree(n1);
    while(haveNextTreeNode())
    {
        TreeNode* res = getNextNode();
        cout<<res->val<<endl;
    }
    */
    
    /*
    string str1 = "notnecessary to go to";
    string str2 = "kind of you not a";
    string res;
    res = longestCommonSubStr(str1, str2);
    cout<<res<<endl;
    
    int commonLength;
    commonLength = longestCommonSubString(str1, str2);
    cout<<commonLength<<endl;
     */
    
    /*
    vector<int> v1 = {1,2,3,4};
    vector<int> v2 = {5,6,7,8};
    vector<int> v3 = {9,10};
    vector<int> v4 = {};
        vector<int> v6 = {};
        vector<int> v7 = {};
    vector<int> v5 = {11, 12};
    
    vector<vector<int>> input = {v1,v2,v3,v4,v6,v7,v5};
    initializerVector(input);
    
    while( hasNextIterator())
    {
        auto item = getNextIterator();
        cout<<*item<<endl;
    }
    */
    
    /*
    int res;
    res = findHowManyWays(3,3);
    cout<<res<<endl;
    res = findHowManyWaysOneArray(3,3);
    cout<<res<<endl;
    */
    
    /*
    vector<int> input = {1,1,1,1,1,0,0,0,0};
    int res;
    res = findChange(input);
    cout<<res<<endl;
    
    input = {0,1};
    res = findChange(input);
    cout<<res<<endl;
    
    input = {0,0,0,0};
    res = findChange(input);
    cout<<res<<endl;
    */
    /*
    TreeNode* n1 = new TreeNode(9);
    TreeNode* n2 = new TreeNode(6);
    TreeNode* n3 = new TreeNode(12);
    TreeNode* n4 = new TreeNode(3);
    TreeNode* n8 = new TreeNode(2);
    TreeNode* n5 = new TreeNode(7);
    TreeNode* n9 = new TreeNode(5);
    TreeNode* n6 = new TreeNode(10);
    TreeNode* n7 = new TreeNode(14);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n3->left = n6;
    n3->right = n7;
    n4->left = n8;
    n4->right = n9;
    
    int res;
    res = printNth(n1, 1);
    cout<<res<<endl;
    res = printNth(n1, 2);
    cout<<res<<endl;
    res = printNth(n1, 3);
    cout<<res<<endl;
    res = printNth(n1, 4);
    cout<<res<<endl;
    res = printNth(n1, 5);
    cout<<res<<endl;
    res = printNth(n1, 6);
    cout<<res<<endl;
    res = printNth(n1, 7);
    cout<<res<<endl;
    res = printNth(n1, 8);
    cout<<res<<endl;
    res = printNth(n1, 9);
    cout<<res<<endl;
    res = printNth(n1, 10);
    cout<<res<<endl;
    */
    
    /*
    TreeNode* n1 = new TreeNode(9);
    TreeNode* n2 = new TreeNode(5);
    TreeNode* n3 = new TreeNode(3);
    TreeNode* n4 = new TreeNode(3);
    TreeNode* n8 = new TreeNode(1);
    TreeNode* n5 = new TreeNode(1);
    TreeNode* n9 = new TreeNode(1);
    TreeNode* n6 = new TreeNode(1);
    TreeNode* n7 = new TreeNode(1);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n3->left = n6;
    n3->right = n7;
    n4->left = n8;
    n4->right = n9;
    
    int res;
    res = printNthWithNumber(n1, 1);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 2);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 3);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 4);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 5);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 6);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 7);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 8);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 9);
    cout<<res<<endl;
    res = printNthWithNumber(n1, 10);
    cout<<res<<endl;
    */
    
    /*
     
    TreeNode* n1 = new TreeNode(9);
    TreeNode* n2 = new TreeNode(5);
    TreeNode* n3 = new TreeNode(3);
    TreeNode* n4 = new TreeNode(3);
    TreeNode* n8 = new TreeNode(1);
    TreeNode* n5 = new TreeNode(1);
    TreeNode* n9 = new TreeNode(1);
    TreeNode* n6 = new TreeNode(1);
    TreeNode* n7 = new TreeNode(1);
    
    n1->left = n2;
    n1->right = n3;
    n2->left = n4;
    n2->right = n5;
    n3->left = n6;
    n3->right = n7;
    n4->left = n8;
    n4->right = n9;
    
    cout<<"sum = 18 :"<<endl;
    FindAllPath(n1,18);
    cout<<"sum = 15 :"<<endl;
    FindAllPath(n1,15);
    cout<<"sum = 13 :"<<endl;
    FindAllPath(n1,13);
    */
    
    /*
    bool res;
    res = isPalindromeX("abcba");
    cout<<res<<endl;
    res = isPalindromeX("aB   Cb   A  ");
    cout<<res<<endl;
    res = isPalindromeX("a&b&c#b!a");
    cout<<res<<endl;
    res = isPalindromeX("axbcba");
    cout<<res<<endl;
     */
    
    /*
    vector<int> input = {2,1};
    int res = findMinX(input);
    cout<<res<<endl;
    */
    
    /*
    ListNode* l1 = new ListNode(1);
    ListNode* l2 = new ListNode(2);
    ListNode* l3 = new ListNode(3);
    ListNode* l4 = new ListNode(4);
    ListNode* l5 = new ListNode(5);
    ListNode* l6 = new ListNode(6);
    
    l1->next = l2;
    l2->next = l3;
    l3->next = l4;
    l4->next = l5;
    l5->next = l6;
    l6->next = l1;
    
    insertCircularList(l1, 3);
    insertCircularList(l4, 10);
    insertCircularList(l4, 7);
    */
    
    /*
    addWord("rat");
    addWord("cat");
    addWord("cat");
    addWord("bat");
    cout<<searchWord("dat")<<endl;
    cout<<searchWord("bat")<<endl;
    cout<<searchWord(".a.")<<endl;
    cout<<searchWord("r.t")<<endl;
    cout<<searchWord("t.t")<<endl;
    cout<<searchWord("...")<<endl;
    */
    
    /*
    addWordTrie("rat");
    addWordTrie("cat");
    addWordTrie("cat");
    addWordTrie("bat");
    cout<<searchWordTrie("dat")<<endl;
    cout<<searchWordTrie("bat")<<endl;
    cout<<searchWordTrie(".a.")<<endl;
    cout<<searchWordTrie("c..")<<endl;
    cout<<searchWordTrie("r.t")<<endl;
    cout<<searchWordTrie("t.t")<<endl;
    cout<<searchWordTrie("...")<<endl;
    */
    
    /*
    addWordTrie("rat");
    addWordTrie("cat");
    addWordTrie("cat");
    addWordTrie("bat");
    cout<<searchWordTrieDFS("dat")<<endl;
    cout<<searchWordTrieDFS("bat")<<endl;
    cout<<searchWordTrieDFS(".a.")<<endl;
    cout<<searchWordTrieDFS("c..")<<endl;
    cout<<searchWordTrieDFS("r.t")<<endl;
    cout<<searchWordTrieDFS("t.t")<<endl;
    cout<<searchWordTrieDFS("...")<<endl;
    */
    
    /*
    vector<string> res;
    res = letterCombinationsx("2");
    for(auto s: res)
    {
        cout<<s<<endl;
    }
    */
    
    /*
    vector<string> input = {"abc", "ba","cb", "aaaa", "cba"};
    vector<pair<string, string>> res;
    res = FindAllPalindromePairs(input);
    for(auto item: res)
    {
        cout<<item.first<<":  "<<item.second<<endl;
    }
    */
    
    /*
    vector<char> input= {'a', 's', 'b', 'b', 'b', 'b', 'a', 'a'};
    vector<char> res;
    res = modify(input);
    for(auto c: res)
    {
        cout<<c<<" ";
    }
    cout<<endl;
    */
    
    /*
    thread t1(func1);
    thread t2(func2);
    t1.join();
    t2.join();
     */
    
    /*
    vector<string> res;
    res = generateAllParentheses(3);
    for(auto s: res)
    {
        cout<<s<<endl;
    }
     */
    /*
    cout<<divideIIII(2147483647,1)<<endl;
     */
    
    /*
    string res;
    res = FindMinwindowsSubStr("a", "aa");
    cout<<res;
    */
    /*
    //vector<int> input = {1, 2, 3, 8, 10, 5, 6, 7, 12, 9, 11, 4, 0};
    vector<int> input = {1,2,3,8,10,4,5,6,7,12,9};
    int res;
    res = MakeAscending(input);
    cout<<res<<endl;
    */
    
    /*
    vector<int> input1 = {2, 3, 5, 8, 13};
    vector<int> input2 = {4, 8, 12, 16};
    
    vector<int> res;
    res = FindMaxSum(input1, input2, 5);
    for(auto i: res)
    {
        cout<<i<<" ";
    }
    cout<<endl;
     */
    
    /*
    ListNodeWithChild* n1 = new ListNodeWithChild(1);
    ListNodeWithChild* n2 = new ListNodeWithChild(2);
    ListNodeWithChild* n3 = new ListNodeWithChild(3);
    ListNodeWithChild* n4 = new ListNodeWithChild(4);
    ListNodeWithChild* n5 = new ListNodeWithChild(5);
    ListNodeWithChild* n6 = new ListNodeWithChild(6);
    ListNodeWithChild* n7 = new ListNodeWithChild(7);
    ListNodeWithChild* n8 = new ListNodeWithChild(8);

    n1->next = n2;
    n2->next = n3;
    n3->next = n7;
    n7->next = n8;
    n4->next = n5;
    n5->next = n6;
    n3->child = n4;
    
    ListNodeWithChild* res;
    res = flatten(n1);
    while(res)
    {
        cout<<res->val<<" ";
        res= res->next;
    }
    cout<<endl;
    
    res = unflatten(n1);
    while(res)
    {
        cout<<res->val<<" ";
                res= res->next;
    }
    cout<<endl;
    */
    
    /*
    vector<vector<int>> res ;
    vector<int> input = {1,2,2};
    res = subsetsWithDupx(input);
    res = subsetsWithDup(input);
    res = subsetsWithDupII(input);
    */
    /*
    string str1 = "notnecessary to go to";
    string str2 = "kind of you not a";
    string res;
    res = longestCommonSubStr(str1, str2);
    cout<<res<<endl;
    
     
    int commonLength;
    commonLength = longestCommonSubString(str1, str2);
    cout<<commonLength<<endl;
    */
    
    /*
    bool res;
    vector<int> input = {1,0,1,0,1};
    res = canJump(input);
    cout<<res<<endl;
    
    int a[] = {1,1,0,0,0,0,0,100};
    res = canJump(a, 8);
    cout<<res<<endl;
    */
    
    /*
    TreeNode* n1 = new TreeNode(100);
    TreeNode* n2 = new TreeNode(2);
    TreeNode* n3 = new TreeNode(300);
    TreeNode* n4 = new TreeNode(64);
    TreeNode* n5 = new TreeNode(5);
    TreeNode* n6 = new TreeNode(60);
    TreeNode* n7 = new TreeNode(7);
    TreeNode* n8 = new TreeNode(8);
    TreeNode* n9 = new TreeNode(9);
    TreeNode* n10 = new TreeNode(10);
    
    n1->left = n2;
    n1->right = n3;
    n3->left = n10;
    n2->left = n4;
    n2->right = n5;
    n4->left = n6;
    n4->right = n7;
    n6->left = n8;
    n6->right = n9;
    int res = INT_MIN;
    TreeNode* root = LCA(n1, n9, n10);
    findMaxB(root, n9, n10,res);
    cout<<res<<endl;
    res = 0;
    findMaxC(n1, n9, n10, res);
    cout<<res<<endl;
    
    
    root = LCA(n1, n9, n5);
    findMaxB(root, n9, n5,res);
    cout<<res<<endl;
    res = 0;
    findMaxC(n1, n9, n5, res);
    cout<<res<<endl;
    */
    /*
    vector<vector<int>> res;
    res = allSubset(3);
    cout<<res.size()<<endl;
     */
    return 0;
}