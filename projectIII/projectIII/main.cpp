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
#include <unordered_map>
#include <stack>
#include <assert.h>
#include <numeric>
#include <unordered_set>
#include <queue>
#include <math.h>
#include <stdlib.h>
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
    TrieNode()
    {
        fill_n(&count[0], 256, 0);
        fill_n(&val[0], 256, nullptr);
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
void workerx(vector<int> &S, int level, vector<int>& tmp, vector<vector<int>>& res)
{
    res.push_back(tmp);
    for(int i = level;i< S.size();i++)
    {
        if(i == level || S[i] != S[i-1])
        {
            tmp.push_back(S[i]);
            workerx(S, level+1, tmp, res);
            tmp.pop_back();
        }
    }
}
vector<vector<int> > subsetsWithDup(vector<int> &S) {
    vector<vector<int>> res;
    if(S.size() == 0) return res;
    
    sort(S.begin(), S.end());
    vector<int> tmp;
    workerx(S, 0, tmp, res);
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
        if(!(tracker.find(num[i])->second))
        {
            tracker.find(num[i])->second = true;
            int tmp = 1;
            
            int runner = num[i]+1;
            while(tracker.find(runner) != tracker.end())
            {
                tmp ++;
                tracker.find(runner)->second = true;
                runner++;
            }
            runner = num[i] - 1;
            while(tracker.find(runner) != tracker.end())
            {
                tmp ++;
                tracker.find(runner)->second = true;
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

int findKthElement(int A[], int m, int B[], int n, int k)
{
    if(m > n)
    {
        return findKthElement(B, n, A, m, k);
    }
    if(m == 0)
        return B[k-1];
    if(k ==1)
        return min(A[0], B[0]);
    
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

//2.1.5
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
        if(cur->left &&
           cur->left != prev &&
           (!(cur->right) || cur->right != prev))
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
    for(int i = level;i<=n;i++)
    {
        tracker.push_back(i);
        combineWorker(n, k, level+1, res, tracker);
        tracker.pop_back();
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
    int n = queens.size();
    
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
int sqrt(int x) {
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
        lastIndex = max(lastIndex, i+ A[i]);
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
        return numDecodings(s.substr(1, s.length()-1));
    }
    else
    {
        if(s[0] == '1' || (s[0] =='2' && s[1] <= '6'))
        {
            return numDecodings(s.substr(1, s.length()-1)) +
            numDecodings(s.substr(2, s.length()-2));
        }
        else
        {
            return numDecodings(s.substr(1, s.length()-1));
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
    for (size_t i = 1; i <= s.length(); ++i) {
        for (int j = i - 1; j >= 0; --j) {
            if (f[j] && dict.find(s.substr(j, i - j)) != dict.end()) {
                f[i] = true;
                prev[i][j] = true;
            }
        } }
    vector<string> result;
    vector<string> path;
    gen_path(s, prev, s.length(), path, result);
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
    int right = intervals.size()-1;
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
    right = intervals.size()-1;
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
            result.push_back(distance(begin(s), i));
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
    vector<int> input = {2,3,1};
    int res = findMin(input);
    cout<<res<<endl;
    
    // insert code here...
    //std::cout << "Hello, World!\n";
    return 0;
}

