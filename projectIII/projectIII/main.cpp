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
void worker(vector<int>& num, int target, vector<int>& tmp, int level, vector<vector<int>>& res)
{
    if(0 == target)
    {
        res.push_back(tmp);
        return;
    }
    else
    {
        int prev = -1;
        if(level < num.size())
        {
            for(int i = level+1;i< num.size();i++)
            {
                if(num[i] == prev)
                {
                    continue;
                }
                prev = num[i];
                tmp.push_back(num[i]);
                worker(num, target - num[i], tmp, i+1, res);
                tmp.pop_back();
            }
        }
    }
}
vector<vector<int> > combinationSum2(vector<int> &num, int target) {
    
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
    
    while(left < right)
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


int main(int argc, const char * argv[])
{
    
    int A[] = {1};
    int B[] = {2,3,4,5,6,7};
    
    double res = findMedianSortedArrays(A, 1, B, 6);
    cout<<res<<endl;

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
    // insert code here...
    //std::cout << "Hello, World!\n";
    return 0;
}

