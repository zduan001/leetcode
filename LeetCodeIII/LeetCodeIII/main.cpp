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

int maximalRectangle(vector<vector<char> > &matrix) {
    if(matrix.size() == 0 || matrix[0].size() == 0) return 0;
    vector<int> tracker(matrix[0].size()+1, 0);
    int maxArea = 0;
    for(int i = 0;i< (int)matrix.size();i++)
    {
        for(int j = 0;j< (int)matrix[i].size();j++)
        {
            if(matrix[i][j] == '1'){
                tracker[j] ++;
            }
            else
            {
                tracker[j] = 0;
            }
        }
        maxArea = max(maxArea, largestRectangleArea(tracker));
    }
    return maxArea;
}

bool isInterleave(string s1, string s2, string s3)
{
    if(s1.length() + s2.length() != s3.length()) return false;
    vector<vector<bool>> tracker(s1.length() + 1, vector<bool>(s2.length()+1, true));
    
    for(int i = 1;i<= s1.length();i++)
    {
        tracker[i][0] = s1[i-1] == s3[i-1] && tracker[i-1][0];
    }
    
    for(int j = 1;j<=s2.length();j++)
    {
        tracker[0][j] = s2[j-1] == s3[j-1] && tracker[0][j-1];
    }
    
    for(int i = 1;i<= s1.length();i++)
    {
        for(int j = 1;j<=s2.length();j++)
        {
            tracker[i][j] = (s1[i-1] == s3[i+j-1] && tracker[i-1][j]) || (s2[j-1] == s3[i+j-1] && tracker[i][j-1]);
        }
    }
    return tracker[s1.length()][s2.length()];
}

string longestPalindrome(string s) {
    int maxlength = 0;
    string res;
    for(int i = 0;i< s.length();i++)
    {
        int left = i;
        int right = i;
        while(left>=0 && right<s.length() && s[left] == s[right])
        {
            int length = right - left +1;
            if(length > maxlength)
            {
                maxlength = length;
                res = s.substr(left, length);
            }
            left--;right++;
        }
        left = i;
        right = i+1;
        while(left>=0 && right<s.length() && s[left] == s[right])
        {
            int length = right - left +1;
            if(length > maxlength)
            {
                maxlength = length;
                res = s.substr(left, length);
            }
            left--;right++;
        }
    }
    return res;
}

ListNode* reverse(ListNode* head)
{
    ListNode* newHead= NULL;
    ListNode* tmp = NULL;
    while(head)
    {
        tmp = head;
        head = head->next;
        tmp->next = newHead;
        newHead = tmp;
    }
    return newHead;
}

ListNode *reverseKGroup(ListNode *head, int k)
{
    ListNode* dummyHead = new ListNode(-1);
    dummyHead->next = head;
    
    ListNode* front = dummyHead;
    ListNode* tmp = dummyHead;
    
    while(front->next)
    {
        for(int i = 0;i< k;i++)
        {
            tmp = tmp->next;
            if(!tmp)
            {
                return dummyHead->next;
            }
        }
        
        head = front->next;
        ListNode* nextHead = tmp->next;
        tmp->next = NULL;
        
        front->next = reverse(head);
        head->next = nextHead;
        front = head;
        tmp = front;
        
    }
    return dummyHead->next;
}

char *strStr(const char *haystack, const char *needle) {
    int n = (int)strlen(needle);
    int m = (int) strlen(haystack);
    if(m<n) return NULL;
    
    int next[n];
    next[0] = -1;
    int j = next[0];
    for(int i = 1;i< n;i++)
    {
        while(j != -1 && needle[j+1] != needle[i])
            j = next[j];
        if(needle[j+1] == needle[i])
            j++;
        next[i] = j;
    }
    
    j = -1;
    int i = 0;
    while(i<m)
    {
        while(j!= -1 && needle[j+1] != haystack[i])
        {
            j = next[j];
        }
        
        if(needle[j+1] == haystack[i]) j++;

        if(j == n-1)
        {
            return (char*)haystack+i-j;
        }
        i++;
    }
    
    return NULL;
}

int divide(int dividend, int divisor) {
    long long dividendlong = dividend;
    long long divisorlong = divisor;
    int sign = 1;
    if (dividendlong < 0)
    {
        sign *= -1;
        dividendlong *= -1;
    }
    if (divisorlong < 0)
    {
        sign *= -1;
        divisorlong *= -1;
    }
    
    if(divisorlong == 0) return -1;
    int res = 0;
    int i = 0;
    while(divisorlong <= dividendlong)
    {
        i = 0;
        long tmp = divisorlong;
        while(tmp <= dividendlong)
        {
            res += 1 << i;
            dividendlong -= tmp;
            i++;
            tmp = tmp<< 1;
        }
    }
    return res * sign;
}

int longestValidParentheses(string s) {
    int last = -1;
    stack<int> left;
    int maxLength = 0;
    for(int i = 0;i<s.length();i++)
    {
        if(s[i] == '(')
        {
            left.push(i);
        }
        else
        {
            if(left.empty())
            {
                last = i;
            }
            else
            {
                left.pop();
                if(left.empty())
                {
                    maxLength = max(maxLength, i - last);
                }
                else
                {
                    maxLength = max(maxLength, i- left.top());
                }
            }
        }
    }
    return maxLength;
}


int searchI(int A[], int n, int target) {
    int mid;
    int left = 0;
    int right = n-1;
    while(left <=right)
    {
        mid = left + (right-left)/2;
        if(A[mid] == target)
        {
            return mid;
        }
        else if(A[left] < A[mid])
        {
            if(target >= A[left] && target < A[mid])
            {
                right = mid-1;
            }
            else
            {
                left = mid +1;
            }
        }
        else if(A[left] > A[mid])
        {
            if(target > A[mid] && target <= A[right])
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
            left++;
        }
    }
    return -1;
}

vector<int> searchRange(int A[], int n, int target) {
    int left = 0;
    int right = n-1;
    int mid;
    while(left<=right)
    {
        mid = left + (right-left)/2;
        if(A[mid] == target)
        {
            break;
        }
        else if(A[mid] > target)
        {
            right = mid -1;
        }
        else if(A[mid] < target)
        {
            left = mid +1;
        }
    }
    
    if(left > right) return vector<int>(2, -1);
    int leftindex = -1;
    int rightindex = -1;
    
    int ll = left;
    int rr = mid;
    while(ll <=rr)
    {
        int midl = ll + (rr-ll) /2;
        if(A[midl] == target)
        {
            if(midl == 0 || A[midl-1] < A[midl])
            {
                leftindex = midl;
                break;
            }
            else
            {
                rr = midl-1;
            }
        }
        else if(A[midl] < target)
        {
            ll = midl+1;
        }
        else
        {
            rr = midl-11;
        }
    }
    
    ll = mid;
    rr = right;
    while(ll<=rr)
    {
        int midl = ll + (rr-ll) /2;
        if(A[midl] == target)
        {
            if(midl == n-1 || A[midl] < A[midl+1])
            {
                rightindex = midl;
                break;
            }
            else
            {
                ll = midl+1;
            }
        }
        else if(A[midl] < target)
        {
            ll = midl+1;
        }
        else
        {
            rr = midl-1;
        }
    }
    
    vector<int> res;
    res.push_back(leftindex);
    res.push_back(rightindex);
    return res;
    
}

void worker(vector<int>& num, int target, vector<int>& tmp, int level, vector<vector<int>>& res)
{
    if(0 == target)
    {
        res.push_back(tmp);
        return;
    } else if(0 > target)
    {
        return;
    }else if(level == num.size())
    {
        return;
    }
    else
    {
        int prev = -1;
        for(int i = level;i< num.size(); i++)
        {
            if(prev == num[i]) continue;
            prev = num[i];
            tmp.push_back(num[i]);
            worker(num, target-num[i], tmp, i+1, res);
            tmp.pop_back();
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

int trap(int A[], int n) {
    if(n<2) return 0;
    int leftMax = A[0];
    int rightMax = A[n-1];
    int maxWater = 0;
    int left = 0;
    int right = n-1;
    while(left<right)
    {
        if(A[left]< A[right])
        {
            if(A[left] > leftMax)
            {
                leftMax = A[left];
            }
            else
            {
                maxWater += leftMax - A[left];
            }
            left ++;
        }
        else
        {
           if(A[right]> rightMax)
           {
               rightMax = A[right];
           }
           else
           {
               maxWater += rightMax- A[right];
               right --;
           }
    
        }
    }
    return maxWater;
}

string multiply(string num1, string num2) {
    // check num1 and num2 to make sure it
    // meet the right
    if(num1 == "0" || num2 == "0") return "0";
    int m = (int)num1.length();
    int n = (int)num2.length();
    vector<int> res(m+n, 0);
    vector<char> charnum1(num1.begin(), num1.end());
    vector<char> charnum2(num2.begin(), num2.end());
    
    reverse(charnum1.begin(), charnum1.end());
    reverse(charnum2.begin(), charnum2.end());
    
    int carry = 0;
    int tmp = 0;
    for(int i = 0;i<m;i++)
    {
        for(int j = 0;j<n;j++)
        {
            tmp = (charnum1[i] - '0') * (charnum2[j] - '0') + carry;
            tmp += res[i+j];
            carry = tmp /10;
            tmp = tmp % 10;
            res[i+j] = tmp;
        }
        res[i+n] += carry;
        carry = 0;
    }
    res[n+m-1] += carry;
    
    if(res[n+m-1] == 0) res.pop_back();
    
    reverse(res.begin(), res.end());
    string s;
    for(int i = 0;i<res.size();i++)
    {
        s.push_back((char)(res[i] + '0'));
    }
    return s;
    
}

bool canJump(int A[], int n) {
    unordered_set<int> visited;
    queue<int> q;
    q.push(n-1);
    while(!q.empty())
    {
        int tmp = q.front();
        q.pop();
        for(int i = 0;i<= tmp;i++)
        {
            if(A[i]>= tmp-i)
            {
                if(i == 0) return true;
                if(visited.find(i) == visited.end())
                {
                    q.push(i);
                    visited.insert(i);
                }
            }
        }
    }
    return false;
    
}

int jump(int A[], int n) {
    if(n == 1) return 0;
    unordered_set<int> visited;
    queue<int> q;
    int steps = 0;
    q.push(n-1);
    q.push(-1);
    while(!q.empty())
    {
        int tmp = q.front();
        q.pop();
        if(tmp == -1)
        {
            steps++;
            if(q.empty())
            {
                break;
            }
            
            q.push(-1);
            continue;
        }
        
        for(int i = 0;i<= tmp;i++)
        {
            if(A[i]>= tmp-i)
            {
                if(i == 0) return steps+1;
                if(visited.find(i) == visited.end())
                {
                    q.push(i);
                    visited.insert(i);
                }
            }
        }
    }
    return steps;
}
/*
void swap(vector<int>& num, int i, int j)
{
    int tmp = num[i];
    num[i] = num[j];
    num[j] = tmp;
}
*/
void worker(vector<int>& num, int level, vector<vector<int>>& res)
{
    if(level == num.size())
    {
        res.push_back(num);
        return;
    }
    else
    {
        unordered_set<int> set;
        for(int i = level;i<num.size(); i++)
        {
            if(set.find(num[i]) == set.end())
            {
                swap(num, level, i);
                worker(num, level+1, res);
                swap(num, level, i);
                set.insert(num[i]);
            }
        }
    }
}

vector<vector<int> > permuteUnique(vector<int> &num) {
    vector<vector<int>> res;
    sort(num.begin(), num.end());
    worker(num, 0, res);
    return res;
}

void rotate(vector<vector<int>> &matrix) {
    
    int n = (int)matrix.size();
    if (n <=1) return;
    
    for(int i = 0;i< n;i++)
    {
        for(int j = i;j< n-i;j++)
        {
            int tmp = matrix[i][j];
            matrix[i][j] = matrix[n-j][i];
            matrix[n-j][i] = matrix[n-i][n-j];
            matrix[n-i][n-j] = matrix[j][n-i];
            matrix[j][n-i] = tmp;
            
            
            //matrix[i][j] = matrix[j][n-i];
            //matrix[j][n-i] = matrix[n-i][n-j];
            //matrix[n-i][n-j] = matrix[n-j][i];
            //matrix[n-j][i] = tmp;
        }
    }
}

bool check(vector<vector<char>> tmp, int n)
{
    for(int i = 0;i< n;i++)
    {
        int count = 0;
        for(int j = 0;j<n;j++)
        {
            if(tmp[i][j] == 'Q')
                count++;
        }
        if(count >1) return false;
    }
    
    for(int j = 0;j<n;j++)
    {
        int count = 0;
        for(int i = 0;i<n;i++)
        {
            if(tmp[i][j] == 'Q')
                count ++;
        }
        if(count>1) return false;
    }
    
    for(int i = 0;i<n;i++)
    {
        int x = i;
        int y = -1;
        
        for(int j = 0;j<n;j++)
        {
            if(tmp[i][j] == 'Q')
            {
                y = j;
                break;
            }
        }
        
        if(x>=0 && y >=0)
        {
            int i = 1;
            while(x+i < n && y +i <n)
            {
                if(tmp[x+i][y+i] == 'Q')
                    return false;
                i++;
            }
            
            i = 1;
            while(x-i >=0 && y-i >= 0)
            {
                if(tmp[x-i][y-i] == 'Q')
                    return false;
                i++;
            }
        }
    }
    return true;
}

void worker(vector<vector<char>> tmp, int level, vector<vector<string>>& res, int n)
{
    if(level == n )
    {
        vector<string> xyz;
        for(int i = 0;i< (int)tmp.size(); i++)
        {
            string s(tmp[0].begin(), tmp[0].end());
            xyz.push_back(s);
            
        }
        res.push_back(xyz);
    }
    else
    {
        for(int i = 0;i<n;i++)
        {
            tmp[level][i] = 'Q';
            if (check(tmp,n))
            {
                worker(tmp, level+1, res, n);
            }
            tmp[level][i] = '.';
        }
    }
}

vector<vector<string> > solveNQueens(int n) {
    vector<vector<char>> tmp(n, vector<char>(n, '.'));
    vector<vector<string>> res;
    worker(tmp, 0, res, n);
    return res;
}

vector<int> spiralOrder(vector<vector<int> > &matrix) {
    vector<int> res;
    int m = (int)matrix.size();
    if(m == 0) return res;
    int n = (int)matrix[0].size();
    if(n ==0) return res;
    int top = 0;
    int right = n-1;
    int bottom = m-1;
    int left = 0;
    
    int direction = 0;
    
    int i = 0;
    int j = 0;
    
    while(i>=top && i<=bottom && j>= left && j <= right)
    {
        res.push_back(matrix[i][j]);
        if(direction == 0)
        {
            j++;
            if(j >right)
            {
                j = right;
                i++;
                direction ++;
                top ++;
            }
        }
        else if(direction == 1)
        {
            i++;
            if(i>bottom)
            {
                i = bottom;
                j--;
                direction ++;
                right --;
            }
        }
        else if(direction ==2)
        {
            j--;
            if(j<left)
            {
                j = left;
                i--;
                direction ++;
                bottom --;
            }
        }
        else{
            i--;
            if(i < top)
            {
                i = top;
                j++;
                direction = 0;
                left ++;
            }
        }
    }
    return res;
}

bool sortIntervalFunctions(Interval i1, Interval i2)
{
    return (i1.start < i2.start);
}

vector<Interval> merge(vector<Interval> &intervals) {
    sort(intervals.begin(), intervals.end(), sortIntervalFunctions);
    vector<Interval> res;
    if (intervals.size() == 0) return res;
    
    Interval current = intervals[0];
    for(int i = 1;i<(int)intervals.size(); i++)
    {
        if(current.end < intervals[i].start)
        {
            res.push_back(current);
            current = intervals[i];
        }
        else if(current.end < intervals[i].end)
        {
            current.end = intervals[i].end;
        }
    }
    return res;
}

vector<Interval> insert(vector<Interval> &intervals, Interval newInterval) {
    vector<Interval> res;
    int n = (int)intervals.size();
    if(n == 0)
    {
        res.push_back(newInterval);
        return res;
    }
    int i = 0;
    while(i<n && newInterval.start > intervals[i].end)
    {
        res.push_back(intervals[i++]);
    }
    
    while(i<n && newInterval.end <= intervals[i].start)
    {
        newInterval.start = min(newInterval.start, intervals[i].start);
        newInterval.end = max(newInterval.end, intervals[i].end);
        i++;
    }
    res.push_back(newInterval);
    
    while(i<n)
        res.push_back(intervals[i++]);
    
    return res;
    
}

vector<string> fullJustify(vector<string> &words, int L) {
    vector<string> res;
    if(words.size() == 0) return res;
    int i = 0;
    vector<string> tmp;
    int totalWordLength =0;
    int totalCount =0;
    
    while(i< (int)words.size())
    {
        if(totalWordLength + totalCount + words[i].length() <=L)
        {
            totalWordLength += words[i].length();
            tmp.push_back(words[i++]);
            totalCount ++;
        }
        else
        {
            if(tmp.size() >1 )
            {
                int space = (L-totalWordLength)/(totalCount-1);
                int extra = L - (space * (totalCount-1)) - totalWordLength;
                string xyz;
                
                for(int k = 0;k<tmp.size()-1;k++)
                {
                    xyz +=tmp[k];
                    for(int j = 0;j<space;j++)
                        xyz += " ";
                    if(extra >0)
                    {
                        xyz+=" ";
                        extra--;
                    }
                }
                xyz += tmp[tmp.size() -1];
                res.push_back(xyz);
            }
            else
            {
                if(tmp.size() == 1)
                {
                    // although there is no need for this padding white space
                    while(tmp[0].size() < L)
                        tmp[0].append(" ");
                    res.push_back(tmp[0]);
                }
            }
            tmp.clear();
            totalCount = 0;
            totalWordLength = 0;
        }
    }
    
    if(tmp.size() > 0){
        string xyz;
        for(int i = 0;i<tmp.size();i++)
        {
            xyz += tmp[i];
            if(i != tmp.size()-1)
            {
                xyz += ' ';
            }
        }
        while(xyz.length() < L)
            xyz.append(" ");
        res.push_back(xyz);
    }
    return res;
}

int sqrt(int x) {
    if(x == 0) return 0;
    if(x == 1) return 1;
    double right = x;
    double left = 0;
    double mid;
    while(abs(right-left) > 0.00000001){
        mid = left + (right-left) / 2;
        if(mid * mid > x){
            right = mid;
        } else if(mid * mid < x){
            left = mid;
        } else {
            return mid;
        }
    }
    
    return right;
    
}

int minDistance(string word1, string word2) {
    int m = (int)word1.length();
    int n = (int)word2.length();
    
    if(m== 0 && n == 0) return 0;
    if(m == 0) return n;
    if(n == 0) return m;
    
    vector<vector<int>> tracker(m+1, vector<int>(n+1));
    
    for(int i = 0;i<n;i++)
    {
        tracker[0][i] = i;
    }
    
    for(int i = 0;i<m;i++)
    {
        tracker[i][0] = i;
    }
    
    for(int i = 1;i<=m;i++)
    {
        for(int j = 1;j<=n;j++)
        {
            tracker[i][j] = min(tracker[i-1][j] +1, tracker[i][j-1]+1);
            int tmp = word1[i-1] == word2[j-1]? tracker[i-1][j-1]: tracker[i-1][j-1] +1;
            tracker[i][j] = min (tracker[i][j] , tmp);
        }
    }
    return tracker[m-1][n-1];
}

void sortColors(int A[], int n)
{
    int i = -1;
    int left = 0;
    int right = n-1;
    while(left<=right)
    {
        if(A[left] ==0)
        {
            swap(A, ++i, left++);
        }
        else if(A[left] == 1)
        {
            left++;
        }
        else if(A[left] == 2)
        {
            swap(A, left, right--);
        }
    }
}

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
        if (TAlphaStat[tc] > 0)
        {
            if (currentAlphaStat[tc] < TAlphaStat[tc])
            {
                containNum++;
            }
            currentAlphaStat[tc]++;
        }
        if (containNum == T.length())
        {
            while (head < tail && (currentAlphaStat[S[head]] > TAlphaStat[S[head]] || TAlphaStat[S[head]] == 0))
            {
                currentAlphaStat[S[head]]--;
                head++;
            }
            if (minWidth > tail - head + 1)
            {
                minWidth = tail - head + 1;
                left = head;
            }
        }
        tail++;
    }
    return left == -1 ? "" : S.substr(left, minWidth);
}

void worker(vector<int> &s, int begin, vector<int> &path, vector<vector<int> > &result) {
    result.push_back(path);
    
    for (int i = begin; i < s.size(); ++i) {
        if (i == begin || s[i] != s[i - 1])
        {
            path.push_back(s[i]);
            worker(s, i + 1, path, result);
            path.pop_back();
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

vector<int> inorderTraversal(TreeNode *root) {
    vector<int> res;
    if(!root) return res;
    stack<TreeNode*> tracker;
    bool backtrack = false;
    tracker.push(root);
    while(!tracker.empty())
    {
        TreeNode* tmp = tracker.top();
        
        if(tmp->left && !backtrack)
        {
            tracker.push(tmp->left);
            continue;
        }
        
        tracker.pop();
        res.push_back(tmp->val);
        backtrack = true;
        
        if( tmp->right)
        {
            tracker.push(tmp->right);
            backtrack = false;
        }
    }
    return res;
    
}

int worker(int left, int right)
{
    if(left >= right) return 1;
    int res = 0;
    
    for(int i = left;i<=right;i++)
    {
        int leftres = worker(left, i-1);
        int rightres = worker(i+1, right);
        res += leftres * rightres;
    }
    return res;
}

int numTrees(int n) {
    return worker(1,n);
}


vector<vector<int> > zigzagLevelOrder(TreeNode *root) {
    stack<TreeNode *> s1;
    stack<TreeNode *> s2;
    vector<vector<int>> res;
    if (!root) return res;
    s1.push(root);
    bool leftToRight = true;
    vector<int>* container  = new vector<int>();
    while(!s1.empty())
    {
        TreeNode* tmp = s1.top();
        s1.pop();
        container->push_back(tmp->val);
        if(leftToRight){
            if(tmp->left) s2.push(tmp->left);
            if(tmp->right) s2.push(tmp->right);
        }else{
            if(tmp->right) s2.push(tmp->right);
            if(tmp->left) s2.push(tmp->left);
        }
        
        if(s1.empty())
        {
            leftToRight = !leftToRight;
            swap(s1, s2);
            res.push_back(*container);
            container = new vector<int>();
        }
    }
    return res;
}


vector<vector<int> > zigzagLevelOrderI(TreeNode *root) {
    vector<vector<int>> res;
    if(!root) return res;
    stack<TreeNode*> s1;
    stack<TreeNode*> s2;
    bool leftToRight = true;
    s1.push(root);
    vector<int>* container = new vector<int>();
    while(!s1.empty()){
        TreeNode* tmp = s1.top();
        container->push_back(tmp->val);
        if(leftToRight){
            if(tmp->left) s2.push(tmp->left);
            if(tmp->right) s2.push(tmp->right);
        }else{
            if(tmp->right) s2.push(tmp->right);
            if(tmp->left) s2.push(tmp->left);
        }
        s1.pop();
        
        if(s1.empty()){
            leftToRight = !leftToRight;
            swap(s1,s2);
            res.push_back(*container);
            container = new vector<int>();
        }
    }
    return res;
}

int numDistinct(string S, string T) {
    int n = (int)S.length();
    int m = (int)T.length();
    vector<vector<int>> tracker(m, vector<int>(n, 0));
    
    if(T[0] == S[0]) tracker[0][0] = 1;
    else tracker[0][0] = 0;
    
    for(int i = 1;i<n;i++)
    {
        if(T[0] == S[n])
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
            if(T[i] == S[j])
            {
                tracker[i][j] = tracker[i-1][j] + tracker[i][j-1];
            }
            else
            {
                tracker[i][j] = tracker[i][j-1];
            }
        }
    }
    return tracker[m-1][n-1];
}

void connect(TreeLinkNode *root) {
    TreeLinkNode* prev= NULL;
    TreeLinkNode* nextStart = root;
    if(!root) return;
    while(nextStart)
    {
        TreeLinkNode* tmp = nextStart;
        nextStart = NULL;
        prev= NULL;
        while(tmp)
        {
            if(!nextStart)
            {
                nextStart = tmp->left? tmp->left:tmp->right;
            }
            
            if(tmp->left)
            {
                if(prev) prev->next = tmp->left;
                prev = tmp->left;
            }
            
            if(tmp->right)
            {
                if(prev) prev->next = tmp->right;
                prev = tmp->right;
            }

            tmp = tmp->next;
        }
    }
}

int maxSum = INT_MIN;

int worker(TreeNode *root) {
    if(!root) return 0;
    
    int left = worker(root->left);
    int right = worker(root->right);
    
    int sum = root->val;
    if(left>0) sum+= left;
    if(right>0) sum+=right;
    maxSum = max(sum, maxSum);
    
    return max(left,right)>0 ? root->val + max(left,right): root->val;
}


int maxPathSum(TreeNode *root) {
    return max(worker(root), maxSum);
}

int longestConsecutive(vector<int> &num) {
    unordered_map<int,int> tracker;
    unordered_set<int>visited;
    
    for(int i = 0;i< num.size();i++)
    {
        tracker.insert(make_pair(num[i],i));
    }
    
    int maxlength = INT_MIN;

    for(int i = 0;i<num.size();i++)
    {
        int count;
        if(visited.find(i) == visited.end())
        {
            int tmp = num[i];
            count = 1;
            visited.insert(i);
            while( tracker.find(++tmp) != tracker.end())
            {
                count ++;
                visited.insert(tracker.find(tmp)->second);
            }
            tmp = num[i];
            while(tracker.find(--tmp) != tracker.end())
            {
                count ++;
                visited.insert(tracker.find(tmp)->second);
            }
        }
        maxlength = max(maxlength, count);
        
    }
    return maxlength;
}


void worker(vector<pair<int,int>>& tmp, int i, int j, vector<vector<char>>& board, unordered_set<pair<int,int>>& visited)
{
    if(i >=0 && i< board.size() &&
       j >=0 && j < board[0].size() &&
       board[i][j] == 'O')//&&
//       visited.find(make_pair(i,j)) == visited.end())
    {
        //tmp.push_back(make_pair(i,j));
        //visited.insert(make_pair(i,j));
        //worker(tmp, i-1, j, board, visited);
        //worker(tmp, i+1, j, board, visited);
        //worker(tmp, i, j-1, board, visited);
        //worker(tmp, i, j+1, board, visited);
    }
}
 
/*
void solve(vector<vector<char>> &board) {
    int m = (int)board.size();
    int n = (int)board[0].size();
    unordered_set<pair<int,int>> visited;
    vector<vector<pair<int,int>>> tracker;
    for(int i = 0;i< m;i++)
    {
        for(int j = 0;j<n;j++)
        {
            vector<pair<int,int>> tmp;
            if(visited.find(make_pair(i,j)) == visited.end() && board[i][j] == 'O')
            {
                worker(tmp, i, j, board, visited);
            }
            tracker.push_back(tmp);
        }
    }
    
    for(int i = 0;i< tracker.size();i++)
    {
        bool inside = true;
        for(auto it : tracker[i])
        {
            if(it.first == 0 ||
               it.first == m-1 ||
               it.second == n-1 ||
               it.second == 0)
            {
                inside = false;
                break;
            }
        }
        
        if(inside)
        {
            for(auto it:tracker[i])
            {
                board[it.first][it.second] = 'X';
            }
        }
    }
}
*/

int minCut(string s) {
    int n = (int)s.length();
    int f[n+1];
    bool tracker[n][n];
    fill_n(&tracker[0][0], n*n, false);
    for (int i = 0; i <= n; i++)
        f[i]=n-1-i;
    for(int i = n-1;i>=0;i++)
    {
        for(int j = i; j<n;j++)
        {
            if(s[i] == s[j] && (j-i <2 || tracker[i+1][j-1]))
            {
                tracker[i][j] = true;
                f[i] = min(f[i], f[j] +1);
            }
        }
    }
    return f[0];
}

int minCutI(string s) {
    const int n = (int)s.size();
    int f[n+1];
    bool p[n][n];
    fill_n(&p[0][0], n * n, false);
    //the worst case is cutting by each char
    for (int i = 0; i <= n; i++)
        f[i]=n-1-i;// 最后一个 f[n]=-1
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i; j < n; j++) {
            if (s[i] == s[j] && (j - i < 2 || p[i + 1][j - 1])) {
                p[i][j] = true;
                f[i] = min(f[i], f[j + 1] + 1);
            }
        }
    }
    return f[0];
}

int lengthOfLongestSubstring(string s) {
    if(s.length() == 0)return 0;
    int res = 1;
    int tracker[256];
    fill_n(&tracker[0], 256, 0);
    int left = 0;
    int right = 0;
    while(right < s.length())
    {
        tracker[s[right]]++;
        if(tracker[s[right]] == 1)
        {
            res = max(res, right-left+1);
        }
        else
        {
            while(tracker[s[left]] == 1)
            {
                tracker[s[left++]]--;
            }
            tracker[s[left]] --;
        }
    }
    return res;
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
    
    /*
    vector<int> height;
    height.push_back(2);
    height.push_back(3);
    
    int res = largestRectangleArea(height);
    cout <<res<<endl;
    */
    /*
    string s1 = "aa";
    string s2 = "ab";
    string s3 = "aaba";
    
    bool res = isInterleave(s1,s2,s3);
    cout <<res;
    */
    
    /*ListNode* l1 = new ListNode(1);
    ListNode* l2 = new ListNode(2);
    ListNode* l3 = new ListNode(3);
    l1->next = l2;
    l2->next = l3;
    
    ListNode* res = reverseKGroup(l1,1);
    cout<<res->val<<endl;
    */
    
    /*
    string haystack="a";
    string needle= "a";
    
    char* res = strStr(haystack.c_str(), needle.c_str());
    */
    /*
    int res = divide(2147483647, 1);
    cout<<res<<endl;
    */
    
    /*
    int res = longestValidParentheses("(()()");
    cout <<res;
    */
    //int input[] = {1,2,3};
    //vector<int> res = searchRange(input, 3, 2);
    /*
    vector<int> input(5,1);
    vector<vector<int>> res = combinationSum2(input, 3);
    */
    /*
    string s1 = "12";
    string s2 = "4";
    string res = multiply(s1, s2);
    cout<<res<<endl;
    */
    /*
    int A[] = {0};
    int res = jump(A, 1);
    cout<<res;
     */
    
    //vector<vector<string>> res = solveNQueens(8);
    
    /*vector<int> v1 = {1,2,3};
    vector<int> v2 = {4,5,6};
    vector<int> v3 = {7,8,9};
    
    vector<vector<int>> input = {v1,v2,v3};
    vector<int> res = spiralOrder(input);
    */
    
    vector<string> input = {"This", "is", "an", "example","of", "text","justification."};
    vector<string> res = fullJustify(input, 16);
    for(int i =0;i< res.size();i++)
    {
        cout<<res[i]<<endl;
    }
    
    //int res = sqrt(2);
    //cout <<res;
    
    //int res = minDistance("ab", "bc");
    //cout<<res<<endl;
    
    //vector<int> input = {2,1,2};
    //vector<vector<int>> res = subsetsWithDup(input);
    //cout<<res.size()<<endl;
    /*
    TreeNode* t1 = new TreeNode(1);
    TreeNode* t2 = new TreeNode(2);
    TreeNode* t3 = new TreeNode(3);
    TreeNode* t4 = new TreeNode(4);
    TreeNode* t5 = new TreeNode(5);
    
    t1->left = t2;
    t2->right = t3;
    t1->right = t4;
    t4->left = t5;
    
    vector<int> res = inorderTraversal(t1);
    cout<<res.size()<<endl;
    */
    //int res = numTrees(3);
    /*
    TreeLinkNode* t1 = new TreeLinkNode(1);
    TreeLinkNode* t2 = new TreeLinkNode(2);
    TreeLinkNode* t3 = new TreeLinkNode(3);
    TreeLinkNode* t4 = new TreeLinkNode(4);
    TreeLinkNode* t5 = new TreeLinkNode(5);
    
    t1->left = t2;
    t2->right = t3;
    t1->right = t4;
    t4->left = t5;
    
    connect(t1);
    */
    //TreeNode* t1 = new TreeNode(-3);
    
    //int res = maxPathSum(t1);
    //cout<< res<<endl;
    
    //vector<int> A = {-339,-524,-826,870,395,-173,805,690,798,992,96,-735};

    //int res = longestConsecutive(A);
    //cout<<res<<endl;
    
    
    // insert code here...
    //std::cout << "Hello, World!\n";
    return 0;
}