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

struct UndirectedGraphNode {
     int label;
     vector<UndirectedGraphNode *> neighbors;
     UndirectedGraphNode(int x) : label(x) {};
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

void swap(vector<int>& num, int i, int j)
{
    int tmp = num[i];
    num[i] = num[j];
    num[j] = tmp;
}

void reverse(vector<int>&num, int i, int j)
{
    while(i<j)
        swap(num, i++,j--);
}

void nextPermutation(vector<int> &num) {
    int i = num.size()-1;
    
    while(i>0 && num[i-1] >= num[i])
    {
        i--;
    }
    
    if(i == 0)
    {
        reverse(num, 0, num.size()-1);
        return;
    }
    
    int j = num.size() -1;
    while(num[j] <= num[i-1] && j >=i)
    {
        j --;
    }
    swap(num, i-1, j);
    reverse(num, i, num.size()-1);
}

void swap(int A[], int i, int j)
{
    int tmp = A[i];
    A[i] = A[j];
    A[j] = tmp;
}

int firstMissingPositive(int A[], int n) {
    if (n == 0)
        return 1;
    int i = 0;
    while(i< n)
    {
        if((A[i] < 1 || A[i] > n) || A[i] == A[A[i]-1])
            i++;
        else
            swap(A, i, A[i] -1);
    }
    
    for(int i = 0;i<n;i++)
    {
        if(A[i] != i+1)
        {
            return i+1;
        }
    }
    return n+1;
}

/*
bool isMatchWildCard(const char *s, const char *p) {
    if(*p == '\0') return *s == '\0';
    
    if(*p == '*' && *s != '\0')
    {
        return isMatchWildCard(s, p+1) || isMatchWildCard(s+1, p);
    }
    else if((*p == '?' && *s != '\0') || *p == *s)
        return isMatchWildCard(s+1, p+1);
    else
        return false;
}
*/

bool isMatchWildCard(const char *s, const char *p) {
    const char* start = NULL;
    const char* ss = s;
    while(*s)
    {
        if(*p == *s || *p == '?'){p++;s++;}
        else if(*p == '*'){start = p++;ss =s;}
        else if(start){p = start+1;s = ++ss;}
    }
    while(*p == '*'){p++;}
    return *p == '\0';
}

bool search(int A[], int n, int target) {
    int left = 0;
    int right = n-1;
    while (left <=right)
    {
        int mid = left + (right-left)/2;
        if(A[mid] == target) return true;
        if(A[left] < A[mid])
        {
            if(target >=A[left] && target <= A[mid])
                right = mid-1;
            else
                left = mid+1;
        }
        else if(A[left] > A[mid])
        {
            if(target <= A[mid]  && target >= A[left])
                right = mid-1;
            else
                left = mid +1;
        }
        else if(A[left] == A[mid])
        {
            left++;
        }
    }
    return false;
}

int largestRectangleArea(vector<int> &height) {
    if(height.size() == 0) return 0;
    height.push_back(0);
    int res = 0;
    stack<int> tracker;
    int i = 0;
    while(i < height.size())
    {
        if(tracker.empty() || height[i] > height[tracker.top()])
        {
            tracker.push(i++);
        }
        else
        {
            int h = height[tracker.top()];
            tracker.pop();
            int index = tracker.empty()?-1:tracker.top();
            res = max(res, h * (i- index-1));
        }
    }
    return res;
}

int maximalRectangle(vector<vector<char> > &matrix) {
    if(matrix.size() == 0 || matrix[0].size() == 0) return 0;
    vector<int> height(matrix[0].size()+1, 0);
    int res = 0;
    for(int i = 0;i< matrix.size();i++)
    {
        for(int j = 0;j<matrix[0].size();j++)
        {
            if(matrix[i][j] == '1') height[j]++;
            else if(matrix[i][j] == '0') height[j] = 0;
        }
        res = max(res, largestRectangleArea(height));
    }
    return res;
}


bool isScramble(string::iterator first1, string::iterator last1, string::iterator first2)
{
    auto length =distance(first1, last1);
    auto last2 = next(first2, length);
    
    if(length ==1) return *first1 == *first2;
    
    for(int i = 1;i<length;i++)
    {
        if((isScramble(first1, first1 + i , first2) && isScramble(first1+i , last1, first2+i))
           ||(isScramble(first1, first1+i, last2-i) && isScramble(first1+i, last1, first2)))
           return true;
    }
    return false;
    
}

bool isScramble(string s1, string s2)
{
    if(s1.length() != s2.length()) return false;
    return isScramble(s1.begin(), s1.end(), s2.begin());
}

bool isScrambleDp(string s1, string s2)
{
    int n = (int)s1.length();
    if(s1.length() != s2.length()) return false;
    
    bool f[n+1][n][n];
    fill_n(&f[0][0][0], (n+1)*n*n, false);
    
    for(int i = 0;i< n;i++)
    {
        for(int j = 0;j<n;j++)
        {
            if(s1[i] == s2[j])
            {
                f[1][i][j] = true;
            }
        }
    }
    
    for(int m  = 2;m<=n;m++)
    {
        for(int i = 0;i< n;i++)
        {
            for(int j = 0;j<n;j++)
            {
                for(int k = 1;k<m;k++)
                {
                    if((f[k][i][j] && f[n-k][i+k][j+k]) ||
                       (f[k][i][j+n-k] && f[n-k][i+k][j]))
                    {
                        f[n][i][j] = true;
                        break;
                    }
                }
            }
        }
    }
    return f[n][0][0];
}

bool isInterleave(string s1, string s2, string s3) {
    int n = (int)s1.length();
    int m = (int)s2.length();
    int k = (int)s3.length();
    if(n+m != k ) return false;
    bool tracker[n+1][m+1];
    fill_n(&tracker[0][0], (n+1)*(m+1), true);
    for(int i = 1;i<=n;i++)
    {
        tracker[i][0] = s1[i-1] == s3[i-1] && tracker[i-1][0];
            }
    for(int j = 1;j<=m;j++)
    {
        tracker[0][j] = s2[j-1] == s3[j-1] && tracker[0][j-1];
            }
    
    for(int i = 1;i<= n;i++)
    {
        for(int j = 1;j<=m;j++)
        {
            tracker[i][j] = (s1[i-1] == s3[i+j-1] && tracker[i-1][j]) ||
            (s2[j-1] == s3[i+j-1] && tracker[i][j-1]);
        }
    }
    return tracker[n][m];
}

string longestPalindrome(string s) {
    int n = (int)s.length();
    if( n==0) return "";
    int mL = -1;
    int mR = -1;
    int maxLength = -1;
    
    for(int i = 0;i<n;i++)
    {
        int left = i;
        int right = i;
        while(left >=0 && right <n && s[left] == s[right]){
            if(right-left > mR-mL){
                mL=left;
                mR = right;
            }
            left--;right++;
        }
        
        
        left = i;
        right = i+1;
        while(left>=0 && right <n && s[left] == s[right])
        {
            if(right-left > mR-mL){
                mL= left;
                mR = right;
            }
            left --;right++;
        }
    }
    return s.substr(mL, mR-mL+1);
}

int maxProfit(vector<int> &prices) {
    if(prices.size() < 2) return 0;
    int res = 0;
    int low = prices[0];
    for(int i = 1;i<prices.size();i++)
    {
        if(prices[i] < low)
        {
            low = prices[i];
        }
        res = max(res, prices[i] - low);
    }
    return res;
}

int maxProfitII(vector<int> &prices) {
    if(prices.size() <2) return 0;
    int res = 0;
    for(int i = 1;i<prices.size();i++)
    {
        if(prices[i] > prices[i-1])
        {
            res+=prices[i]-prices[i-1];
        }
    }
    return res;
}

int maxProfitIII(vector<int> &prices) {
    int n = (int)prices.size();
    if(n == 0) return 0;
    vector<int> f(n,0);
    vector<int> g(n,0);
    int low = prices[0];
    for(int i = 1;i<n;i++)
    {
        if(prices[i] < low)
        {
            low = prices[i];
        }
        f[i] = max(prices[i] - low, f[i-1]);
    }
    int high = prices[n-1];
    for(int i = n-2;i>=0;i--)
    {
        if(prices[i] > high)
        {
            high = prices[i];
        }
        g[i] = max(high-prices[i], g[i+1]);
    }
    int res = 0;
    for(int i = 0;i<n;i++)
    {
        res = max(res, f[i] + g[i]);
    }
    return res;
}

int maxSubArray(int A[], int n) {
    int res = 0;
    for(int i = 0;i< n;i++)
    {
        int f = max(res+A[i], A[i]);
        res = max(f, res);
    }
    return res;
}

int minimumTotal(vector<vector<int> > &triangle) {
    int n = (int)triangle.size();
    if(n == 0)return 0;
    for(int i = n-2;i>=0;i--)
    {
        for(int j = 0;j< (int)triangle[i].size();j++)
        {
            triangle[i][j] = triangle[i][j] + min(triangle[i+1][j] , triangle[i+1][j+1]);
        }
    }
    return triangle[0][0];
}

int maxArea(vector<int> &height) {
    int n = (int)height.size();
    if(n<2) return 0;
    int left = 0;
    int right = n-1;
    int res = 0;
    while(left< right)
    {
        res = max(res, (right-left)* min(height[left], height[right]));
        if(height[left] < height[right]) left++;
        else right--;
    }
    return res;
    
}

int lengthOfLongestSubstring(string s) {
    int start = 0;
    int end = 0;
    int res = 0;
    unordered_map<char, bool> t;
    while(start<s.size() && end<s.size())
    {
        if(!t[s[end]])
        {
            t[s[end]] = true;
            res = max(res, end-start+1);
            end++;
        }
        else
        {
            t[s[start]] = false;
            start++;
        }
    }
    return res;
}

bool canJump(int A[], int n) {
    if(n ==0) return true;
    int reach = 0;
    for(int i = 0;i<= reach;i++)
    {
        reach = max(reach, i+A[i]);
        if(reach >=n-1) return true;
    }
    return false;
}

int jump(int A[], int n) {
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

int minCut(string s) {
    int n = (int)s.size();
    int f[n+1];
    bool p[n][n];
    for(int i = 0;i<= n;i++)
    {
        f[i] = n-i -1;
    }
    
    fill_n(&p[0][0], n*n, false);
    
    for(int i = n-1;i>=0;i--)
    {
        for(int j =i;j<n;j++)
        {
            if(s[i] == s[j] && (j-i < 2 || p[i+1][j-1]))
            {
                p[i][j] = true;
                f[i] = min(f[i], f[j] + 1);
            }
        }
    }
    return f[0];
}

bool isInterleaveII(string s1, string s2, string s3) {
    int m = (int)s1.size();
    int n = (int)s2.size();
    int k = (int)s3.size();
    if(m+n != k ) return false;
    bool f[m+1][n+1];
    fill_n(&f[0][0], (m+1)*(n+1), false);
    f[0][0] = true;
    for(int i = 1;i<=m;i++)
    {
        if(s1[i-1] == s3[i-1] && f[i-1][0])
            f[i][0] = true;
    }
    
    for(int j = 1;j<=n;j++)
    {
        if(s2[j-1] == s3[j-1] && f[0][j-1])
            f[0][j] = true;
    }
    
    for(int i = 1;i<=m;i++)
    {
        for(int j = 1;j<=n;j++)
        {
            f[i][j] = (s1[i-1] == s3[i+j-1] && f[i-1][j])||
            (s2[j-1] == s3[i+j-1] && f[i][j-1]);
        }
    }
    return f[m][n];
}

int minPathSum(vector<vector<int>> &grid) {
    int m = (int)grid.size();
    int n = (int)grid[0].size();
    vector<vector<int>>f(m, vector<int> (n,0));
    
    for(int i = 0;i<m;i++)
    {
        for(int j = 0;j<n;j++)
        {
            f[i][j] = min( i-1>=0? f[i-1][j]: INT_MAX, j-1>=0? f[i][j-1]: INT_MAX) + grid[i][j];
        }
    }
    
    return f[m-1][n-1];
}

int minDistance(string word1, string word2) {
    int m = (int)word1.size();
    int n = (int)word2.size();
    if(m == 0) return n;
    if(n == 0) return m;
    
    int p[m+1][n+1];
    
    for(int i = 0;i<= m;i++)
    {
        p[i][0] = i;
    }
    for(int j = 0;j<=n;j++)
    {
        p[0][j] = j;
    }
    
    for(int i = 1;i<=m;i++)
    {
        for(int j = 1;j<=n;j++)
        {
            int tmp = word1[i-1] == word2[j-1]? p[i-1][j-1]: p[i-1][j-1] + 1;
            p[i][j] = min(min(p[i-1][j]+1, p[i][j-1]+1), tmp);
        }
    }
    return p[m][n];
}

int numDecodings(string s) {
    if(s.size() == 0) return 0;
    int prev = 0;
    int cur = 1;
    
    for(int i = 0;i<s.size();i++)
    {
        if(s[i] == '0') cur = 0;
        if(i<1 || !(s[i-1] == '1' || (s[i-1] == '2' && s[i] <= '6')))
        {
            prev = 0;
        }
        int tmp = cur;
        cur = cur+prev;
        prev = tmp;
    }
    return cur;
}

int numDistinct(string S, string T) {
    int n = (int)T.size();
    int m = (int)S.size();
    if(m <n) return false;
    
    int f[n+1];
    fill_n(f, n+1, 0);
    f[0] = 1;
    for(int i = 0;i< S.size();i++)
    {
        for(int j = n-1;j>=0;j--)
        {
            f[j+1] = f[j+1] + (S[i] == T[j]? f[j]: 0);
        }
    }
    return f[n];
}

bool wordBreak(string s, unordered_set<string> &dict)
{
    int n = (int)s.size();
    int f[n];
    fill_n(f, n+1, false);
    f[0] = true;
    for(int i = 0;i<n;i++)
    {
        for(int j = i;j>=0;j--)
        {
            if(f[j] && dict.find(s.substr(j, i-j+1)) != dict.end())
            {
                f[i] = true;
                break;
            }
        }
    }
    return f[n-1];
}

UndirectedGraphNode* worker(unordered_map<UndirectedGraphNode*, UndirectedGraphNode*>& tracker,UndirectedGraphNode* node)
{
    if(tracker.find(node) == tracker.end())
    {
        UndirectedGraphNode* tmp = new UndirectedGraphNode(node->label);
        tracker[node] = tmp;
        for(auto neighbor : node->neighbors)
        {
            (tmp->neighbors).push_back(worker(tracker, neighbor));
        }
        return tmp;
    }
    else
    {
        return tracker[node];
    }
}


UndirectedGraphNode *cloneGraph(UndirectedGraphNode *node) {
    if(!node) return NULL;
    unordered_map<UndirectedGraphNode*, UndirectedGraphNode*> tracker;
    return worker(tracker, node);
}


void gen_path(const string &s, const vector<vector<bool> > &prev,
    int cur, vector<string> &path, vector<string> &result) {
    if (cur == 0)
    {
        string tmp;
        for (auto iter = path.crbegin(); iter != path.crend(); ++iter)
            tmp += *iter + " ";
        tmp.erase(tmp.end() - 1);
        result.push_back(tmp);
    }
    
    for (int i = 0; i < s.size(); ++i) {
        if (prev[cur][i])
        {
            path.push_back(s.substr(i, cur - i));
            gen_path(s, prev, i, path, result);
            path.pop_back();
        }
    }
}

vector<string> wordBreakII(string s, unordered_set<string> &dict) {
    // 长度为 n 的字符串有 n+1 个隔板
    vector<bool> f(s.length() + 1, false);
    // prev[i][j] 为 true,表示 s[j, i) 是一个合法单词,可以从 j 处切开 // 第一行未用
    vector<vector<bool> > prev(s.length() + 1, vector<bool>(s.length()));
    f[0] = true;
    for (int i = 1; i <= s.length(); ++i)
    {
        for (int j = i - 1; j >= 0; --j)
        {
            if (f[j] && dict.find(s.substr(j, i - j)) != dict.end()) {
                f[i] = true;
                prev[i][j] = true;
            }
        }
    }
    vector<string> result;
    vector<string> path;
    gen_path(s, prev, (int)s.length(), path, result);
    return result;
}

int sqrt(int x)
{
    if (x<0) return -1;

    double prev = 0;
    double current = x;
    while(abs(prev - current) > 0.001)
    {
        prev = current;
        current = prev - (prev * prev - x)/(2*prev);
    }
    return (int) current;
}

double pow(double x, int n) {
    bool negative = false;
    if(n<0)
    {
        negative = true;
        n = -n;
    }
    if(n == 0) return 1;
    if(n == 1) return negative ? 1/x: x;
    double res = pow(x, n/2);
    res *=res;
    if(n & 1)
    {
        res = res*x;
    }
    
    return negative? 1/res: res;
}

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

void worker(string s, vector<vector<string>>& res, int start, vector<string>& path)
{
    if(start == s.length())
    {
        res.push_back(path);
        return;
    }
    else
    {
        for(int i = start; i<s.length();i++)
        {
            if(isPalindrome(s, start, i))
            {
                path.push_back(s.substr(start, i-start+1));
                worker(s, res, i+1, path);
                path.pop_back();
            }
        }
    }
}

vector<vector<string>> partition(string s) {
    vector<vector<string>> res;
    vector<string> path ;
    if(s.length() == 0) return res;
    worker(s, res, 0, path);
    return res;
}

int uniquePathsWithObstacles(vector<vector<int> > &obstacleGrid) {
    int m = obstacleGrid.size();
    int n = obstacleGrid[0].size();
    
    if(m == 0 && n == 0) return 0;
    
    vector<vector<int>> f(m, vector<int>(n, 0));
    if(obstacleGrid[0][0] == 1)
    {
        f[0][0] = 0;
    }
    else
    {
        f[0][0] =1;
    }
    
    for(int i = 1;i<m;i++)
    {
        if(obstacleGrid[i][0] == 0 && f[i-1][0] != 0)
        {
            f[i][0] = 1;
        }
        else
        {
            f[i][0] = 0;
        }
    }
    
    for(int j = 1;j<n;j++)
    {
        if(obstacleGrid[0][j] ==0 && f[0][j-1] != 0)
        {
            f[0][j] = 1;
        }
        else
        {
            f[0][j] = 0;
        }
    }
    
    for(int i = 1;i<m;i++)
    {
        for(int j = 1;j<n;j++)
        {
            if(obstacleGrid[i][j] == 1)
            {
                f[i][j] =0;
            }
            else
            {
                f[i][j] = f[i-1][j] + f[i][j-1];
            }
        }
    }
    
    return f[m-1][n-1];
}

int main(int argc, const char * argv[])
{
/*    int A[] = {1000};
    int B[] = {1001};
    double res = findMedianSortedArrays(A, 1, B, 1);
    cout<<res<<endl;
 */
    /*
    cout<<isMatch("aaa", "a*a")<<endl;
    cout<<isMatch("aa","a")<<endl;
    cout<<isMatch("aa","aa")<<endl;
    cout<<isMatch("aaa","aa")<<endl;
    cout<<isMatch("aa", "a*")<<endl;
    cout<<isMatch("aa", ".*")<<endl;
    cout<<isMatch("ab", ".*")<<endl;
    cout<<isMatch("aab", "c*a*b")<<endl;
*/
/*    vector<int> num = {1,5,1};
    nextPermutation(num);
    num = {5,1,1};
    nextPermutation(num);
 */
    
    /*int A[] = {1,2,0};
    int res = firstMissingPositive(A, 3);
    cout<<res<<endl;
    int B[] = {3,4,-1,1};
    res = firstMissingPositive(B, 4);
    cout<<res<<endl;
    int C[] = {2,3,1,4};
    res = firstMissingPositive(C, 4);,
    cout<<res<<endl;
    int D[] = {1,1};
    res = firstMissingPositive(D, 2);
    cout<<res<<endl;
     */
    
    /*bool res = false;
    res = isMatchWildCard("aa","a");
    cout<<res<<endl;
    res = isMatchWildCard("aa","aa");
    cout<<res<<endl;
    res = isMatchWildCard("aaa","aa");
    cout<<res<<endl;
    res = isMatchWildCard("aa", "*");
    cout<<res<<endl;
    res = isMatchWildCard("aa", "a*");
    cout<<res<<endl;
    res = isMatchWildCard("ab", "?*");
    cout<<res<<endl;
    res = isMatchWildCard("aab", "c*a*b");
    cout<<res<<endl;
    res = isMatchWildCard("", "****");
    cout<<res<<endl;
     */
    
    /*int A[] = {5,1,3};
    bool res = search(A,3,3);
    cout<<res<<endl;
    int B[] = {3,1,1};
    res = search(B, 3,3);
    cout<<res<<endl;
    */
    /*vector<int> height = {2,1,5,6,2,3};
    int res = largestRectangleArea(height);
    cout<<res<<endl;
    */
    
    /*
    string s1 = "aabcc";
    string s2 = "dbbca";
    
    bool res = isInterleave(s1, s2, "aadbbbaccc");
    cout<<res<<endl;
    */
    /*
    string res = longestPalindrome("ccc");
    cout<<res<<endl;
     */
    /*
    vector<int> input = {1,2,4};
    int res = maxProfitIII(input);
    cout<<res<<endl;
     */
    /*
    int res = lengthOfLongestSubstring("");
    cout<<res<<endl;
     */
    /*
    int A[] = {1,2,3};
    int res = jump(A, 3);
    cout<<res<<endl;
     */
    /*
     int res = numDecodings("10");
    cout<<res<<endl;
     */
    
    /*int res = numDistinct("d", "d");
    cout<<res<<endl;
     */
    /*unordered_set<string> dict = {"a"};
    bool res = wordBreak("a", dict);
    cout<<res<<endl;
    */
    /*unordered_set<string> dict;
    dict.insert("cat");
    dict.insert("cats");
    dict.insert("and");
    dict.insert("sand");
    dict.insert("dog");
    
    
    vector<string> res = wordBreakII("catsanddog", dict);
     */
    //double res = pow(34.00515, -3);
    vector<int> r1 = {0,0};
    vector<int> r2 = {0,0};
    vector<vector<int>> input;
    input.push_back(r1);
    input.push_back(r2);
    int res = uniquePathsWithObstacles(input);
    cout<<res<<endl;
    return 0;
}

