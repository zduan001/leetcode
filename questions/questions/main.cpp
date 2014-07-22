//
//  main.cpp
//  questions
//
//  Created by Duan, David on 6/10/14.
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

/*
 http://www.mitbbs.com/article_t/JobHunting/32466847.html
 */
int findCeiling(TreeNode* root, int target)
{
    int res = INT_MAX;
    while(root)
    {
        if(target == root->val)
        {
            root = root->right;
        }
        else if(target <root->val)
        {
            res = min(res, root->val);
            root = root->left;
        }
        else if(target > root->val)
        {
            root = root->right;
        }
    }
    return res;
}

bool worker(unordered_set<int> & visited, unordered_set<int>& tmp, int input[],int i, int n)
{
    visited.insert(i);
    tmp.insert(i);
    if(i>=0 && i<n)
    {
        i = input[i];
        if(tmp.find(i) != tmp.end())
        {
            return true;
        }
        if(visited.find(i) == visited.end())
        {
            bool x = worker(visited,tmp, input, i, n);
            if(x) return x;
        }
    }
    return false;
}

/*
 [Google] Given an array of integers where each element points to the index of 
 the next element how would you detect if there is a cycle in this array?
 */
bool isCircle(int input[], int n)
{
    unordered_set<int> visited;
    for(int i = 0;i< n;i++)
    {
        if(visited.find(i) == visited.end())
        {
            unordered_set<int> tmp;
            bool x = worker(visited, tmp, input, i, n);
            if(x)
                return true;
        }
    }
    return false;
}

/*
 http://www.mitbbs.com/article_t/JobHunting/32440405.html
 */
int longestSubStr(string s)
{
    if( s.length() ==0) return 0;
    unordered_set<char> tracker;
    int maxLength = 0;
    int left = 0;
    int right = 0;
    int p[26];
    fill_n(&p[0], 26, 0);
    p[s[0]-'a'] = 1;
    tracker.insert(s[0]);
    for(int i = 1;i<s.length();i++)
    {
        if(p[s[i]-'a'] >0 && tracker.size() <=2)
        {
            p[s[i]-'a']++;
            right = i;
            maxLength = max(maxLength, right-left +1);
        }
        else if(tracker.size() ==2)
        {
            while(tracker.size() >1 && left<=right)
            {
                p[s[left]-'a']--;
                if(p[s[left]-'a'] == 0)
                {
                    tracker.erase(s[left]);
                }
                left++;
            }
        }
        else if(p[s[i]-'a'] ==0)
        {
            tracker.insert(p[s[i]]);
            p[s[i]-'a'] ++;
            right = i;
            maxLength = max(maxLength, right-left + 1);
        }

    }
    return maxLength;
}

//    int p[] = {1,5,8,9,10,17,17,20,24,30};
int CutRod(int p[], int n)
{
    int tracker[n+1];
    fill_n(&tracker[0], n, 0);
    int tmp;
    for(int i = 1;i<= n;i++)
    {
        tmp = INT_MIN;
        for(int j = 1;j<=i;j++)
        {
            tmp = max(tmp, p[j-1] + tracker[i-j]);
        }
        tracker[i] = tmp;
    }
    
    return tracker[n];
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

int sqrts(int x) {
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

vector<int> inorderTraversal(TreeNode *root) {
    vector<int> res;
    if(!root) return res;
    bool backTrack = false;
    stack<TreeNode* > tracker;
    tracker.push(root);
    while(!tracker.empty())
    {
        TreeNode* tmp = tracker.top();
        if(tmp->left && !backTrack)
        {
            tracker.push(tmp->left);
        }
        else
        {
            res.push_back(tmp->val);
            tracker.pop();
            backTrack = true;
            
            if(tmp->right)
            {
                tracker.push(tmp->right);
                backTrack = false;
            }
        }
    }
    return res;
}

int main(int argc, const char * argv[])
{
    /*
    TreeNode* t8 = new TreeNode(8);
    TreeNode* t3 = new TreeNode(3);
    TreeNode* t12 = new TreeNode(12);
    TreeNode* t2 = new TreeNode(2);
    TreeNode* t6 = new TreeNode(6);
    TreeNode* t4 = new TreeNode(4);
    TreeNode* t10 = new TreeNode(10);
    TreeNode* t15 = new TreeNode(15);
    
    t8->left = t3;
    t8->right = t12;
    t3->left = t2;
    t3->right = t6;
    t6->left = t4;
    t12->left = t10;
    t12->right = t15;
    
    int res;
    res = findCeiling(t8, 13);
    cout<< "13=>"<<res<<endl;
    res = findCeiling(t8, 4);
    cout<<"4=>"<<res<<endl;
    res = findCeiling(t8, 8);
    cout<<"8=>"<<res<<endl;
    res = findCeiling(t8, 2);
    cout<<"2=>"<<res<<endl;
    */
    
    /*
    int input[] = {1,2,3,4,5,6};
    bool res = isCircle(input, 6);
    cout<<res<<endl;
    int input1[] = {1,2,0,4,5,3};
    res = isCircle(input1, 6);
    cout<<res<<endl;
    */
    
    //int res = longestSubStr("abbcccccce");
    //cout<<res<<endl;
    
    //int p[] = {1,5,8,9,10,17,17,20,24,30};
    //int res = CutRod(p, 4);
    //cout<<res<<endl;

    
    // insert code here...
//    std::cout << "Hello, World!\n";
    
    //vector<string> input = {"Listen","to","many,","speak","to","a","few."};
    //vector<string> input = {"a", "b", "c", "d", "e"};
    //vector<string> res = fullJustify(input, 3);
    //cout<<res.size()<<endl;
    //int res =sqrts(9);
    //cout<<res<<endl;
    
    vector<int> S = {1,2,2};
    vector<vector<int>> res = subsetsWithDup(S);
    
    return 0;
}

