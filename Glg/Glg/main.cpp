//
//  main.cpp
//  Glg
//
//  Created by Duan, David on 11/8/14.
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

/*
 http://www.mitbbs.com/article_t/JobHunting/32748027.html
 有个followup的问题就是：因为我不想过多移动这些变量，所以怎么才能设计一个算法
 所需要移动的object最少。
 比如如果变量的size一次是4, 4, 1, 1, 8, 8, 1, 1最好的排法是4, 4, 8, 8, 1, 1,
 1, 1.而不是8 8 4 4 1 1 1 1因为前一种所需要移动的cost最小。这个没想出来了。。
 应该用divide and conquer？
 
 [Quick sort to bring 4, 8 byte element to the front]
 */
void swap(vector<int>& input, int i, int j)
{
    int tmp = input[i];
    input[i] = input[j];
    input[j] = tmp;
}

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
    int index = 0;
    int j = 0;
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

char* readRandomStream(int n)
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
    return res;
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
 find longest substring which contains n distinct characters.
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
 
 Write a function which, given two integers (a numerator and a denominator), prints the decimal representation of the rational number "numerator/denominator".
 Since all rational numbers end with a repeating section, print the repeating section of digits inside parentheses; the decimal printout will be/must be
 
 Example:
 1 , 3 = 0.(3)
 2 , 4 = 0.5(0)
 22, 7 = 3.(142857)
 */
string divide(int a, int b)
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
int* merge(int a[], int n, int b[], int m)
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
    int* res = merge(left, k, right, n-k);
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
void swap(int a[], int i, int j)
{
    int tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
}

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
void worker(vector<int> input, int level, vector<int> tmp, vector<vector<int>>& res )
{
    if(level == input.size())
    {
        res.push_back(tmp);
        return;
    }
    else
    {
        worker(input, level+1, tmp, res);
        tmp.push_back(input[level]);
        worker(input, level+1, tmp, res);
        tmp.pop_back();
    }
}

vector<vector<int>> allMulti(vector<int> input)
{
    vector<vector<int>>res;
    if(input.size() ==0) return res;
    sort(input.begin(), input.end());
    vector<int> tmp;
    worker(input, 0, tmp, res);
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


int main(int argc, const char * argv[]) {
    /*
    vector<int> input = {4, 4, 1, 1, 8, 8, 1, 1};
    sort(input);
    for(int i = 0;i< input.size();i++)
    {
        cout<<input[i]<<" ";
    }
     */
    /*vector<string> res = findAllRotate(5);
    for(auto item: res)
    {
        cout<<item<<endl;
    }
     */
    //string res = findSub("aabbccccccd", 3);
    //cout<<res<<endl;
    
    /*string res = divide(22, 4);
    cout<<res<<endl;
     */
    //int a[] = {4,3,2,1};
    //int res = countInversion(a, 4);
    //cout<<res<<endl;
    
    /*int a[] = { 3,1,2,4,5};
    int b[] = { 2,4,5,1,3};
    match(a, b, 5);
    for(int i = 0;i<5;i++)
    {
        cout<< a[i] << "  "<<b[i] << endl;
    }
     */
    /*
    push(4);
    push(3);
    push(2);
    push(1);
    cout<<start->key<<endl;
    cout<<findMid()<<endl;
    pop();
    cout<<start->key<<endl;
    cout<<findMid()<<endl;
     */
    vector<int> input = {1,2,3};
    vector<vector<int>> res = allMulti(input);
    for(auto item: res)
    {
        cout<<item.size()<<endl;
    }
    
    
    return 0;
}
