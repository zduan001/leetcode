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

int main(int argc, const char * argv[]) {
    // insert code here...
    vector<int> input = {4, 4, 1, 1, 8, 8, 1, 1};
    sort(input);
    for(int i = 0;i< input.size();i++)
    {
        cout<<input[i]<<" ";
    }
    return 0;
}
