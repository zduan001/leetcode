//
//  main.cpp
//  practice
//
//  Created by Duan, David on 3/31/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <queue>
#include <ctime>
#include <assert.h>
using namespace std;

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

struct TreeNodeWithParent {
    int val;
    TreeNodeWithParent *left;
    TreeNodeWithParent *right;
    TreeNodeWithParent *parent;
    TreeNodeWithParent(int x) : val(x), left(NULL), right(NULL), parent(NULL) {}
};

struct TreeListNode {
    int val;
    TreeListNode *left;
    TreeListNode *right;
    vector<TreeListNode*> afterIt;
    TreeListNode* next;
    TreeListNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
};

struct Interval {
    int start;
    int end;
    Interval() : start(0), end(0) {}
    Interval(int s, int e) : start(s), end(e) {}
};

/*
 一个链表1->2->3->4->5转换成1->5->2->4->3
 http://www.mitbbs.com/article_t/JobHunting/32468029.html
 */
ListNode* interLeave(ListNode* input){
    // 1. split link list in to 2 parts
    ListNode* fast = input;
    ListNode* slow = input;
    while(fast && fast->next){
        fast = fast->next->next;
        slow = slow->next;
    }
    ListNode* second = slow->next;
    slow->next = NULL;
    
    
    // 2. reverse second part.
    ListNode* newHead = NULL;
    ListNode* tmp = second;;
    while (second) {
        tmp = second;
        second = second->next;
        tmp->next = newHead;
        newHead = tmp;
    }
    
    // 3. interleave 2 link list.
    ListNode* dummyHead = new ListNode(0);
    tmp = dummyHead;
    
    while(input && newHead){
        tmp->next = input;
        input = input->next;
        tmp = tmp->next;
        tmp->next = newHead;
        newHead = newHead->next;
        tmp = tmp->next;
    }
    
    if(input){
        tmp->next = input;
    }else{
        tmp->next = newHead;
    }
    
    return dummyHead->next;
}

/*
 FB
 1) 给一棵树，tree node里面增加定义一项 vector<TreeNode*> afterIt,也就是保存
 这一层中此节点之后（右侧）的所有节点指针。求如何实现这个vector。
 */
void populateRight(TreeListNode* n){
    TreeListNode* prev = NULL;
    TreeListNode* next = NULL;
    //TreeListNode* n = root;
    
    while(n){
        for(;n;n= n->next){
            if(!next) next = n->left? next = n->left:n->right;
            if(n->left){
                if(prev) prev->next = n->left;
                prev = n->left;
            }
            
            if(n->right){
                if(prev) prev->next = n->right;
                prev = n->right;
            }
        }
        n = next;
    }
}

void populateVector(TreeListNode* root){
    queue<TreeListNode*> q;
    queue<TreeListNode*> next;
    TreeListNode* tmp;
    
    q.push(root);
    while(!q.empty()){
        vector<TreeListNode*> v;
        while(!q.empty()){
            tmp = q.front();
            q.pop();
            tmp->afterIt = v;
            v.push_back(tmp);
            if(tmp->right) next.push(tmp->right);
            if(tmp->left) next.push(tmp->left);
        }
        swap(q, next);
    }
}


/*
int Read(int Size, char * buffer) using int Read4(char * buffer):
*/
int read4(char* buffer){
    int count = 0;
    for(count = 0;count < 4;count++)
    {
        if(*buffer!='\0')
        {
            count ++;
            buffer++;
        }
    }
    return count;
}

int Read(int size, char* buffer){
    int count = 0;
    while(count< size){
        int tmp = read4(buffer);
        
    }
    return 0;
}

struct elementWithParent{
    int x;
    int y;
    int val;
    elementWithParent(int i, int j, int v): x(i), y(j), val(v){}
};

int knapsack(vector<int> input, int capacity){
    int n = (int)input.size();
    elementWithParent* tracker[n][capacity];
    for(int i = 0;i< n;i++){
        for(int j = 0;j< capacity;j++){
            int x1, x2;
            x1 = j-input[i] >= 0 ? tracker[i-1][j-input[i]]->val: 0;
            x2 = tracker[i-1][j]->val;
            elementWithParent* tmp;
            if(x1+ input[i] > x2){
                tmp = new elementWithParent(i-1,j-input[i], x1+ input[i]);
            }else{
                tmp = new elementWithParent(i-1, j, x2);
            }
            tracker[i][j] = tmp;
        }
    }
    
    return tracker[n-1][capacity-1]->val;
}

/*
 Given a array of numbers, pick numbers not next to each other to achieve max sum :)
 */
int maxSum(int input[], int n)
{
    if(n<1) return 0;
    int tracker[n];
    tracker[0] = input[0];
    int sum;
    for(int i = 1;i< n;i++){
        sum = INT_MIN;
        for(int j = 0;j<i-1;j++){
            sum = max(sum, tracker[j] + input[i]);
        }
        tracker[i] = max(input[i], sum);
    }
    return tracker[n-1];
}

/*
 给一个二叉树,找到与给定节点距离为N的所有节点(没有parent link,有parent link),
 两个节点间隔着几条边,就是距离为几.
 */
queue<TreeNodeWithParent*> findNode(TreeNodeWithParent *node, int d){
    
    unordered_set<TreeNodeWithParent*> tracker;
    queue<TreeNodeWithParent*> q1;
    if(d == 0) return q1;
    queue<TreeNodeWithParent*> q2;
    
    q1.push(node);
    tracker.insert(node);
    int step = 0;
    while(!q1.empty()){
        TreeNodeWithParent* tmp = q1.front();
        if(tmp->parent && tracker.find(tmp->parent) == tracker.end()){
            q2.push(tmp->parent);
            tracker.insert(tmp->parent);
        }
        if(tmp->left && tracker.find(tmp->left) == tracker.end()){
            q2.push(tmp->left);
            tracker.insert(tmp->left);
        }
        if(tmp->right && tracker.find(tmp->right) == tracker.end()){
            q2.push(tmp->right);
            tracker.insert(tmp->right);
        }
        q1.pop();
        if(q1.empty()){
            step ++;
            if(step == d){
                return q2;
            }
            swap(q1, q2);
        }
        
    }
    
    // if we travelled all. return empty set.
    return q1;
    
}

void nextPermutation(vector<int> &num) {
    vector<int>::iterator first = num.begin();
    vector<int>::iterator last = num.end();
    
    if(first == last) return;
    vector<int>::iterator i = first;
    i++;
    if(i == last) return;
    i = last;
    i --;
    
    while(true){
        vector<int>::iterator ii = i;
        --i;
        if(*i < *ii){
            vector<int>::iterator j = last;
            while(!(*i < *--j)){}
            iter_swap(i,j);
            reverse(ii, last);
            return;
        }
        
        // the whole vector is in a reverse order.
        if(i == first)
        {
            reverse(first, last);
            return;
        }
    }
}

bool prevpermutation(vector<int>::iterator first, vector<int>::iterator last)
{
    if (first == last) return false;
    vector<int>::iterator i = last;
    if (first == --i) return false;
    
    while (1) {
        vector<int>::iterator i1 = i;
        --i;
        if (*i1 < *i) {
            vector<int>::iterator j = last;
            while (!(*--j < *i)){}
            std::iter_swap(i, j);
            std::reverse(i1, last);
            return true;
        }
        if (i == first) {
            std::reverse(first, last);
            return false;
        }
    }
}

int secArray[3600];
int secCount;

int lastsecond;
int lastmin;
int lasthour;
int lastday;

void processRequest(time_t t){
    tm *ltm = localtime(&t);
    int second = ltm->tm_sec;
    int minute = ltm->tm_min;
    int hour = ltm->tm_hour;
    
    if(second == lastsecond &&
       minute == lastmin &&
       hour == lasthour ){
        secCount++;
        return;
    }else{
        int secDiff = (hour-lasthour) * 3600 + (minute - lastmin) * 60 + (second - lastsecond);
        int lastIndex = lastmin * 60 + lastsecond;
        int currentIndex = minute * 60 + second;
        
        if(secDiff >=3600){
            //set whole array to 0;
        }
        else
        {
            secArray[lastIndex] = secCount;
            if(currentIndex > lastIndex){
                for(int i = currentIndex + 1;i< lastIndex;i++){
                    secArray[i] = 0;
                }
            }else{
                for(int i = currentIndex + 1;i< 3600;i++){
                    secArray[i] = 0;
                }
                for(int i = 0;i< lastIndex;i++){
                    secArray[i] = 0;
                }
            }
        }
        lasthour = hour;
        lastmin = minute;
        lastsecond = second;
        secCount ++;
        
    }
}

int findKthElement(int A[], int m, int B[], int n, int k)
{
    assert(A && B);
    if(m<=0) return B[k-1];
    if(n<=0) return A[k-1];
    
    if(k<=0) return min(A[0], B[0]);
    
    if(m/2 + n/2 +1 >=k)
    {
        // k must be in left half..
        if(A[m/2] > B[n/2]){
            // remove second half of A. A[m/2] > B[n/2] we know
            // the second half of A will not contain the kth element.
            return findKthElement(A, m/2, B, n, k);
        }else{
            return findKthElement(A, m, B, n/2, k);
        }
    }else{
        // kth element shoudl be in right half..
        if(A[m/2] >= B[n/2]){
            return findKthElement(A, m, B+n/2 +1, n - n/2 -1, k - n/2-1);
        }else{
            return findKthElement(A+m/2+1, m-m/2-1, B, n, k-m/2-1);
        }
    }
}

/*
 bool isInterleave(string s1, string s2, string s3) {
 if(s3.length() != s1.length() + s2.length()) return false;
 if(s1 == "") return s2 == s3;
 if(s2 == "") return s1 == s3;
 bool tracker[s1.length()+1][s2.length()+1];
 for(int i = 0; i< s1.length()+1;i++){
 for(int j = 0;j<s2.length()+1; j++){
 tracker[i][j] =false;
 }
 }
 tracker[0][0] = true;
 for(int i = 1;i< s1.length() +1; i++)
 {
 tracker[i][0] = s1.substr(s1.length() -i, i) == s3.substr(s3.length()-i, i);
 //cout<<i<<", "<<0<<": "<<tracker[i][0]<<endl;
 }
 for(int j = 1; j< s2.length()+1; j++)
 {
 tracker[0][j] = s2.substr(s2.length() - j, j) == s3.substr(s3.length()-j, j);
 //cout<<0<<", "<<j<<": "<<tracker[0][j]<<endl;
 }
 
 for(int i = 1;i< (int)s1.length()+1; i++){
 for(int j= 1;j< (int)s2.length()+1; j++){
 if(s1[s1.length() -i] != s3[s3.length() - i-j] && s2[s2.length()-j] != s3[s3.length()-i-j] ){
 tracker[i][j] = false;
 }else
 {
 if(s1[s1.length() -i] == s3[s3.length() - i-j])
 {
 tracker[i][j] = tracker[i][j] || tracker[i-1][j];
 }
 if(s2[s2.length()-j] == s3[s3.length()-i-j])
 {
 tracker[i][j] = tracker[i][j] || tracker[i][j-1];
 }
 }
 //cout<<i<<", "<<j<<": "<<tracker[i][j]<<endl;
 }
 }
 return tracker[s1.length()][s2.length()];
 }
 */


bool isInterleave(string s1, string s2, string s3) {
    if(s1.length() + s2.length() != s3.length()) return false;
    bool tracker[s1.length()+1][s2.length()+1];
    for(int i = 0;i<s1.length()+1;i++){
        for(int j = 0;j< s2.length()+1;j++){
            tracker[i][j] = false;
        }
    }
    
    tracker[0][0] = true;
    for(int i = 1;i<s1.length()+1;i++){
        if(s1[i-1] == s3[i-1]){
            tracker[i][0] = tracker[i-1][0];
            
        }
        cout<< i << "  " << "0" << "   " << tracker[i][0] <<endl;
    }
    for(int i = 1;i< s2.length()+1;i++){
        if(s2[i-1] == s3[i-1]){
            tracker[0][i] = tracker[0][i-1];
            
        }
        cout<< "0" << "  " << i << "   " << tracker[0][i] <<endl;
    }
    for(int i = 1;i< s1.length()+1;i++){
        for(int j = 1;j<s2.length()+1;j++){
            if(s1[i-1] == s3[i+j-1] && s2[j-1] == s3[i+j-1]){
                tracker[i][j] = tracker[i-1][j] || tracker[i][j-1];
            }else if(s1[i-1]== s3[i+j-1]){
                tracker[i][j] = tracker[i-1][j];
            }else if(s2[j-1] == s3[i+j-1]){
                tracker[i][j] = tracker[i][j-1];
            }
            cout<< i << "  " << j << "   " << tracker[i][j] <<endl;
        }
        
    }
    
    return tracker[s1.length()][s2.length()];
    
}

#define ONE_YEAR (365 * 24 * 3600)

struct stockPrice{
    time_t time;
    float price;
    stockPrice(time_t t, float p) : time(t), price(p){}
};

deque<stockPrice> minDeque;

void insertStockPrice(stockPrice s){
    time_t currentTime;
    time(&currentTime);
    while(!minDeque.empty() ||
          difftime(currentTime, (minDeque.back()).time) > ONE_YEAR ||
          (minDeque.back()).price > s.price)
    {
        minDeque.pop_back();
    }
    minDeque.push_back(s);
    
}

float fetchMin()
{
    time_t currentTime;
    time(&currentTime);
    while(difftime(currentTime, (minDeque.front()).time) > ONE_YEAR)
    {
        minDeque.pop_front();
    }
    if(minDeque.empty()){
        return (minDeque.front()).price;
    }
    return -1.0;
}


int main(int argc, const char * argv[])
{
    /**
    ListNode* l1 = new ListNode(1);
    ListNode* l2 = new ListNode(2);
    ListNode* l3 = new ListNode(3);
    ListNode* l4 = new ListNode(4);
    ListNode* l5 = new ListNode(5);
    
    l1->next = l2;
    l2->next = l3;
    l3->next = l4;
    l4->next = l5;
    
    ListNode* res = interLeave(l1);
    cout<<res<<endl;
*/
    //int A[] = {1,2,3,4,5};
    //int res = maxSum(A, 5);
    //cout<<res<<endl;
    
    /*
    char xs[] = "AAADKRRV";
    do
    {
        std::puts(xs);
    }
    while (next_permutation(xs, xs + sizeof(xs) - 1));
    //return 0;
    
    vector<int> input = {1,2,4,3};
    prevpermutation(input.begin(), input.end());
*/
    //string s1 = "ef";
    //string s2 = "gh";
    //string s3 = "ehgf";
    
    //bool res = isInterleave(s1, s2,s3);
    //cout<<res<<endl;

    // insert code here...
    // std::cout << "Hello, World!\n";
    printf("%d", 4&3);
    return 0;
}

