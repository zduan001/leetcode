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
    ListNode* fast = input;
    ListNode* slow = input;
    while(fast && fast->next){
        fast = fast->next->next;
        slow = slow->next;
    }
    
    ListNode* second = slow->next;
    slow->next = NULL;
    
    ListNode* newHead = NULL;
    ListNode* tmp = second;;
    while (second) {
        tmp = second;
        second = second->next;
        tmp->next = newHead;
        newHead = tmp;
    }
    
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

int main(int argc, const char * argv[])
{
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


    // insert code here...
    // std::cout << "Hello, World!\n";
    return 0;
}

