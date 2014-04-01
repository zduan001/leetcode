//
//  additional.cpp
//  LeetCodeXCode
//
//  Created by Duan, David on 3/31/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include "additional.h"

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

struct Interval {
    int start;
    int end;
    Interval() : start(0), end(0) {}
    Interval(int s, int e) : start(s), end(e) {}
};

/*
 一个链表1->2->3->4->5转换成1->5->2->4->3
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




















