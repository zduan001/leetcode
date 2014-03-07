//
//  main.cpp
//  LeetCode
//
//  Created by Duan, David on 3/4/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//
#include "container.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <set>
#include <stack>

using namespace std;

struct ListNode {
    int val;
    ListNode *next;
    ListNode(int x) : val(x), next(NULL) {}
};


ListNode *insertionSortList(ListNode *head) {
	ListNode *dummyHead = new ListNode(-1);
	ListNode *tmp;
	ListNode *tail;
    
	while(head){
		tmp = head;
		head = head->next;
		tail = dummyHead;
		while(tail->next != NULL && tail->next->val < tmp->val){
			tail = tail->next;
		}
		tmp->next = tail->next;
		tail->next = tmp;
	}
    
	return dummyHead->next;
}

struct RandomListNode {
    int label;
    RandomListNode *next, *random;
	RandomListNode(int x) : label(x), next(NULL), random(NULL) {}
};

RandomListNode *copyRandomList(RandomListNode *head) {
	unordered_map<RandomListNode *, RandomListNode*> tracker;
    
	RandomListNode *tmp = head;
	while(tmp){
		tracker.insert(make_pair(tmp, new RandomListNode(tmp->label)));
		tmp = tmp->next;
	}
    
	tmp = head;
	for(unordered_map<RandomListNode*, RandomListNode*>::iterator iterator = tracker.begin(); iterator != tracker.end(); iterator++){
		iterator->second->next = tracker[iterator->first->next];
		iterator->second->random = tracker[iterator->first->random];
	}
    
	return tracker[head];
}

RandomListNode *copyRandomListII(RandomListNode *head){
    if(!head) return NULL;
    
    RandomListNode *tmp = head;
    while (tmp) {
        RandomListNode *newNode = new RandomListNode(tmp->label);
        newNode->next = tmp->next;
        tmp->next = newNode;
        tmp = tmp->next->next;
    }
    
    tmp = head;
    while(tmp ){
        if(tmp->random){
            tmp->next->random = tmp->random->next;
        }
        tmp = tmp->next->next;
    }
    
    RandomListNode *res = head->next;
    tmp = head;
    RandomListNode *n = tmp->next;
    while(tmp && tmp->next){
        tmp->next = tmp->next->next;
        if(n->next){
            n->next = n->next->next;
        }
        tmp = tmp->next;
        n = n->next;
        
    }
    return res;
}

/**
* http://oj.leetcode.com/problems/word-break/
 */
bool wordBreak(string s, unordered_set<string> dict) {
    vector<bool> tracker;
    
    if(dict.find(s.substr(0,1)) != dict.end()){
        tracker.push_back(true);
    }
    else
    {
        tracker.push_back(false);
    }
    
    bool found= false;
    for(int i = 1;i<s.length();i++){
        found = false;
        for(int k = i-1;k>=-1;k--){
            if(dict.find(s.substr(k+1,i-k)) != dict.end() && (k == -1 || tracker[k])){
                found = true;
                break;
            }
        }
        tracker.push_back(found);
        
    }
    return tracker[s.length()-1];
    
}

/*
 http://oj.leetcode.com/problems/word-break-ii/
 */
vector<string> wordBreakII(string s, unordered_set<string> &dict) {
    vector<string> res;
    vector<bool> tracker;
    vector<vector<int>> length;
    vector<int> first;
    
    if(dict.find(s.substr(0,1)) != dict.end()){
        tracker.push_back(true);
        first.push_back(1);
    } else {
        tracker.push_back(false);
    }
    
    for(int i = 1;i< s.length(); i++){
        bool found = false;
        for(int k = i - 1; k>=-1;k++){
            if(dict.find(s.substr(k+1, i-k)) != dict.end() && (k == -1 || tracker[k])){
                found = true;
                vector<int> tmp;
                if(k >=0){
                    length[k].push_back(i-k);
                }else{
                    first.push_back(i-k);
                }
                tracker.push_back(found);
                length.push_back(tmp);
               
            }
        }
    }
    
    for(int index: first){
        string tmp;
        tmp += s.substr(0, index);
        //combineString(length, tmp, s, index);
        
    }
    
    
    return res;
}

void combinestirng(vector<vector<int>> length, int index, string s, vector<string> res){
    
}

struct classcomp {
    bool operator() (const char& lhs, const char& rhs) const
    {return lhs>rhs;}
};

//http://oj.leetcode.com/problems/candy/
int candy(vector<int> &ratings) {
    
    map<int,vector<int>> tracker;
    for(int i = 0;i< ratings.size(); i++){

        if(tracker.find(ratings[i]) != tracker.end()){
            (tracker.find(ratings[i])->second).push_back(i);

        }
        else
        {
            vector<int> tmp = {i};
            tracker.insert(make_pair(ratings[i], tmp));
        }
    }
    
    vector<int> *res = new vector<int>(ratings.size());
    
    for(map<int,vector<int>>::iterator it = tracker.begin(); it != tracker.end();it++){
        for(vector<int>::iterator vit = it->second.begin(); vit != it->second.end(); vit ++){
            if((*vit) -1 >= 0 && tracker[*vit] > tracker[(*vit)-1] && (*res)[*vit] <= (*res)[*vit-1]){
                (*res)[*vit] = (*res)[*vit-1]+1;
            }
            if((*vit)+1 < ratings.size() && tracker[*vit] > tracker[(*vit)+1] && (*res)[*vit] <= (*res)[*vit+1]){
                (*res)[*vit] = (*res)[*vit+1]+1;
            }
            
            
        }
    }
    int sum = 0;
    for(int i = 0;i<res->size(); i++){
        sum += (*res)[i];
    }
    sum += (*res).size();
    return sum;
    
}

//http://oj.leetcode.com/problems/candy/
int candyII(vector<int> &ratings){
    if(ratings.size() == 0){
        return 0;
    }
    int res = 0;
    int candy = 1;
    int seqLen = 0;
    int maxCntInSeq = candy;
    
    res += candy;
    seqLen++;
    for(int i = 1;i< ratings.size(); i ++){
        
        if(ratings[i] < ratings[i-1]){
            seqLen ++;
            if(maxCntInSeq == seqLen){
                seqLen++;
            }
            res += seqLen;
            candy = 1;
            
        }else{
            if(ratings[i] > ratings[i-1]){
                candy++;
            }else{
                candy = 1;
            }
            res += candy;
            seqLen = 0;
            maxCntInSeq = candy;
        }
    }
    
    return res;
}

//http://oj.leetcode.com/problems/reverse-words-in-a-string/
/*
 
 What constitutes a word?
 A sequence of non-space characters constitutes a word.
 Could the input string contain leading or trailing spaces?
 Yes. However, your reversed string should not contain leading or trailing spaces.
 How about multiple spaces between two words?
 Reduce them to a single space in the reversed string.

 */
bool isWhiteSpace(char c){
    return c == ' ' || c== '\n' || c == '\t';
}

void reverseWords(string &s){
    stack<string> st;
    int startindex = -1;
    if(s[0] != ' '){
        startindex = 0;
    }
    for(int i = 1;i<s.length();i++){
        if(startindex == -1){
            if(!isWhiteSpace(s[i])){
                startindex = i;
            }
        }else{
            if(isWhiteSpace(s[i])){
                st.push(s.substr(startindex, i-startindex));
                startindex = -1;
            }
        }
    }
    if(startindex != -1){
        st.push(s.substr(startindex, s.length()-startindex));
    }
    
    string res;
    while(!st.empty()){
        if(res.length() != 0){
            res += " ";
        }
        res += st.top();
        st.pop();
    }
    
    s = res;
    
}

void* memcopy(void *dest, const void *source, int lengh){
    char *to = (char*) dest;
    char *from = (char*) source;
    while(lengh){
        *to++ = *from++;
        lengh--;
    }
    return dest;
}

/*
 http://oj.leetcode.com/problems/clone-graph/
 */
struct UndirectedGraphNode {
    int label;
    vector<UndirectedGraphNode *> neighbors;
    UndirectedGraphNode(int x) : label(x) {};
};

void worker(unordered_map<UndirectedGraphNode*, UndirectedGraphNode*>* map, UndirectedGraphNode* node){
    if(map->find(node) == map->end()){
        UndirectedGraphNode* newNode = new UndirectedGraphNode(node->label);
        map->insert(make_pair(node, newNode));
        for(vector<UndirectedGraphNode*> :: iterator it = node->neighbors.begin() ; it != node->neighbors.end(); it ++){
            worker(map, *it);
            newNode->neighbors.push_back((map->find(*it))->second);
        }
    }
}

UndirectedGraphNode *cloneGraph(UndirectedGraphNode *node) {
    if(!node) return NULL;
    unordered_map<UndirectedGraphNode*, UndirectedGraphNode*>* tracker = new unordered_map<UndirectedGraphNode*, UndirectedGraphNode*>();
    worker(tracker, node);
    unordered_map<UndirectedGraphNode*, UndirectedGraphNode*>::iterator it = tracker->find(node);
    if(it != tracker->end()){
        return it->second;
    }
    else{
        return NULL;
    }
}

int main(int argc, const char * argv[])
{
    /*
	cout<<"start of the program"<<endl;
	RandomListNode *n1 = new RandomListNode(1);
	RandomListNode *n2 = new RandomListNode(2);
	RandomListNode *n3 = new RandomListNode(3);
	RandomListNode *n4 = new RandomListNode(4);
    
	n1->next = n2;
	n2->next = n3;
	n3->next = n4;
	n1->random = n3;
	n2->random = n2;
	n3->random = n2;
    n4->random = n1;
    
    RandomListNode *n5 = new RandomListNode(5);
    
	RandomListNode *res = copyRandomList(n1);
    res = copyRandomListII(n5);
	cout<<"end of the program"<<endl;
    
    string input = "dogs";
    unordered_set<string> dict = {"dog", "s", "gs"};
    bool breaker = wordBreak(input, dict);
     */
    /*
    set<int> myset = {2, 3, 4, 5};
    vector<int> vec;
    
    //int x = multiplies<int>(3, 4);
    transform(myset.begin(), myset.end(), back_inserter(vec), bind(multiplies<int>(), placeholders::_1,10));
    
    for(auto x: vec){
        cout<<x <<endl;
    }
    
    string res;
    string begin = "string to copy";
    cout<< res <<endl<<begin<<endl;
    memcopy(&res, &begin, 10);
    cout<< res <<endl<<begin<<endl;
    */
    //container *c = new container();
    //c->testMath();
    //vector<int> input = {2,3,2};
    //int res = candyII(input);
	//cout<<res<<endl;
    string s= "1 ";
    reverseWords(s);
    
    return 0;
}

