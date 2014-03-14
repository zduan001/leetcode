//
//  main.cpp
//  LeetCode
//
//  Created by Duan, David on 3/4/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <set>
#include <stack>
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

/*
 http://oj.leetcode.com/problems/palindrome-partitioning/
 */
bool isPlindrome(string s, int startindex, int lastIndex){
    //single character is palindrome.
    if(startindex == lastIndex){
        return true;
    }
    int leftedge;
    int rightedge;
    if((lastIndex- startindex+1)%2 == 0){
        leftedge = startindex + (lastIndex - startindex)/2;
        rightedge = leftedge +1;
        while(leftedge >= startindex && rightedge <=lastIndex){
            if(s[leftedge] != s[rightedge]){
                return false;
            }
            leftedge --;
            rightedge ++;
        }
    }else{
        leftedge = startindex + (lastIndex - startindex) / 2 -1;
        rightedge = startindex + (lastIndex - startindex) / 2 +1;
        while(leftedge >= startindex && rightedge <=lastIndex){
            if(s[leftedge] != s[rightedge]){
                return false;
            }
            leftedge --;
            rightedge ++;
        }
    }
    return true;
}

void workfunction(vector<vector<string>>* res,vector<string>* inprogress, int startIndex, int lastIndex, string s){
    
    if(isPlindrome(s, startIndex, lastIndex)){
        if(lastIndex == s.length()-1){
            inprogress->push_back(s.substr(startIndex, lastIndex-startIndex +1));
            vector<string>* tmp = new vector<string>(*inprogress);
            //copy(inprogress->begin(), inprogress->end(), tmp->begin());
            res->push_back(*tmp);
            inprogress->pop_back();
        }else{
            
            for(int i = 0;i < s.length() - lastIndex; i++){
                inprogress->push_back(s.substr(startIndex, lastIndex-startIndex +1));
                workfunction(res, inprogress, lastIndex + 1, lastIndex +1 + i, s);
                inprogress->pop_back();
            }
        }
    }
}

vector<vector<string>> partition(string s) {
    vector<vector<string>>* res = new vector<vector<string>>();
    vector<string>* inprogress = new vector<string>();
    for(int i = 0; i < s.length(); i++){
        workfunction(res, inprogress, 0, i, s);
    }
    return *res;
}

/*
 http://oj.leetcode.com/problems/palindrome-partitioning-ii/
 */
int minPalindromeCut(string s){
    if((int)    s.length() < 2 || isPlindrome(s, 0, (int)s.length()-1)) {
        return 0;
    }
    vector<int> tracker;
    tracker.push_back(0);

    for(int i = 1;i< s.length();i++){
        tracker.push_back(INT_MAX);
        if(isPlindrome(s,0,i)){
            tracker[i] = 0;
        }
        for(int j = i;j>0;j--){
            if(tracker[j-1] != INT_MAX && isPlindrome(s,j, i)){
                tracker[i] = tracker[i] < tracker[j-1]+1? tracker[i]: tracker[j-1] +1;
                if(tracker[i] == 1) break;
            }
        }
    }
    if(tracker[s.length()-1] == INT_MAX){
        return (int)s.length()-1;
    }
    return tracker[s.length()-1];
}

/*
 http://oj.leetcode.com/problems/surrounded-regions/
 */
void findArea(vector<vector<char>>& board, int k, int l, int m, int n, char old, char newchar){

    stack<pair<int,int>> tracker;
    
    if(board[k][l] == old){
        tracker.push(make_pair(k, l));
    }
    
    while(!tracker.empty()){
        pair<int,int> tmp = tracker.top();
        tracker.pop();
        board[tmp.first][tmp.second] = newchar;
        k = tmp.first;
        l = tmp.second;
        
        if(k>0 && board[k-1][l] == old){
            tracker.push(make_pair( k-1, l));
        }
        if(k+1 <m && board[k+1][l] == old){
            tracker.push(make_pair( k+1, l));
        }
        if(l>0 && board[k][l-1] == old){
            tracker.push(make_pair( k, l-1));
        }
        if(l+1<n && board[k][l+1] == old){
            tracker.push(make_pair( k, l+1));
        }
    }
}

void solve(vector<vector<char>> &board) {
    if(&board == NULL) return;
    
    const int m = (int)board.size();
    if(m == 0) return;
    const int n = (int)board[0].size();
    if(n == 0) return;
    
    for(int i = 0;i< m;i++){
        for(int j = 0;j < n;j ++){
            if(board[i][j] == 'O'){
                //board[i][j] = 'A';
                findArea(board, i, j,m,n, 'O', 'A');
            }
        }
    }
    
    for(int j = 0; j<n;j++){
        if(board[0][j] == 'A'){
            findArea(board, 0, j, m, n, 'A', 'O');
        }
        if(board[m-1][j] == 'A'){
            findArea(board, m-1, j, m, n, 'A', 'O');
        }
    }
    
    for(int i = 0;i<m;i++){
        if(board[i][0] == 'A'){
            findArea(board, i, 0, m, n, 'A', 'O');
        }
        if(board[i][n-1] == 'A'){
            findArea(board, i, n-1, m, n, 'A', 'O');
        }
    }
    
    for(int i = 0;i< m;i++){
        for(int j = 0;j < n;j ++){
            if(board[i][j] == 'A'){
                board[i][j] = 'X';
                //findArea(board, i, j,m,n, 'A', 'X');
                
            }
        }
    }
}

/*
 http://oj.leetcode.com/problems/sum-root-to-leaf-numbers/
 */
int sumNumbers(TreeNode *root) {
    int res = 0;
    int tmp = 0;
    
    if(!root) {return 0;}
    
    stack<TreeNode*> tracker;
    TreeNode *cur;
    TreeNode *prev = NULL;
    tracker.push(root);
    tmp = root->val;
    while(!tracker.empty()){
        cur = tracker.top();
        if(cur->left && cur->left != prev &&
           (!prev || prev != cur->right)){
            tracker.push(cur->left);
            tmp = tmp * 10 + cur->left->val;
            prev = NULL;
        } else if(cur->right != NULL && cur->right != prev){
            tracker.push(cur->right);
            tmp = tmp * 10 + cur->right->val;
        } else {
            tracker.pop();
            prev = cur;

            if(!(cur->left) && !(cur->right)){
                res += tmp;
            }
            tmp -= cur->val;
            tmp = tmp /10;
        }
    }
    
    return res;
}
/*
 http://oj.leetcode.com/problems/path-sum/
 */
bool hasPathSum(TreeNode *root, int sum) {
    int tmp = 0;
    
    if(!root) {return false;}
    
    stack<TreeNode*> tracker;
    TreeNode *cur;
    TreeNode *prev = NULL;
    tracker.push(root);
    tmp = root->val;
    while(!tracker.empty()){
        cur = tracker.top();
        if(cur->left && cur->left != prev &&
           (!prev || prev != cur->right)){
            tracker.push(cur->left);
            tmp += cur->left->val;
            prev = NULL;
        } else if(cur->right != NULL && cur->right != prev){
            tracker.push(cur->right);
            tmp += cur->right->val;
        } else {
            tracker.pop();
            prev = cur;
            
            if(!(cur->left) && !(cur->right)){
                if(sum == tmp){
                    return true;
                }
            }
            tmp -= cur->val;
        }
    }
    
    return false;
}

/*
 http://oj.leetcode.com/problems/path-sum-ii/
 */
vector<int> worker(stack<TreeNode*> s){
    vector<int> res;
    while(!s.empty()){
        
        res.push_back(s.top()->val);
        s.pop();
    }
    reverse(res.begin(), res.end());
    return res;
}

vector<vector<int> > pathSum(TreeNode *root, int sum) {
    vector<vector<int>> res;
    int tmp = 0;
    
    if(!root) {return res;}
    
    stack<TreeNode*> tracker;
    TreeNode *cur;
    TreeNode *prev = NULL;
    tracker.push(root);
    tmp = root->val;
    while(!tracker.empty()){
        cur = tracker.top();
        if(cur->left && cur->left != prev &&
           (!prev || prev != cur->right)){
            tracker.push(cur->left);
            tmp += cur->left->val;
            prev = NULL;
        } else if(cur->right != NULL && cur->right != prev){
            tracker.push(cur->right);
            tmp += cur->right->val;
        } else {
            tracker.pop();
            prev = cur;
            
            if(!(cur->left) && !(cur->right)){
                if(sum == tmp){
                    res.push_back(worker(tracker));
                }
            }
            tmp -= cur->val;
        }
    }
    
    return res;
}

/*
 http://oj.leetcode.com/problems/longest-consecutive-sequence/
 */
int longestConsecutive(vector<int> &num){
    if(num.size() == 0) return 0;
    unordered_map<int, int> tracker;
    int maxLength = 1;
    
    for(auto it : num){
        if(tracker[it] != 0){
            continue;
        }
        tracker[it] = 1;
        int left = tracker[it-1];
        int right = tracker[it+1];
        tracker[it-left] = tracker[it+right] = 1+ left + right;
        maxLength = maxLength > 1+ left + right? maxLength: 1+ left + right;
    }
    return maxLength;
}

/*
 http://oj.leetcode.com/problems/word-ladder/
 */
bool isOneCharDiff(const string s1, const string s2){
    if(s1.length() != s2.length()){
        return false;
    }
    int count = 0;
    for(int i = 0;i< s1.length(); i++)
    {
        if(s1[i] != s2[i]){
            count ++;
            if(count >1){
                return false;
            }
        }
    }
    return count ==1;
}

int ladderLength(string start, string end, unordered_set<string> &dict) {
    if(dict.size() == 0 || start == end){
        return 0;
    }
    
    queue<string> tracker;
    unordered_set<string> visited;
    int distance = 1;
    
    tracker.push(start);
    tracker.push("-1");
    visited.insert(start);
    while(!tracker.empty()){
        if(tracker.front() == "-1"){
            distance ++;
            tracker.pop();
            tracker.push("-1");
        }
        string tmp = tracker.front();
        if(isOneCharDiff(tmp, end)){
            return distance + 1;
        }
        
        for(auto it : dict){
            if(visited.find(it) == visited.end() && isOneCharDiff(it, tmp)){
                tracker.push(it);
                visited.insert(it);
            }
        }
        tracker.pop();
    }
    return 0;
}

/*
 http://oj.leetcode.com/problems/valid-palindrome/
 */
bool isAlphaNumeric(char ch){
    if((ch>=65)&&(ch<=90))
        return true;
        
    if((ch>=97) &&(ch<=122))
        return true;
            
    if((ch>=48) && (ch<=57))
        return true;
    
    return false;
}

bool isPalindrome(string s) {
    int left = 0;
    int right = (int)s.length()-1;
    
    while(left <=right){
        while(!isAlphaNumeric(s[left]) && left < s.length()-1) left++;
        while(!isAlphaNumeric(s[right]) && right >= 0) right--;
        if(left <= right){
            if(tolower(s[left]) == tolower(s[right])){
                left++;
                right--;
            }
            else
            {
                return false;
            }
        }else{
            return true;
        }
    }
    return true;
}

/*
 http://oj.leetcode.com/problems/binary-tree-maximum-path-sum/
 */
int maxPath(TreeNode* root, int& sum){
    if(!root) return INT_MIN;
    
    int left = maxPath(root->left, sum);
    int right = maxPath(root->right, sum);
    if(root->val >0){
        sum = sum > left+ right + root->val ? sum : left + right + root->val;
        sum = sum > root->val? sum: root->val;
    }
    else{
        sum = sum > left? sum : left;
        sum = sum > right? sum : right;
        sum = sum > root->val ? sum : root->val;
    }
    int max;
    if(left == INT_MIN && right == INT_MIN){
        max = root->val;
    }else{
        max = left > right? left + root->val : right + root->val;
    }
    
    return max;
}

int maxPathSum(TreeNode *root) {
    int sum = INT_MIN;
    
    int max = maxPath(root, sum);
    max = max > sum? max: sum;
    return max;
}

/*
 http://oj.leetcode.com/problems/two-sum/
 */
vector<int> twoSum(vector<int> &numbers, int target){
    unordered_map<int, int> tracker;
    for(int i = 0;i< numbers.size();i++){
        if(tracker.find(target - numbers[i]) != tracker.end())
        {
            int a =tracker.find(target-numbers[i])->second;
            int first = min(a, i);
            int second = max (a,i);
            vector<int> res;
            res.push_back(first+1);
            res.push_back(second+1);
            return res;
        }else{
            tracker.insert(make_pair(numbers[i], i));
        }
    }
    vector<int> res;
    return res;
    
}

/*
 http://oj.leetcode.com/problems/longest-substring-without-repeating-characters/
 */
int lengthOfLongestSubstring(string s) {
    if(s.length() == 0) {
        return 0;
    }
    
    bool tracker[256];
    for(int i = 0;i< 256;i++){
        tracker[i] = false;
    }
    int maxLength = 0;
    int currentLength = 0;
    
    int startindex = 0;
    int endindex = 0;
    
    while(endindex < (int)s.length()){
        if(tracker[s[endindex]] == false){
            tracker[s[endindex]] = true;
            currentLength ++;
            maxLength = maxLength > currentLength ? maxLength : currentLength;
            endindex ++;
        }else{
            while(s[startindex] != s[endindex] && startindex < endindex){
                tracker[s[startindex]] = false;
                startindex++;

            }
            tracker[s[endindex]] = true;
            startindex ++;
            endindex ++;

            currentLength = endindex - startindex;
        }
    }
    return maxLength;
}

/*
 http://oj.leetcode.com/problems/add-two-numbers/
 */
ListNode *addTwoNumbers(ListNode *l1, ListNode *l2) {
    ListNode *dummyHead = new ListNode(-1);
    ListNode *tmp = dummyHead;
    int carry = 0;
    int res = 0;
    
    while(l1 && l2){
        res = l1->val + l2->val + carry;
        carry = res  / 10;
        res = res %10;
        ListNode *cur = new ListNode(res);
        tmp->next = cur;
        tmp = cur;
        l1 = l1->next;
        l2 = l2->next;
    }
    
    if(!l1 && !l2){
        if(carry != 0){
            tmp->next = new ListNode(carry);
        }
    }else{
        ListNode *rest;
        if(l1){ rest = l1;}
        if(l2){ rest = l2;}
        while(rest){
            if(carry == 0){break;}
            res = rest->val + carry;
            carry = res / 10;
            res = res % 10;
            ListNode *cur = new ListNode(res);
            tmp->next = cur;
            tmp = tmp->next;
            rest =rest->next;
        }
        if(carry == 0){
            tmp->next = rest;
        }else{
            tmp->next = new ListNode(carry);
        }
    }

    return dummyHead ->next;
}

/*
 http://oj.leetcode.com/problems/longest-palindromic-substring/
 */
string longestPalindrome(string s) {
    if(s.length() <1){return s;}
    int maxlength = 1;
    string res;
    
    for(int i = 1;i<s.length(); i++){
        int left= i;
        int right = i;
        while(left >=0 && right <s.length() && s[left] == s[right]){
            if(right - left + 1 > maxlength){
                maxlength = right-left +1;
                res = s.substr(left, maxlength);
            }
            left--;
            right++;
        }
        left = i -1;
        right = i;
        while(left >=0 && right <s.length() && s[left] == s[right]){
            if(right - left + 1 > maxlength){
                maxlength = right-left +1;
                res = s.substr(left, maxlength);
            }
            left--;
            right++;
        }
    }
    return res;
}

/*
 http://oj.leetcode.com/problems/reverse-integer/
 */
int reverse(int x) {
    int res = 0;
    int tmp;
    int sign = x >=0? 1: -1;
    x = x>=0?x:-x;
    while(x>0){
        tmp = x %10;
        x = x/10;
        res = res * 10 + tmp;
    }
    res = res * sign;
    
    return res;
}

/*
 http://oj.leetcode.com/problems/best-time-to-buy-and-sell-stock/
 */
int maxProfit(vector<int> &prices) {
    if(prices.size() < 1) return 0;
    int max = 0;
    int cur;
    bool buying = false;
    
    for(int i = 1;i< prices.size();i++){
        if(buying){
            cur += prices[i] - prices[i-1];
            if(cur <0){
                buying = false;
                cur = 0;
            }
        }else if(prices[i] > prices[i-1]){
            buying = true;
            cur = prices[i] - prices[i-1];
        }
        
        max = max > cur? max: cur;
    }
    return max;
}

/*
 http://oj.leetcode.com/problems/best-time-to-buy-and-sell-stock-ii/
 */
int maxProfitII(vector<int> &prices) {
    if(prices.size() == 0) return 0;
    int max = 0;
    int cur = 0;
    for(int i = 0;i< prices.size()-1; i ++)
    {
        if(prices[i+1]> prices[i]){
            cur += prices[i+1]- prices[i];
        }
        
        if(i == prices.size() -2 || prices[i+1] < prices[i]){
            max += cur;
            cur = 0;
        }
    }
    return max;
}

/*
 http://oj.leetcode.com/problems/minimum-depth-of-binary-tree/
 */
int minDepth(TreeNode *root) {
    if(!root) return 0;
    if(!root->left && !root->right) return 1;

    int res = 0;
    if(root->left && root->right){
        int left = minDepth(root->left);
        int right = minDepth(root->right);
        res = min(left, right);
    } else if(root->left){
        res =minDepth(root->left);
    } else {
        res = minDepth(root->right);
    }
    
    return 1 + res;
}
/*
http://oj.leetcode.com/problems/sqrtx/
 */
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

bool IsSpace(char c)
{
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') { return true; } else { return false;}
}

bool IsDigit(char c)
{
    if (c - '0' <= 9) { return true; } else { return false;}
}


int atoi(const char *str) {
    long res = 0; int sign = 1; int c = 0;
    
    int i = 0;
    
    while (str[i] != '\0' && IsSpace(str[i]))
    {
        i++;
    }
    
    if (str[i] == '+' || str[i] == '-')
    {
        
        if (str[i] == '-')
        {
            sign = -1;
        }
        i++;
    }
    
    while (str[i] != '\0')
    {
        if (IsDigit(str[i]))
        {
            //throw new FormatException(str + "is not a valid number"); }
            c = str[i] - '0';
            if (sign > 0 && (res >LONG_MAX / 10 || (res == LONG_MAX / 10 && c > LONG_MAX % 10)))
            {
                //throw new OverflowException(input + "is larger than maximum long");
            }else if (sign < 0 && (res > LONG_MAX / 10 || (res == LONG_MAX / 10 && c > LONG_MAX % 10 + 1)))
            {
                //throw new OverflowException(input + "is less than minmum long.");
            }
            res = (res * 10) + c;
            i++;
        }
        else{
            break;
        }


    }
    res = sign > 0 ? res : -res;
    if(res > (long)INT_MAX) res = INT_MAX;
    if(res < (long)INT_MIN) res = INT_MIN;
    return (int)res;
}

/*
 http://oj.leetcode.com/problems/palindrome-number/
 */
bool isPalindrome(int x) {
    long origin = x;
    long tmp = 0;
    int sign = 1;
    if(origin <0){
        sign = -1;
        origin = -origin;
    }
    while(origin >0){
        int i = origin %10;
        origin = origin/10;
        tmp = tmp * 10 + i;
    }
    return (int)tmp == x;
}

/*
 http://oj.leetcode.com/problems/regular-expression-matching/
 */
bool isMatch(const char *s, const char *p) {
    if( *s == '\0'){
        if(*p == '\0')
        {
            return true;
        }else if(*(p+1) == '*')
        {
            return isMatch(s, p+2);
        }
        return false;
    }else if(*p == '\0')
    {
        return false;
    }
    
    if((*p == '.' && *(p+1) != '*') || (*s == *p && *(p+1) != '*')){
        return isMatch(s+1, p+1);
    }
    
    if( *(p+1) == '*') {
        if( *p == *s || *p == '.'){
            return isMatch(s+1, p) || isMatch(s, p+2);
        }else{
            return isMatch(s, p+2);
        }
    }
    
    return false;
}

/*
 http://oj.leetcode.com/problems/integer-to-roman/
 */
string intToRoman(int num) {
    struct romandata_t { int value; char const* numeral; };
    static romandata_t const romandata[] =
    { 1000, "M",
        900, "CM",
        500, "D",
        400, "CD",
        100, "C",
        90, "XC",
        50, "L",
        40, "XL",
        10, "X",
        9, "IX",
        5, "V",
        4, "IV",
        1, "I",
        0, NULL }; // end marker

    std::string result;
    for (romandata_t const* current = romandata; current->value > 0; ++current)
    {
        while (num >= current->value)
        {
            result += current->numeral;
            num  -= current->value;
        }
    }
    return result;
}

/*
 http://oj.leetcode.com/problems/roman-to-integer/
 */
int romanToInt(string s) {
    struct romandata_t { char const* numeral;int value;  };
    unordered_map<string, int> tracker;
    tracker.insert(make_pair("M",1000));
    tracker.insert(make_pair("CM",900));
    tracker.insert(make_pair("D",500));
    tracker.insert(make_pair("CD",400));
    tracker.insert(make_pair("C",100));
    tracker.insert(make_pair("XC",90));
    tracker.insert(make_pair("L",50));
    tracker.insert(make_pair("XL",40));
    tracker.insert(make_pair("X",10));
    tracker.insert(make_pair("IX",9));
    tracker.insert(make_pair("V",5));
    tracker.insert(make_pair("IV",4));
    tracker.insert(make_pair("I",1));
    tracker.insert(make_pair("\0",0));
    
    int res = 0;
    int i = 0;
    while(i<s.length()-1){
        if(tracker.find(s.substr(i, 2)) != tracker.end()){
            res += tracker.find(s.substr(i,2))->second;
            i +=2;
        }else if(tracker.find(s.substr(i, 1)) != tracker.end()){
            res += tracker.find(s.substr(i, 1))->second;
            i++;
        }
    }
    if(i == s.length() -1){
        if(tracker.find(s.substr(s.length()-1,1)) != tracker.end()){
            res +=tracker.find(s.substr(s.length()-1,1))->second;
        }
    }
    return res;

}

/*
 http://oj.leetcode.com/problems/longest-common-prefix/
 */
struct TrieNode{
    TrieNode *val[256]; bool isEnd[256];
};

string longestCommonPrefix(vector<string> &strs) {
    if(strs.size() == 1){return strs[0];}
    if(strs.size() == 0){return "";}
    TrieNode *root = new TrieNode();
    TrieNode *tmp = root;
    int maxlength = INT_MAX;
    int index = 0;
    for(int i = 0;i< strs.size(); i++){
        index = 0;
        tmp = root;
        string s = strs[i];
        maxlength = maxlength < (int)s.length()? maxlength : (int)s.length();
        if(s.length() == 0) {return s;}
        while(index != strs[i].length()){
            if(!tmp->val[strs[i][index]]){
                tmp->val[strs[i][index]] = new TrieNode();
            }
            tmp = tmp->val[strs[i][index]];
            index++;
        }
    }
    
    tmp = root;
    string res;
    while(tmp){
        int count = 0;
        int index = -1;
        for(int i = 0;i < 256; i++)
        {
            if(tmp->val[i]){
                count++;
                index = i;
            }
        }
        if(count==1 && res.length() < maxlength){
            res += (char)index;
            tmp = tmp->val[index];
        }else{
            break;
        }
        
    }
    return res;
}

/*
 http://oj.leetcode.com/problems/3sum/
 */
vector<vector<int> > threeSum(vector<int> &num) {
    vector<vector<int>> tracker;
    vector<vector<int>> res;
    sort(num.begin(), num.end(), [](int a, int b){return a<b;});
    
    for(int i= 0;i< num.size();i++){
        int left = 0;
        int right = (int)(num.size()) -1;
        while(left< right){
            if(left == i){
                left ++;
                continue;
            }
            if(right == i){
                right--;
                continue;
            }
            if(num[i] + num[left] + num[right] == 0){
                vector<int> x = {i, left, right};
                sort(x.begin(), x.end(), [](int a, int b)->bool {return a<b;});
                bool contain = false;
                for(vector<vector<int>>::iterator it = tracker.begin(); it != tracker.end(); it++){
                    if((*it)[0] == x[0] && (*it)[1] == x[1] && (*it)[2] == x[2]){
                        contain = true;
                        break;
                    }
                }
                if(!contain){
                    tracker.push_back(x);
                }
                left++;
                right--;
                
            }else if(num[i] + num[left] + num[right] > 0){
                right--;
            }else {
                left ++;
            }
        }
    }
    
    for(vector<vector<int>>::iterator it = tracker.begin(); it != tracker.end(); it++){
        vector<int> tmp = {num[(*it)[0]], num[(*it)[1]], num[(*it)[2]]};
        res.push_back(tmp);
    }
    
    return res;

}

/*
 http://oj.leetcode.com/problems/remove-nth-node-from-end-of-list/
 */

ListNode *removeNthFromEnd(ListNode *head, int n) {
    ListNode* front = head;
    ListNode *behind = head;
    for(int i = 0;i<n;i++){
        front = front->next;
    }
    if(!front) return head->next;
    while(front->next){
        front = front->next;
        behind = behind ->next;
    }
    behind ->next= behind->next->next;
    return head;
}

/*
 http://oj.leetcode.com/problems/letter-combinations-of-a-phone-number/
 */
void worker(unordered_map<char, string> tracker, vector<string>& res,const string digits, int level, string tmp){
    if(level == digits.length()){
        res.push_back(tmp);
        return;
    }
    
    string letters = tracker.find(digits[level])->second;
    for(int j = 0; j < letters.length(); j++){
        tmp = tmp + letters[j];
        worker(tracker, res, digits, level+1, tmp);
        tmp = tmp.substr(0, tmp.length()-1);
    }

}
vector<string> letterCombinations(string digits) {
    unordered_map<char, string> numbermap;
    numbermap.insert(make_pair('1', ""));
    numbermap.insert(make_pair('2', "abc"));
    numbermap.insert(make_pair('3', "def"));
    numbermap.insert(make_pair('4', "ghi"));
    numbermap.insert(make_pair('5', "jkl"));
    numbermap.insert(make_pair('6', "mno"));
    numbermap.insert(make_pair('7', "pqrs"));
    numbermap.insert(make_pair('8', "tuv"));
    numbermap.insert(make_pair('9', "wxyz"));
    numbermap.insert(make_pair('0', " "));
    
    vector<string> res;
    worker(numbermap, res, digits, 0, "");
    return res;
}

/*
 http://oj.leetcode.com/problems/valid-parentheses/
 */
bool isValid(string s) {
    unordered_map<char, char> mapper;
    mapper.insert(make_pair('}', '{'));
    mapper.insert(make_pair(')', '('));
    mapper.insert(make_pair(']', '['));
    
    stack<char> tracker;
    
    for(int i = 0;i<s.length();i++){
        if(s[i] == '{' || s[i] == '(' || s[i] == '['){
            tracker.push(s[i]);
        }else {
            if(!tracker.empty()){
                char tmp = tracker.top();
                if(tmp == mapper.find(s[i])->second){
                    tracker.pop();
                }else{
                    return false;
                }
            }else{
                return false;
            }
        }
    }
    if(tracker.empty()){
        return true;
    }else{
        return false;
    }
}

/*
 http://oj.leetcode.com/problems/generate-parentheses/
 */
void worker(vector<string>& res, int n , int left, int right, string tmp){
    if(left == n && right == n){
        res.push_back(tmp);
        return;
    }
    
    if(left <n){
        worker(res,n, left+1,right,tmp+'(');
    }
    if(right < left){
        worker(res, n, left, right+1, tmp+')');
    }
    
}
vector<string> generateParenthesis(int n) {
    vector<string> res;
    worker(res, n, 0,0,"");
    return res;
}

/*
 http://oj.leetcode.com/problems/merge-k-sorted-lists/
 */
ListNode* getSmallest(vector<ListNode*>& lists){
    int min = INT_MAX;
    int index = -1;
    for(int i = 0;i< lists.size(); i++){
        if(lists[i] != NULL && min > lists[i]->val){
            min = lists[i]->val;
            index = i;
        }
    }
    ListNode* tmp = NULL;
    if(index != -1){
        tmp = lists[index];
        lists[index] = lists[index] -> next;
    }
    return tmp;
}

ListNode *mergeKLists(vector<ListNode *> &lists) {
    ListNode* dummyHead = new ListNode(0);
    ListNode* tmp = dummyHead;
    ListNode* smallest = getSmallest(lists);
    
    while(smallest){
        tmp->next = smallest;
        tmp = tmp->next;
        smallest = getSmallest(lists);
    }
    return dummyHead->next;
}

/*
 http://oj.leetcode.com/problems/swap-nodes-in-pairs/
 */
ListNode *swapPairs(ListNode *head) {
    ListNode* dummyHead = new ListNode(0);
    dummyHead->next = head;
    
    ListNode* firstNode = dummyHead;
    ListNode* lastNode = NULL;
    ListNode* nextHead = NULL;
    
    while(firstNode->next){
        if(firstNode->next){
            lastNode = firstNode->next->next;
        }
        if(lastNode){
            nextHead = lastNode->next;
        }
        
        if(!lastNode){
            break;
        }
        
        lastNode->next = firstNode->next;
        firstNode->next= lastNode;
        lastNode->next->next = nextHead;
        firstNode = firstNode->next->next;
        lastNode = NULL;
        nextHead = NULL;
    }
    return dummyHead->next;
    
}

/*
 http://oj.leetcode.com/problems/remove-duplicates-from-sorted-array/
 */
int removeDuplicates(int A[], int n) {
    if(A==NULL || n == 0) return NULL;
    if(n == 1) return 1;
    int index = 0;
    int runner = 1;
    while(runner < n){
        if(A[index] != A[runner]){
            A[index+1] = A[runner];
            index ++;
            runner++;
        }else{
            runner++;
        }
    }
    return index+1;
}

/*
 http://oj.leetcode.com/problems/remove-element/
 */
int removeElement(int A[], int n, int elem) {
    if(A==NULL || n == 0) return NULL;
    int index = -1;
    for(int i = 0;i< n;i++){
        if(A[i] != elem){
            A[index+1] = A[i];
            index++;
        }
    }
    return index+1;
}

/*
 http://oj.leetcode.com/problems/divide-two-integers/
 */
int divide(int dividend, int divisor) {
    long a = abs(dividend);
    long b = abs(divisor);
    
    long c = b;
    int res = 0;
    while(a>=b){
        c = b;
        for(int i = 0;a>=c;i++, c<<=1){
            res += 1<<i;
            a -= c;
        }
    }
    
    return((dividend^divisor)>>31)? (-res): res;
}

/*
 http://oj.leetcode.com/problems/next-permutation/
 */
void nextPermutation(vector<int> &num) {
    int lastindex = (int)num.size()-1;
    int firstindex = -1;
    while(lastindex >0 &&  num[lastindex-1] >= num[lastindex]){
        lastindex --;
    }
    
    if(lastindex == 0){
        reverse(num.begin(), num.end());
        return;
    }
    
    firstindex = lastindex -1;
    while(lastindex < num.size() && num[firstindex]<num[lastindex] ){
        lastindex++;
    }
    lastindex--;
    
    int tmp = num[lastindex];
    num[lastindex] = num[firstindex];
    num[firstindex] = tmp;
    
    int left = lastindex + 1;
    int right = (int)num.size()-1;
    while(left < right){
        int tmp = num[left];
        num[left] = num[right];
        num[right] = tmp;
        left++;
        right--;
    }
}

/*
 http://oj.leetcode.com/problems/triangle/
 */
int minimumTotal(vector<vector<int> > &triangle) {
    if(triangle.size() == 0) return 0;
    vector<vector<int>>::iterator it = triangle.end() -1;
    //vector<int> tracker(*it);
    vector<int> tracker;
    for(int i = 0;i< (*it).size();i++){
        tracker.push_back((*it)[i]);
    }
    it --;
    while(it != triangle.begin()-1){

        for(int i = 0;i< (*it).size(); i++)
        {
            tracker[i] = tracker[i] < tracker[i+1]? (*it)[i] + tracker[i] : (*it)[i] + tracker[i+1];
        }
        it --;
    }
    return tracker[0];
}

/*
 http://oj.leetcode.com/problems/pascals-triangle/
 */
vector<vector<int> > generate(int numRows) {
    vector<vector<int>> res;
    for(int i = 0;i< numRows;i++){
        if(i == 0){
            vector<int> tmp = {1};
            res.push_back(tmp);
        }else if(i == 1){
            vector<int>tmp = {1,1};
            res.push_back(tmp);
        }else{
            vector<int>tmp = {1};
            vector<int> last = res[i-1];
            for(int j = 0;j<i-1;j++){
                tmp.push_back(last[j] + last[j+1]);
            }
            tmp.push_back(1);
            res.push_back(tmp);
        }

    }
    return res;
}

/*
 http://oj.leetcode.com/problems/pascals-triangle-ii/
 */
vector<int> getRow(int rowIndex) {
    
    vector<vector<int>> res;
    for(int i = 0;i< rowIndex;i++){
        if(i == 0){
            vector<int> tmp = {1};
            res.push_back(tmp);
        }else if(i == 1){
            vector<int>tmp = {1,1};
            res.push_back(tmp);
        }else{
            vector<int>tmp = {1};
            vector<int> last = res[i-1];
            for(int j = 0;j<i-1;j++){
                tmp.push_back(last[j] + last[j+1]);
            }
            tmp.push_back(1);
            res.push_back(tmp);
        }
        
    }
    return res[rowIndex-1];
}

/*
 http://oj.leetcode.com/problems/populating-next-right-pointers-in-each-node-ii/
 */
struct TreeLinkNode {
    int val;
    TreeLinkNode *left, *right, *next;
    TreeLinkNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
};

void worker(TreeLinkNode *root, TreeLinkNode *parent){
    if(!root) return;
    if(parent){
        if(root == parent->left && parent->right){
            root->next = parent->right;
        }else{
            if(parent->next){
                if(parent->next->left){
                    root->next = parent->next->left;
                }else{
                    root->next = parent->next->right;
                }
            }
        }
    }
    worker(root->right, root);
    worker(root->left, root);
}

void connect(TreeLinkNode *root) {
    worker(root, NULL);
}

/*
 http://oj.leetcode.com/problems/balanced-binary-tree/
 */
bool depth(TreeNode *root, int& maxdepth){
    if(!root){
        maxdepth = 0;
        return true;
    }
    int left = 0;
    int right = 0;
    bool lefttrue = depth(root->left, left);
    bool righttrue = depth(root->right, right);
    maxdepth = left>right? left+1: right+1;
    return abs(left-right)<=1 && lefttrue && righttrue;
}

bool isBalanced(TreeNode *root) {
    int maxdepth = 0;
    return depth(root, maxdepth);
}

/*
 http://oj.leetcode.com/problems/convert-sorted-list-to-binary-search-tree/
 */
TreeNode *sortedListToBST(ListNode *head) {
    if(!head) return NULL;
    ListNode* slow = head;
    ListNode* fast = head->next;
    if(!fast) {return new TreeNode(head->val);}
    fast = fast->next;
    while(fast && fast->next){
        fast = fast->next->next;
        slow = slow->next;
    }
    
    ListNode* newHead = slow->next->next;
    ListNode* treeHeadNote = slow->next;
    slow->next = NULL;
    TreeNode *root = new TreeNode(treeHeadNote->val);
    root->left = sortedListToBST(head);
    root->right = sortedListToBST(newHead);
    return root;
    
    
}

/*
 http://oj.leetcode.com/problems/convert-sorted-array-to-binary-search-tree/
 */
TreeNode* worker(vector<int>& num, int left, int right){
    //if(left == right) return new TreeNode(num[left]);
    if(left > right) return NULL;
    
    int mid = left + (right-left) / 2;
    
    TreeNode *root = new TreeNode(num[mid]);
    root->left = worker(num, left, mid-1);
    root->right = worker(num, mid+1, right);
    return root;
    
}

TreeNode *sortedArrayToBST(vector<int> &num) {
    return worker(num, 0, (int)num.size()-1);
}

/*
 http://oj.leetcode.com/problems/binary-tree-level-order-traversal-ii/
 */
vector<vector<int> > levelOrderBottom(TreeNode *root) {
    vector<vector<int>> res;
    if(!root) return res;
    queue<TreeNode*> q;
    stack<vector<int>*> s;
    q.push(root);
    q.push(NULL);
    vector<int> *tmp = new vector<int>();
    TreeNode* topNode;
    while(!q.empty()){
        topNode = q.front();
        if(!topNode){
            s.push(tmp);
            q.pop();
            if(!q.empty()){
                tmp = new vector<int>();
                q.push(NULL);
            }else{
                break;
            }
        }else{
            tmp->push_back(topNode->val);
            q.pop();
            if(topNode->left) q.push(topNode->left);
            if(topNode->right) q.push(topNode->right);
        }
    }

    while(!s.empty()){
        res.push_back(*(s.top()));
        s.pop();
    }
    return res;
}

/*
 http://oj.leetcode.com/problems/maximum-depth-of-binary-tree/
 */
int maxDepth(TreeNode *root) {
    if(!root) return 0;
    int left = maxDepth(root->left);
    int right = maxDepth(root->right);
    return left> right? left+1: right+1;
}

/*
 http://oj.leetcode.com/problems/binary-tree-zigzag-level-order-traversal/
 */
vector<vector<int> > zigzagLevelOrder(TreeNode *root) {
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

/*
 http://oj.leetcode.com/problems/binary-tree-level-order-traversal/
 */
vector<vector<int> > levelOrder(TreeNode *root) {
    vector<vector<int>> res;
    if(!root) return res;
    queue<TreeNode*> tracker;
    TreeNode *tmp;
    tracker.push(root);
    tracker.push(NULL);
    vector<int>* container = new vector<int>();
    while(!tracker.empty()){
        tmp = tracker.front();
        if(!tmp){
            res.push_back(*container);
            tracker.pop();
            if(!tracker.empty()){
                tracker.push(NULL);
                container = new vector<int>();
            }else{
                break;
            }
        }else{
            container->push_back(tmp->val);
            if(tmp->left) tracker.push(tmp->left);
            if(tmp->right) tracker.push(tmp->right);
            tracker.pop();
        }
    }
    return res;
}

/*
 http://oj.leetcode.com/problems/symmetric-tree/
 */
bool isSymmetric(TreeNode* left, TreeNode* right){
    if(!left && !right) return true;
    if(!left || !right) return false;
    if(left->val != right->val) return false;
    
    return isSymmetric(left->left, right->right) && isSymmetric(left->right, right->left);
    
}

bool isSymmetric(TreeNode *root) {
    if(!root) return true;
    return isSymmetric(root->left, root->right);
    
}

/*
 http://oj.leetcode.com/problems/same-tree/
 */
bool isSameTree(TreeNode *p, TreeNode *q) {
    if(!p && !q) return true;
    if(!p || !q) return false;
    if(p->val != q->val) return false;
    
    return isSameTree(p->left, q->left) && isSameTree(q->right, q->right);
}

/*
 http://oj.leetcode.com/problems/validate-binary-search-tree/
 */
bool isValidBST(TreeNode *root, int min, int max){
    if(!root) return true;
    if(root->val <= min || root->val >= max) return false;
    return isValidBST(root->left, min, root->val) && isValidBST(root->right, root->val, max);
}

bool isValidBST(TreeNode *root) {
    return isValidBST(root, INT_MIN, INT_MAX);
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
    //string s= "1 ";
    //reverseWords(s);
    //string input = "aab";
    //vector<vector<string>> res = partition(input);
    
    //for(vector<vector<string>> :: iterator it = res.begin(); it != res.end(); it ++){
    //    cout<<endl;
    //    for(vector<string>::iterator stit = (*it).begin(); stit != (*it).end(); stit++){
    //        cout << *stit<<endl;
    //    }
    //}
    
    /*
    string input = "aab";
    int res = minPalindromeCut(input);
    cout<<res<<endl;
    
    string s = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    res = minPalindromeCut(s);
    cout<<res<<endl;
    */
    /*
    vector<char> r1 = {'O','X','O'};
    vector<char> r2 = {'X','O','X'};
    vector<char> r3 = {'O','X','O'};
//    vector<char> r4 = {'X','O','X','X'};
    
    vector<vector<char>> board = {r1,r2,r3};
    solve(board);
    
    vector<char> rr1 = {'O','O','O'};
    vector<char> rr2 = {'O','O','O'};
    vector<char> rr3 = {'O','O','O'};
    
    board = {rr1,rr2,rr3};
    solve(board);
    */
    
    //TreeNode* r = new TreeNode(9);
    //TreeNode* right = new TreeNode(1);
    //r->right = right;
    //int res =sumNumbers(r);
    //cout<<res<<endl;
    
    /*
    int res = ladderLength("nanny", "aloud", input);
    cout << res;
     */
    
    /*string input = "A man, a plan, a canal: Panama";
    bool res = isPalindrome(input);
    cout<< res<<endl;;
    
    input = " ";
    res = isPalindrome(input);
    cout<<res<<endl;
    */
    
    /*
     TreeNode *t1 = new TreeNode(-3);
    int res = maxPathSum(t1);
    cout << res << endl;
    
    t1->val = 0;
    TreeNode *t2 = new TreeNode(1);
    TreeNode *t3 = new TreeNode(1);
    t1->left = t2;
    t2->left = t3;
    
    res = maxPathSum(t1);
    cout<< res <<endl;
    */
    
    /*
    string s = "abbcd";
    int res = lengthOfLongestSubstring(s);
    cout<<res<<endl;
    
    s = "qopubjguxhxdipfzwswybgfylqvjzhar";
    res = lengthOfLongestSubstring(s);
    cout<<res<<endl;
    */
    /*
    ListNode *l1 = new ListNode(9);
    l1->next = new ListNode(9);
    ListNode *l2 = new ListNode(1);
    ListNode *res = addTwoNumbers(l1, l2);
    cout << res->val <<endl;
    cout << res->next->val<<endl;
    */
    /*
    vector<int> input = {1,2,4};
    int res = maxProfitII(input);
    cout<< res <<endl;
     */
    /*int input = 2;
    int res = sqrt(input);
    cout<< res<<endl;
     */
    //int res = atoi("2147483648");
    //cout<<res;
    /*bool res = isMatch("aa", ".*c");
    cout <<res<<endl;
    res = isMatch("aab", "c*a*b");
    cout <<res<<endl;
    res = isMatch("aaa", "a*a");
    cout <<res<<endl;
    res = isMatch("a", ".*a*a");
    cout <<res<<endl;
    res = isMatch("", "c*c*");
    cout <<res<<endl;
    */
    //std::cout << intToRoman(123) << std::endl;
    
    //for (int i = 1; i <= 4000; ++i)
    //{
    //    std::cout << intToRoman(i) << std::endl;
    //}
    /*
    vector<string> input = {"", ""};
    string res = longestCommonPrefix(input);
    cout << res <<endl;
    */
    
    /*vector<int> input = {0,0,0};
    vector<vector<int>> res = threeSum(input);
     */
    /*vector<string> res = letterCombinations("22");
    for(vector<string>::iterator it = res.begin(); it != res.end(); it++){
        cout<< (*it)<<endl;
    }
     */
    /*ListNode* t1 = new ListNode(1);
    ListNode* t2 = new ListNode(2);
    t1->next = t2;
    ListNode* res = swapPairs(t1);
     */
    /*
    TreeNode* root = new TreeNode(1);
    vector<vector<int>> res = pathSum(root, 1);
     */
    
    /*
    int res = divide(2147483647, 1);
    cout<< res <<endl;
    
    res = divide(16, 4);
    cout<< res <<endl;
    
    res = divide(16, 3);
    cout<< res <<endl;
    
    res = divide(2147483647, 3);
    cout<< res <<endl;
    */
    /*
    vector<int> input = {1,2};
    nextPermutation(input);
     */
    /*vector<int> t1 = {1};
    vector<int> t2 = {2,3};
    vector<vector<int>> input = {t1,t2};
    int res = minimumTotal(input);
    cout<<res<<endl;
     */
    //vector<vector<int>> res = generate(3);
    /*
    TreeNode *t1 = new TreeNode(1);
    TreeNode *t2 = new TreeNode(2);
    TreeNode *t3 = new TreeNode(3);
    t1->right = t2;
    t2->right = t3;
    
    bool res = isBalanced(t1);
    cout<<res<<endl;
    
    t1->left = t2;
    t2->right = NULL;
    t1->right = t3;
    res = isBalanced(t1);
    cout<<res<<endl;
    */
    
    return 0;
}