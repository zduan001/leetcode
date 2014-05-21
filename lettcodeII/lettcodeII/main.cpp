//
//  main.cpp
//  lettcodeII
//
//  Created by Duan, David on 4/24/14.
//  Copyright (c) 2014 Duan, David. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <stack>
#include <assert.h>

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

int findKthElement(int A[], int m, int B[], int n, int k)
{
    if(m>n)
        return findKthElement(B, n, A,m, k);
    if(m== 0)
        return B[k-1];
    if(k == 1)
        return min(A[0], B[0]);
    
    int ia = min(m, k/2);
    int ib = k-ia;
    
    if(B[ib-1] > A[ia-1])
    {
        return findKthElement( A + ia, m-ia, B, n, k-ia);
    }
    else if(B[ib-1] < A[ia-1])
    {
        return findKthElement( A, m, B+ib, n-ib, k- ib);
    }
    else
    {
        return A[ia-1];
    }
}

double findMedianSortedArrays(int A[], int m, int B[], int n) {
    if((n+m)%2 == 0){
        return (findKthElement(A, m, B, n, (n+m)/2) +
                findKthElement(A, m, B, n, (n+m)/2 +1)) /2.0;
    }else{
        return (findKthElement(A, m, B, n, (n+m)/2+1));
    }
}

// assume there exactly one solution
vector<int> twoSum(vector<int> &n, int target) {
    unordered_map<int,int> tracker;
    vector<int> res;
    for(int i = 0;i< (int)n.size();i++){
        tracker.insert(make_pair(n.at(i), i));
    }
    
    for(int i = 0;i<(int)n.size();i++){
        auto it = tracker.find(target-n.at(i));
        if(it != tracker.end())
        {
            if(i != it->second)
            {
                res.push_back(i+1);
                res.push_back(it->second+1);
                break;
            }

        }
    }
    return res;
}

int lengthOfLongestSubstringII(string s) {
    if(s.length() == 0) return 0;
    int max = 1;
    int length = 1;
    int start = 0;
    int end = 1;
    bool tracker[256];
    for(int i = 0;i< 256;i++){
        tracker[i] = false;
    }
    tracker[s.at(0)] = true;
    
    while(end != s.length()){
        if(!tracker[s.at(end)])
        {
            length++;
            max = max > length? max : length;
            tracker[s.at(end++)] = true;
        } else {
            while(s.at(start) != s.at(end)){
                tracker[s.at(start++)] = false;
                length--;
            }
            start++;
            end++;
        }
    }
    return max;
}


ListNode *addTwoNumbers(ListNode *l1, ListNode *l2) {
    
    ListNode* dummyHead = new ListNode(0);
    int carry = 0;
    ListNode *tmp = dummyHead;
    
    while(l1 || l2)
    {
        int sum = carry;
        if(l1)
            sum += l1->val;
        if(l2)
            sum += l2->val;

        carry = sum /10;
        sum = sum %10;
        
        ListNode* nextNode = new ListNode(sum);
        tmp->next = nextNode;
        tmp = nextNode;
        l1 = l1? l1->next: l1;
        l2 = l2 ? l2->next: l2;
    }
    
    if(carry > 0)
    {
        ListNode* lastNode = new ListNode(carry);
        tmp->next = lastNode;
    }
    
    return dummyHead->next;
    
}

string longestPalindrome(string s) {
    string res = "";
    int left = -1;
    int right = -1;
    
    for(int i = 0;i< s.length();i++){
        left = i;
        right = i;
        while(left>=0 && right <s.length()){
            if(s[left] == s[right])
            {
                if(right-left +1> res.length()){
                    res = s.substr(left, right-left+1);
                }
                left --;
                right ++;
            }
            else
            {
                break;
            }
        }
        left = i;
        right = i+1;
        while(left>=0 && right <s.length()){
            if(s[left] == s[right])
            {
                if(right-left +1> res.length()){
                    res = s.substr(left, right-left+1);
                }
                left --;
                right ++;
            }
            else
            {
                break;
            }
        }
    }
    return res;
}

void flattenI(TreeNode *root) {
    stack<TreeNode*> s;
    s.push(root);
    
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        s.pop();
        if(tmp->right)
            s.push(tmp->right);
        if(tmp->left)
            s.push(tmp->left);
        
        tmp->left = NULL;
        if(!s.empty()){
            tmp->right = s.top();
        }
            
    }
}

void worker(TreeNode* root, TreeNode*& tail)
{
    if(!(root->left) && !(root->right))
    {
        tail = root;
        return;
    }
    TreeNode* left =NULL;
    TreeNode* middle = NULL;
    if(root->left){
        left = root->left;
        worker(root->left, middle);
    }
    
    TreeNode* end = NULL;
    TreeNode* right =NULL;
    if(root->right){
        right = root->right;
        worker(root->right, end);
    }
    root->left = NULL;
    if(left){
        root->right = left;
        if(right){
            middle->right = right;
            tail = end;
        }
        else
        {
            tail = middle;
        }
    }
}

void flatten(TreeNode* root){
    TreeNode* tmp = NULL;
    worker(root, tmp);
    
}

int divide(int dividend, int divisor) {
    long long a = (long) abs(dividend);
    long long b = (long) abs(divisor);
    int res = 0;
    while(a>=b)
    {
        long c = b;
        for(int i = 0; a>=c; i++, c <<= 1)
        {
            a -=c;
            res += 1<<i;
        }
    }
    return (dividend & divisor) >> 31? -res: res;
}

void worker (string&s, int left, int right){
    while(left< right){
        char c = s[left];
        s[left] = s[right];
        s[right] = c;
        left ++;
        right--;
    }
    return;
}

void reverseWords(string &s) {
    worker(s, 0, (int)s.length()-1);
    
    int left = 0;
    string res = "";
    stack<char> tracker;
    while(left < s.length()){
        while(left <s.length() && s[left] == ' ')
        {
            left++;
        }
        
        while(left < s.length() && s[left] != ' '){
            tracker.push(s[left]);
            left++;
        }
        if(!tracker.empty())
        {
            if(res.length() > 0){
                res += " ";
            }
            while(!tracker.empty())
            {
                res += tracker.top();
                tracker.pop();
            }
        }
    }
    s = res;
}

int evalRPN(vector<string> &tokens) {
    stack<int> tracker;
    for( int i = 0;i< tokens.size(); i ++){
        if(tokens[i] != "+" &&
           tokens[i] != "-" &&
           tokens[i] != "*" &&
           tokens[i] != "/")
        {
            tracker.push(atoi(tokens[i].c_str()));
            continue;
        }
        
        int op2 = tracker.top();
        tracker.pop();
        int op1 = tracker.top();
        tracker.pop();
        
        if( tokens[i] == "+")
        {
            tracker.push(op1 + op2);
        }else if(tokens[i] == "-")
        {
            tracker.push(op1 - op2);
        }else if(tokens[i] == "*")
        {
            tracker.push(op1 * op2);
        }else if(tokens[i] == "/")
        {
            tracker.push(op1 / op2);
        }
    }
    return tracker.top();
}

struct TreeLinkNode {
    int val;
    TreeLinkNode *left, *right, *next;
    TreeLinkNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
};

void connect(TreeLinkNode *root) {
    if(!root) return;
    if(root->left && root->right)
    {
        root->left->next = root->right;
    }
    
    TreeLinkNode* tmp = NULL;
    if(root->right){
        tmp = root->right;
    }else if(root->left){
        tmp = root->left;
    }
    
    if(tmp){
        if(root->next)
        {
            if(root->next->left)
            {
                tmp->next = root->next->left;
            }
            else if(root->next->right)
            {
                tmp->next = root->next->right;
            }
            else
            {
                tmp->next = NULL;
            }
        }else
        {
            tmp->next = NULL;
        }
    }
    
    if(root->right)
    {
        connect(root->right);
    }
    if(root->left)
    {
        connect(root->left);
    }
}

int work(int left, int right)
{
    if(left >= right)
    {
        return 1;
    }
    
    int sum = 0;
    for(int i = left; i <= right; i++)
    {
        int leftcount = work(left, i-1);
        int rightcount = work(i+1, right);
        sum += (leftcount * rightcount);
    }
    
    return sum;
}

int numTrees(int n) {
    
    return work(1,n);
    
}

int reverse(int x) {
    int sign = 1;
    if(x <0) sign = -1;
    x = abs(x);
    int res = 0;
    while(x >0)
    {
        int tmp = x % 10;
        res = res *10 + tmp;
        x = x /10;
    }
    return res;
}

int firstMissingPositive(int A[], int n) {

    for(int i = 0;i< n;i++){
        while(A[i] >0 && A[i] <n && A[i] != A[A[i]-1] )
        {
            swap(A[i], A[A[i]-1]);
        }
    }
    for(int i = 0;i< n;i++)
    {
        if(A[i] != i+1)
            return i+1;
    }
    return n;
}

bool isPalindrome(string s, int start, int end){
    while(start <= end){
        if(s[start] == s[end])
        {
            start ++;
            end --;
        }
        else
        {
            return false;
        }
    }
    return true;
}

void worker(vector<string>& tracker, vector<vector<string>>& res, string s, int start){
    if(start == s.length())
    {
        res.push_back(tracker);
        return;
    }
    
    for(int i = start;i< s.length();i++)
    {
        if(isPalindrome(s, start, i))
        {
            tracker.push_back(s.substr(start,i-start +1));
            worker(tracker, res, s, i+1);
            tracker.pop_back();
        }
    }
}

vector<vector<string>> partition(string s) {
    vector<string> tracker;
    vector<vector<string>> res;
    if(s.length() == 0)
    {
        return res;
    }
    worker(tracker, res, s, 0);
    return res;
    
}

int maximalRectangle(vector<vector<char> > &matrix) {
    if (matrix.empty())  return 0;
    const int m = (int)matrix.size();
    const int n = (int)matrix[0].size();
    vector<int> H(n, 0);
    vector<int> L(n, 0);
    vector<int> R(n, n);
    int ret = 0;
    for (int i = 0; i < m; ++i) {
        int left = 0, right = n;
        // calculate L(i, j) from left to right
        for (int j = 0; j < n; ++j)
        {
            if (matrix[i][j] == '1')
            {
                ++H[j];
                L[j] = max(L[j], left);
            }
            else
            {
                left = j+1;
                H[j] = 0;
                L[j] = 0;
                R[j] = n;
            }
        }
        // calculate R(i, j) from right to left
        for (int j = n-1; j >= 0; --j)
        {
            if (matrix[i][j] == '1')
            {
                R[j] = min(R[j], right);
                ret = max(ret, H[j]*(R[j]-L[j]));
            }
            else
            {
                right = j;
            }
        }
    }
    return ret;
}

int lengthOfLongestSubstring(string s) {
    int last[256];
    fill(last, last+ 256, -1);
    int start = 0;
    int maxLength = 0;
    for(int i = 0;i< s.length(); i++)
    {
        if(last[s[i] -'a'] >= start)
        {
            maxLength = max(maxLength, i-start);
            start = last[s[i] - 'a']+1;
            
        }
        last[s[i] - 'a'] = i;
    }
    return max(maxLength, (int)s.length() - start);
}

vector<int> preorderTraversal(TreeNode *root) {
    vector<int> res;
    if(!root) return res;
    stack<TreeNode*> s;
    s.push(root);
    while(!s.empty())
    {
        TreeNode* tmp = s.top();
        s.pop();
        if(tmp->right) s.push(tmp->right);
        if(tmp->left) s.push(tmp->left);
        res.push_back(tmp->val);
            
    }
    return res;
}

int maxPoints(vector<Point> &points) {
    int maxCount = 0;
    
    for(int i = 0; i< points.size();i++)
    {
        for(int j = i+1;j< points.size();j++)
        {
            int count = 2;
            int xDist = points[j].x - points[i].x;
            int yDist = points[j].y - points[i].y;
            
            double slope;
            if(yDist != 0)
                slope = xDist/yDist;
            else
                slope = INT_MAX;
            
            for(int k = j+1;k< points.size();k++)
            {
                int xd = points[k].x - points[i].x;
                int yd = points[k].y - points[i].y;
                double s;
                if(yDist != 0)
                    s = xd/yd;
                else
                    s = INT_MAX;
                
                if(s == slope)
                {
                    count ++;
                }
            }
            maxCount = max(maxCount, count);
        }
    }
    return maxCount;
}

int minimumTotal(vector<vector<int> > &triangle) {
    if(triangle.size() == 0) return 0;
    vector<int> first;
    vector<int> second;
    first = triangle[0];
    for(int i = 1;i< triangle.size(); i ++){
        for(int j = 0;j< triangle[i].size();j++)
        {
            int left = j-1 < 0 ? INT_MAX : first[j-1];
            int right = j >= first.size()? INT_MAX: first[j];
            second.push_back(triangle[i][j] + min(left, right));
        }
        swap(first, second);
        second.clear();
    }
    
    int res = INT_MAX;
    for(int i = 0;i< first.size(); i ++)
    {
        res = min(res, first[i]);
    }
    return res;
}

vector<vector<int> > generate(int numRows) {
    vector<vector<int>> res;
    if(numRows == 0 )return res;
    vector<int> tmp(1,1);
    res.push_back(tmp);
    for(int i = 1;i< numRows;i++)
    {
        vector<int> cur;
        for(int j = 0;j< i+1;j ++)
        {
            int left = j-1 <0? 0 : res[i-1][j-1];
            int right = j >= res[i-1].size() ? 0 : res[i-1][j];
            cur.push_back(left+right);
        }
        res.push_back(cur);
    }
    return res;
}

vector<int> getRow(int rowIndex) {
    vector<int> res;
    res.push_back(1);
    if(rowIndex == 0) return res;
    int i = 1;
    vector<int> tmp;
    while(i <= rowIndex)
    {
        swap(tmp, res);
        res.clear();
        for (int j = 0;j< tmp.size()+1;j++)
        {
            int left = j-1<0? 0: tmp[j-1];
            int right = j>= tmp.size()? 0 : tmp[j];
            res.push_back(left+right);
        }
        i++;
    }
    return res;
}

TreeNode* worker(vector<int> &preorder, int preStart, int preEnd, vector<int>& inorder, int inStart, int inEnd)
{
    if(preEnd < preStart && inEnd < inStart) return nullptr;
    if(preStart == preEnd && inStart == inEnd) return new TreeNode(preorder[preStart]);
    int val = preorder[preStart];
    int index = 0;
    
    for(int i = inStart;i<=inEnd; i++)
    {
        if(inorder[i] == preorder[preStart])
        {
            index = i;
            break;
        }
    }
    
    int leftLength = index -inStart;
    TreeNode* root = new TreeNode(val);
    root->left = worker(preorder, preStart + 1, preStart + leftLength, inorder, inStart, inStart + leftLength-1);
    root->right = worker(preorder, preStart +leftLength +1, preEnd, inorder, inStart + leftLength+1, inEnd);
    return root;
}

TreeNode *buildTree(vector<int> &preorder, vector<int> &inorder)
{
    return worker(preorder, 0, (int)preorder.size()-1, inorder, 0, (int)inorder.size()-1);
}

int numTreesII(int n)
{
    int tracker[n+1];
    tracker[0] = 0;
    tracker[1] = 1;
    for(int i = 2;i < n;i++)
    {
        for(int j = 1;j<=2;j++)
        {
            tracker[i] += tracker[j-1] * tracker[i-j];
        }
    }
    return tracker[n];
}

vector<TreeNode*> workerMe(int left, int right)
{
    vector<TreeNode*> res;
    if(left > right)
    {
        res.push_back(nullptr);
        return res;
    }
    
    if(left == right)
    {
        res.push_back(new TreeNode(left));
        return res;
    }
    
    for(int i= left; i<=right;i++)
    {
        vector<TreeNode*> leftTree = workerMe(left, i-1);
        vector<TreeNode*> rightTree = workerMe(i+1, right);
        
        for(TreeNode* leftNode: leftTree){
            for(TreeNode* rightNode: rightTree){
                TreeNode* root = new TreeNode(i);
                root->left = leftNode;
                root->right = rightNode;
                res.push_back(root);
            }
        }
    }
    return res;
}

vector<TreeNode *> generateTrees(int n) {
    return workerMe(1, n);
}

bool IsValidBST(TreeNode* root, int left, int right)
{
    if(!root) return true;
    if(root->val >=left && root->val <= right){
        if(IsValidBST(root->left, left, root->val) && IsValidBST(root->right, root->val, right))
        {
            return true;
        }
    }
    return false;
}

bool isValidBST(TreeNode *root) {
    return IsValidBST(root, INT_MIN, INT_MAX);
}


int main(int argc, const char * argv[])
{
    /*
    int A[] = {100};
    int B[] = {101};
    
    double res = findMedianSortedArrays(A, 1, B, 1);
    cout<<res;
    */
    //vector<int> input = {3,2,4};
    //vector<int> res = twoSum(input, 6);
    
    //string s = "wlrbbmqbhcdarzowkkyhiddqscdxrjmowfrxsjybldbefsarcbynecdyggxxpklorellnmpapqfwkhopkmco";
    //int res = lengthOfLongestSubstring(s);
    //cout<<res;
    
    //ListNode* l1 = new ListNode(0);
    //ListNode* l2 = new ListNode(1);
    
    //ListNode* res = addTwoNumbers(l1, l2);
    
    //cout<<res->val <<endl;
    
    //string input = "abb";
    //string res = longestPalindrome(input);
    //cout<<res<<endl;
    
    /*TreeNode* n1 = new TreeNode(3);
    TreeNode* n2 = new TreeNode(1);
    TreeNode* n3 = new TreeNode(4);
    TreeNode* n4 = new TreeNode(2);
    n1->left = n2;
    n2->left = n3;
    n3->right = n4;
    flatten(n1);
    */
    
    //string res = "   one.   +two three?   ~four   !five- ";
    //reverseWords(res);
    //cout<<res.length()<<endl;
    //cout <<res<<endl;
    //vector<string> input = {"10","6","9","3","+","-11","*","/","*","17","+","5","+"};
    //int res = evalRPN(input);
    
    //TreeLinkNode* root = new TreeLinkNode(0);
    //connect(root);
    //int res = numTrees(2);
    //cout<<res<<endl;
    //vector<vector<string>> res =partition("a");
    //cout<<res.size() <<endl;
    /*vector<int> test (5,5);
    for(int i = 0;i< test.size();i++)
    {
        cout<<test[i]<<endl;
    }*/
    /*TreeNode* n1 = new TreeNode(1);
     TreeNode* n2 = new TreeNode(2);
     TreeNode* n3 = new TreeNode(3);
     TreeNode* n4 = new TreeNode(4);
     n1->left = n2;
     n1->right= n3;
     n2->right = n4;
    vector<int> res = preorderTraversal(n1);
    */
    vector<int> pre = {1,2,3};
    vector<int> ino = {1,2,3};
    
    TreeNode* res = buildTree(pre, ino);
    
    return 0;
}

