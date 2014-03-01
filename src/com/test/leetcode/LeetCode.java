package com.test.leetcode;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Stack;

import javax.swing.text.html.HTMLDocument.Iterator;

import com.test.leetcode.util.LinkedListNode;
import com.test.leetcode.util.ListNode;
import com.test.leetcode.util.Point;
import com.test.leetcode.util.TreeNode;

public class LeetCode {

	public LinkedListNode<Integer> _001_reverseLinkList(
			LinkedListNode<Integer> input) {
		if (input == null) {
			return null;
		}

		LinkedListNode<Integer> newHead = null;
		LinkedListNode<Integer> temp;

		while (input != null) {
			temp = input;
			input = input.next;
			temp.next = newHead;
			newHead = temp;
		}

		return newHead;

	}

	/**
	 * http://www.mitbbs.com/article_t/JobHunting/32633319.html
	 * 有一种String,是把一个更短的String重复n次而构成的，那个更短的String长度至少为 2，输入一个String写代码返回T或者F
	 * 例子： "abcabcabc" Ture 因为它把abc重复3次构成 "bcdbcdbcde" False 最后一个是bcde
	 * "abcdabcd" True 因为它是abcd重复2次构成 "xyz" False 因为它不是某一个String重复 "aaaaaaaaaa"
	 * False 重复的短String长度应至少为2（这里不能看做aa重复5次)
	 * 
	 * 要求算法复杂度为O(n)
	 * 
	 * @param input
	 * @return
	 */
	public boolean _002_CanbeRepeat(String input) {
		if (input == null || input.length() < 4) {
			return false;
		}

		StringBuilder sb = new StringBuilder();
		int tracker[] = new int[input.length()];

		int indexInSb = 0;
		sb.append(input.charAt(0));
		tracker[0] = 1;
		for (int i = 1; i < input.length(); i++) {
			if (input.charAt(i) == sb.charAt(indexInSb)) {
				tracker[indexInSb]++;
				indexInSb = (indexInSb + 1) % sb.length();
			} else {
				String tmp;
				if (indexInSb > 0) {
					tmp = sb.substring(0, indexInSb - 1);
				} else {
					tmp = sb.toString();
				}
				for (int j = 0; j < tracker[0] - 1; j++) {
					sb.append(tmp);
				}
				sb.append(input.charAt(i));
				for (int k = 0; k < sb.length(); k++) {
					tracker[k] = 1;
				}
				indexInSb = 0;
			}
		}
		if (sb.length() < 2 || sb.length() == input.length()) {
			return false;
		}
		int first = tracker[0];
		for (int i = 1; i < tracker.length; i++) {
			if (tracker[i] == 0) {
				break;
			}
			if (tracker[i] != first) {
				return false;
			}
		}
		return true;
	}

	/**
	 * http://oj.leetcode.com/problems/max-points-on-a-line/
	 * 
	 * @param points
	 * @return
	 */
	public int _003_maxPoints(Point[] points) {
		if (points.length == 1) {
			return 1;
		}

		int count = 0;
		HashMap<Double, Integer> tracker = new HashMap<Double, Integer>();
		int identicalcount = 0;

		for (int i = 0; i < points.length; i++) {
			tracker.clear();
			identicalcount = 1;
			for (int k = i + 1; k < points.length; k++) {

				double dy = points[k].y - points[i].y;
				double dx = points[k].x - points[i].x;
				if (dy == 0 && dx == 0) {
					identicalcount++;
					continue;
				}

				if (dy == 0) {
					if (tracker.containsKey(Double.MAX_VALUE)) {
						tracker.put(Double.MAX_VALUE,
								tracker.get(Double.MAX_VALUE) + 1);
					} else {
						tracker.put(Double.MAX_VALUE, 1);
					}
				} else {
					Double slope = dx / dy;
					if (slope == -0.0d) {
						slope = 0.0d;
					}
					if (tracker.containsKey(slope)) {
						tracker.put(slope, tracker.get(slope) + 1);
					} else {
						tracker.put(slope, 1);
					}
				}
			}
			java.util.Iterator<Entry<Double, Integer>> it = tracker.entrySet()
					.iterator();

			if (!it.hasNext()) {
				count = count > identicalcount ? count : identicalcount;
			} else {

				while (it.hasNext()) {
					Entry<Double, Integer> e = it.next();
					count = count > e.getValue() + identicalcount ? count : e
							.getValue() + identicalcount;
				}
			}
		}

		return count;
	}

	/**
	 * http://oj.leetcode.com/problems/container-with-most-water/ Given n
	 * non-negative integers a1, a2, ..., an, where each represents a point at
	 * coordinate (i, ai). n vertical lines are drawn such that the two
	 * endpoints of line i is at (i, ai) and (i, 0). Find two lines, which
	 * together with x-axis forms a container, such that the container contains
	 * the most water. Note: You may not slant the container.
	 * 
	 * @param height
	 * @return
	 */
	public int _004_maxArea(int[] height) {
		if (height == null || height.length < 2) {
			return 0;
		}

		int max = 0;
		int left = 0, right = height.length - 1;
		int tmp;

		while (left < right) {
			tmp = Math.min(height[left], height[right]) * (right - left);
			max = max > tmp ? max : tmp;
			if(height[left] >= height[right])
			{
				right --;
			}
			else
			{
				left ++;
			}
		}
		return max;
	}
	
	/**
	 * http://oj.leetcode.com/problems/trapping-rain-water/
	 *  Given n non-negative integers representing an elevation map where the width of each bar is 1, 
	 *  compute how much water it is able to trap after raining.
	 *  For example,
	 *  Given [0,1,0,2,1,0,1,3,2,1,2,1], return 6. 
	 *  {6,8,5,0,0,6,5} return 13
	 */
	public int _005_rainWater(int[] A) {
		if (A == null || A.length < 3) {
			return 0;
		}
		int max = 0;
		int left = 0, right = A.length - 1;
		int leftMax = A[left], rightMax = A[right];
		
		while(left < right){
			if(A[left] <= A[right]){
				if(leftMax > A[left+1]){
					max += leftMax - A[left+1];
				} else{
					leftMax = A[left + 1];
				}
				left ++;
			}
			else{
				if(rightMax > A[right-1]){
					max += rightMax - A[right-1];
				} else {
					rightMax = A[right-1];
				}
				right --;
			}
		}
		return max;
	}
	
	/**
	 * http://oj.leetcode.com/problems/search-in-rotated-sorted-array/
	 * Suppose a sorted array is rotated at some pivot unknown to you beforehand.
	 * (i.e., 0 1 2 4 5 6 7 might become 4 5 6 7 0 1 2).
	 * You are given a target value to search. If found in the array return its index, otherwise return -1.
	 * You may assume no duplicate exists in the array.
	 * @param A
	 * @param target
	 * @return
	 */
	public int _006_search(int[] A, int target) {
		if(A == null || A.length == 0){
			return -1;
		}
		return _006_searchRotateArray(A,  target, 0, A.length-1);
	}
	
	private int _006_searchRotateArray(int[] A, int target, int left, int right)
	{
		if(left > right){
			return -1;
		}
		
		if(A[left] == target){
			return left;
		}else if(A[right] == target){
			return right;
		}
		
		int mid = left + (right - left) / 2;
		if(A[mid] == target){
			return mid;
		}
		
		if(A[0] < A[mid]){
			if(A[0] <= target && target < A[mid]){
				return _006_searchRotateArray(A, target, left, mid-1);
			}else {
				return _006_searchRotateArray(A, target, mid + 1, right);
			}
		} else {
			if(A[mid] < target && target <= A[right]){
				return _006_searchRotateArray(A, target, mid +1, right);
			} else {
				return _006_searchRotateArray(A, target, left, mid -1);
			}
		}
	}

	/**
	 *  Evaluate the value of an arithmetic expression in Reverse Polish Notation.
	 *  Valid operators are +, -, *, /. Each operand may be an integer or another expression.
	 *  Some examples:
	 *  ["2", "1", "+", "3", "*"] -> ((2 + 1) * 3) -> 9
	 *  ["4", "13", "5", "/", "+"] -> (4 + (13 / 5)) -> 6
	 * @param tokens
	 * @return
	 */
	public int _007_evalPRN(String[] tokens){
		if(tokens == null || tokens.length == 0){
			return 0;
		}
			
		Integer first, second;
		Stack<Integer> numbers = new Stack<Integer>();
		for(int i = 0;i< tokens.length;i ++){
			
			if(tokens[i].equals("+")){
				second = numbers.pop();
				first = numbers.pop();
				numbers.push(first + second);
			}else if(tokens[i].equals("-")){
				second = numbers.pop();
				first = numbers.pop();
				numbers.push(first - second);
			}else if(tokens[i].equals("*")){
				second = numbers.pop();
				first = numbers.pop();
				numbers.push(first * second);
			}else if(tokens[i].equals("/")){
				second = numbers.pop();
				first = numbers.pop();
				numbers.push(first / second);
			} else {
				numbers.push(Integer.parseInt(tokens[i]));
			}
		}
		return numbers.pop();
	}
	
	/**
	 * http://oj.leetcode.com/problems/sort-list/
	 * Sort a linked list in O(n log n) time using constant space complexity.
	 * @param head
	 * @return
	 */
	public ListNode _008_sortList(ListNode head){
		return _008_sortList_worker( head);
	}
	
	private ListNode _008_sortList_worker(ListNode head){
		if(head == null || head.next == null) {
			return head;
		}
		
		if(head.next.next == null){
			ListNode tail = null;
			if(head.val > head.next.val){
				tail = head;
				head = head.next;
				tail.next = null;
				head.next = tail;
			}
			return head;
		}
		
		ListNode fast = head, slow = head;
		while(fast != null){
			if(fast.next == null){
				fast = null;
				break;
			}else{
				fast = fast.next.next;
			}
			slow = slow.next;
		}
		
		ListNode secondHead = slow.next;
		slow.next = null;
		head = _008_sortList_worker(head);
		secondHead = _008_sortList_worker(secondHead);
		return _008_mergeTwoList(head, secondHead);
	}
	
	private ListNode _008_mergeTwoList(ListNode head1, ListNode head2){
		ListNode dummyHead = new ListNode(-1);
		ListNode tmp;
		ListNode tail = dummyHead;
		
		while(head1 != null || head2 != null){
			if(head1 == null){
				tail.next = head2;
				break;
			}
			if(head2 == null){
				tail.next = head1;
				break;
			}
			
			if(head1.val <= head2.val){
				tmp = head1;
				head1 = head1.next;
			} else {
				tmp = head2;
				head2 = head2.next;
			}
			
			tail.next = tmp;
			tail = tail.next;
		}
		
		return dummyHead.next;
		
	}
	
	/**
	 *  Given a singly linked list L: L0→L1→…→Ln-1→Ln,
	 *  reorder it to: L0→Ln→L1→Ln-1→L2→Ln-2→…
	 *  
	 *  You must do this in-place without altering the nodes' values.
	 *  
	 *  For example,
	 *  Given {1,2,3,4}, reorder it to {1,4,2,3}. 
	 * @param head
	 * @return
	 */
	public ListNode _009_reorderList(ListNode head){
		if(head == null || head.next == null || head.next.next == null){
			return head;
		}
		
		ListNode fast = head.next, slow = head;
		while(fast != null){
			if(fast.next == null){
				fast = null;
				break;
			}else{
				fast = fast.next.next;
			}
			slow = slow.next;
		}
		
		ListNode secondHead = slow.next;
		slow.next = null;
		
		ListNode reversedHead = null;
		ListNode tmp = secondHead;
		
		while(secondHead != null){
			tmp = secondHead;
			secondHead = secondHead.next;
			tmp.next = reversedHead;
			reversedHead = tmp;
		}
		
		tmp = head;
		ListNode nextHead;
		while(reversedHead != null && tmp != null){
			nextHead = tmp.next;
			tmp.next = reversedHead;
			reversedHead = reversedHead.next;
			tmp.next.next = nextHead;
			tmp = nextHead;
		}
		
		return head;
	}

	/**
	 * http://oj.leetcode.com/problems/insertion-sort-list/
	 * Sort a linked list using insertion sort.
	 * @param head
	 * @return
	 */
	public ListNode _010_insertionSortList(ListNode head){
		if(head == null || head.next == null){
			return head;
		}
		
		ListNode dummyHead = new ListNode(-1);
		ListNode dHead1 = new ListNode(-1);
		dHead1.next = head;
		
		ListNode tmp;
		ListNode count = dummyHead;
		
		while(dHead1.next != null){
			tmp = dHead1.next;
			dHead1.next = tmp.next;
			count = dummyHead;
			while(count.next != null && count.next.val <= tmp.val){
				count = count.next;
			}
			tmp.next = count.next;
			count.next = tmp;
		}
		return dummyHead.next;
	}

    public ArrayList<Integer> _011_postorderTraversal(TreeNode root) {
        ArrayList<Integer> res= new ArrayList<Integer>();
        
        Stack<TreeNode> tracker = new Stack<TreeNode>();
        TreeNode previous = null;
        
        tracker.push(root);
        TreeNode tmp;
        while(!tracker.empty()){
        	tmp = tracker.peek();
        	if(tmp.left != null && ( previous == null || previous != tmp.right)){
        		tracker.push(tmp.left);
        		previous = null;
        		continue;
        	}
        	
        	if(tmp.right != null && previous ==null){
        		tracker.push(tmp.right);
        		continue;
        	}
        	
        	if(tmp.right == null || tmp.right == previous){
        		res.add(tmp.val);
        		tracker.pop();
        		previous = tmp;
        	}
        }
        
        return res;
    }
}
