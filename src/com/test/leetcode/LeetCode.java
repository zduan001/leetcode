package com.test.leetcode;


import com.test.leetcode.util.LinkedListNode;

public class LeetCode {
	
	public LinkedListNode<Integer> _001_reverseLinkList(LinkedListNode<Integer> input)
	{
		if(input == null){
			return null;
		}
		
		LinkedListNode<Integer> newHead = null;
		LinkedListNode<Integer> temp;
		
		while(input != null){
			temp = input;
			input = input.next;
			temp.next = newHead;
			newHead = temp;
		}
		
		
		return newHead;
			
	}

}
