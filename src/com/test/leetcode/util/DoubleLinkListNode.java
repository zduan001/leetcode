package com.test.leetcode.util;

public class DoubleLinkListNode {

	public int val;
	public DoubleLinkListNode next;
	public DoubleLinkListNode previous;

	public DoubleLinkListNode(int x) {
		val = x;
		previous = null;
		next = null;
	}

}
