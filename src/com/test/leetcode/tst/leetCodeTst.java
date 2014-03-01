package com.test.leetcode.tst;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

import com.test.leetcode.LeetCode;
import com.test.leetcode.util.LinkedListNode;
import com.test.leetcode.util.ListNode;
import com.test.leetcode.util.Point;
import com.test.leetcode.util.TreeNode;

public class leetCodeTst {

	@Test
	public void test_001_reverseLinkedList() {
		LinkedListNode<Integer> input = new LinkedListNode<Integer>();
		input.value = 1;
		LinkedListNode<Integer> input2 = new LinkedListNode<Integer>();
		input2.value = 2;
		LinkedListNode<Integer> input3 = new LinkedListNode<Integer>();
		input3.value = 3;
		LinkedListNode<Integer> input4 = new LinkedListNode<Integer>();
		input4.value = 4;
		LinkedListNode<Integer> input5 = new LinkedListNode<Integer>();
		input5.value = 5;
		LinkedListNode<Integer> input6 = new LinkedListNode<Integer>();
		input6.value = 6;
		LinkedListNode<Integer> input7 = new LinkedListNode<Integer>();
		input7.value = 7;
		
		input.next = input2;
		input2.next = input3;
		input3.next = input4;
		input4.next = input5;
		input5.next = input6;
		input6.next = input7;
		
		LeetCode code = new LeetCode();
		
		LinkedListNode res = code._001_reverseLinkList(input);
		
		assertSame(res, input7);
		
		for(Integer i = 7;i>0;i--){
			assertEquals(i, res.value);
			res = res.next;
		}
	}
	
	@Test
	public void test_001_1_reverseLinkedList(){
		LinkedListNode<Integer> input = new LinkedListNode<Integer>();
		input.value = 1;
		
		LeetCode code = new LeetCode();
		
		LinkedListNode res = code._001_reverseLinkList(input);
		
		assertEquals(1,res.value);
		assertNull(res.next);
		assertSame(res, input);
	}
	
	@Test
	public void test_002_CanbeRepeat(){
		LeetCode code = new LeetCode();
		assertTrue(code._002_CanbeRepeat("abcabcabc"));
		assertFalse(code._002_CanbeRepeat("bcdbcdbcde"));
		assertTrue(code._002_CanbeRepeat("abcdabcd"));
		assertFalse(code._002_CanbeRepeat("aaaaaaaaaa"));
		assertFalse(code._002_CanbeRepeat("xyz"));
		//assertTrue(code._002_CanbeRepeat("abaaaabaaa"));
	}

	@Test
	public void test_003_maxPoints(){
		LeetCode code = new LeetCode();
		
		Point[] points = new Point[]{new Point(0,0), new Point(0,0)};
		
		int res;
		res = code._003_maxPoints(points);
		assertEquals(2, res);
		
		points = new Point[]{new Point(0,0), new Point(0,0), new Point(0,1)};
		res = code._003_maxPoints(points);
		assertEquals(3, res);
		
		points = new Point[]{new Point(0,0), new Point(1,1), new Point(1,-1)};
		res = code._003_maxPoints(points);
		assertEquals(2, res);
		
		points = new Point[]{new Point(4,0), new Point(4,-1), new Point(4,5)};
		res = code._003_maxPoints(points);
		assertEquals(3, res);
	}
	
	@Test
	public void test_004_maxArea()
	{
		LeetCode code = new LeetCode();
		
		int heights[] = {10,9,8,7,6,5,4,3,2,1};
		
		int res;
		res = code._004_maxArea(heights);
		
		assertEquals(25, res);
	}
	
	@Test
	public void test_005_rainWater()
	{
		LeetCode code = new LeetCode();
		int[] A = {0,1,0,2,1,0,1,3,2,1,2,1};
		
		int res;
		res = code._005_rainWater(A);
		
		assertEquals(6, res);
		
		A = new int[]{6,8,5,0,0,6,5};
		res = code._005_rainWater(A);
		assertEquals(13,res);
	}
	
	@Test
	public void test_006_searchRotateArray()
	{
		LeetCode code = new LeetCode();
		int A[] = {4, 5, 6, 7, 0, 1, 2};
		
		int res;
		res = code._006_search(A, 1);
		assertEquals(5, res);
		
		res = code._006_search(A, 6);
		assertEquals(2, res);
		
		A = new int[] {1,3};
		res = code._006_search(A, 3);
		assertEquals(1, res);
		
		A = new int[] {3, 1};
		res = code._006_search(A, 1);
		assertEquals(1, res);
	}
	
	@Test
	public void test_009_reOrderList()
	{
		LeetCode code = new LeetCode();
		ListNode n1 = new ListNode(1);
		ListNode n2 = new ListNode(2);
		ListNode n3 = new ListNode(3);
		ListNode n4 = new ListNode(4);
		n1.next = n2;
		n2.next = n3;
		n3.next = n4;
		
		ListNode res;
		res = code._009_reorderList(n1);
		
		assertEquals(1, n1.val);
		assertEquals(4, n1.next.val);
		assertEquals(2, n1.next.next.val);
		assertEquals(3, n1.next.next.next.val);
		
	}
	
	@Test
	public void test_010_insertionSortList()
	{
		LeetCode code = new LeetCode();
		ListNode n1 = new ListNode(1);
		ListNode n2 = new ListNode(1);
		n1.next = n2;
		
		ListNode res;
		res = code._010_insertionSortList(n1);
	}
	
	@Test
	public void test_011_postorderTreeTravel(){
		LeetCode code = new LeetCode();
		TreeNode t1 = new TreeNode(1);
		TreeNode t2 = new TreeNode(2);
		TreeNode t3 = new TreeNode(3);
		t1.right = t2;
		t2.left = t3;
		
		ArrayList<Integer> res = code._011_postorderTraversal(t1);
		assertEquals(new Integer(3), res.get(0));
		
	}
	

}
