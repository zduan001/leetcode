package com.test.leetcode.tst;

import static org.junit.Assert.*;

import org.junit.Test;

import com.test.leetcode.LeetCode;
import com.test.leetcode.util.LinkedListNode;

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

}
