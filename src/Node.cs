using System; 
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Threading;

namespace BVH
{
	public class Node
	{
		public AABB aabb;
		public Node parent;
		public Node left;
		public Node right;
		public int height;
		
		public Node(AABB aabb_, Node parent_)
		{
			aabb = aabb_;
			parent = parent_;
			left = null;
			right = null;
			height = parent_ is null? 0: parent_.height + 1;
		}
		public bool isLeaf()
		{
			return (left is null && right is null);
		}
	}

}
