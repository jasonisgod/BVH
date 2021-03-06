using System; 
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
//using System.Drawing.Common;
using System.Linq;

namespace BVH
{
	public class Tree
	{
		public Node root;
		private int dim;
		private double gap;
		
		public Tree(int dim_ = 2, double gap_ = 2)
		{
			root = null;
			dim = dim_;
			gap = gap_;
		}

		public List<Node> Query(AABB aabb, bool touch = false)
		{
			List<Node> ans = new List<Node>();
			Stack<Node> stack = new Stack<Node>();
			stack.Push(root);

			while (stack.Count > 0)
			{
				Node node = stack.Pop();
				if (node is null)
				{
					continue;
				}

				if (node.aabb.Overlaps(aabb, touch))
				{
					if (node.isLeaf())
					{
						ans.Add(node);
					}
					else
					{
						stack.Push(node.left);
						stack.Push(node.right);
					}
				}
			}

			return ans;
		}
		
		public void Insert(Node leaf)
		{
			if (leaf is null) 
			{
				return;
			}

			if (root is null)
			{
				root = leaf;
				root.parent = null;
				return;
			}
			
			Node node = root;

			while (!node.isLeaf())
			{
				Node left = node.left;
				Node right = node.right;
				double cost = 2 * node.aabb.Area();
				double costLeft = (leaf.aabb | left.aabb).Area() - (left.isLeaf()? 0: left.aabb.Area());
				double costRight = (leaf.aabb | right.aabb).Area() - (right.isLeaf()? 0: right.aabb.Area());

				if (cost < costLeft && cost < costRight)
				{
					break;
				}

				node = (costLeft < costRight? left: right);
			}

			Node sibling = node;
			Node oldParent = sibling.parent;
			Node parent = new Node(leaf.aabb | sibling.aabb, oldParent);
			parent.left = sibling;
			parent.right = leaf;
			sibling.parent = parent;
			leaf.parent = parent;

			if (oldParent is not null)
			{
				(oldParent.left == sibling? ref oldParent.left: ref oldParent.right) = parent;
			}
			else
			{
				root = parent;
			}
			Update(parent);
			Balance(parent);
		}

		public void Remove(Node leaf)
		{
			if (leaf is null)
			{
				return;
			}

			if (leaf == root)
			{
				root = null;
				return;
			}

			Node parent = leaf.parent;
			Node grandParent = parent.parent;
			Node sibling = (parent.left == leaf? parent.right: parent.left);

			if (grandParent is not null)
			{
				(grandParent.left == parent? ref grandParent.left: ref grandParent.right) = sibling;
				sibling.parent = grandParent;
			}
			else
			{
				root = sibling;
				sibling.parent = null;
			}
			Update(grandParent);
			Balance(grandParent);
		}

		private void Update(Node node)
		{
			if (node is null) 
			{
				return;
			}
			if (!node.isLeaf())
			{
				Node left = node.left;
				Node right = node.right;
				node.height = 1 + Math.Max(left.height, right.height);
				node.aabb = (left.aabb | right.aabb);
			}
			Update(node.parent);
		}

		private void Swap(Node node1, Node node2)
		{
			if (node1 is null || node2 is null) 
			{
				return;
			}

			Node par1 = node1.parent;
			Node par2 = node2.parent;
			if (par1 is null || par2 is null)
			{
				return;
			}

			node1.parent = par2;
			node2.parent = par1;
			(par1.left == node1? ref par1.left: ref par1.right) = node2;
			(par2.left == node2? ref par2.left: ref par2.right) = node1;

			Update(node1);
			Update(node2);
		}

		private void Balance(Node node)
		{
			if (node is null)
			{
				return;
			}

			if (node.height < 2)
			{
				Balance(node.parent);
				return;
			}

			var isBalanced = (int height1, int height2) => Math.Abs(height1 - height2) < 2;
			var group = (int height1, int height2) => Math.Max(height1, height2) + 1;
			
			// L-RL-RR & R-LL-LR
			{
				Node nL = node.left;
				Node nR = node.right;
				Node nLL = nL.left;
				Node nLR = nL.right;
				Node nRL = nR.left;
				Node nRR = nR.right;

				var list = new[]
				{
					new {index = 0, nL = nR , nRL = nLL, nRR = nLR},
					new {index = 1, nL = nLL, nRL = nR , nRR = nLR},
					new {index = 2, nL = nLR, nRL = nR , nRR = nLL},
					new {index = 0, nL = nL , nRL = nRL, nRR = nRR},
					new {index = 3, nL = nRL, nRL = nL , nRR = nRR},
					new {index = 4, nL = nRR, nRL = nL , nRR = nRL},
				}.ToList()
					.Where(_ => (_.nL is not null && _.nRL is not null && _.nRR is not null))
					.Where(_ => isBalanced(_.nRL.height, _.nRR.height))
					.Where(_ => isBalanced(group(_.nRL.height, _.nRR.height), _.nL.height))
					.Select(_ => new {aabb1 = (_.nRL.aabb | _.nRR.aabb), aabb2 = _.nL.aabb, index = _.index})
					.OrderBy(_ => (_.aabb1 & _.aabb2).Area())
					// .ThenBy(_ => _.aabb1.Area() + _.aabb2.Area());
					.ThenBy(_ => Math.Max(_.aabb1.Area(), _.aabb2.Area()));
					// .ThenBy(_ => Math.Abs(_.aabb1.Area() - _.aabb2.Area()))

				var index = list.Any()? list.First().index: 0;
				switch (index)
				{
					case 1: Swap(nR, nLL); break;
					case 2: Swap(nR, nLR); break;
					case 3: Swap(nL, nRL); break;
					case 4: Swap(nL, nRR); break;
					default: break;
				}
				// Console.WriteLine("3 nodes : " + index.ToString());
			}

			// LL-LR-RL-RR
			{
				Node nL = node.left;
				Node nR = node.right;
				Node nLL = nL.left;
				Node nLR = nL.right;
				Node nRL = nR.left;
				Node nRR = nR.right;

				var list = new[]
				{
					new {index = 0, nLL = nLL, nLR = nLR, nRL = nRL, nRR = nRR},
					new {index = 1, nLL = nLL, nLR = nRL, nRL = nLR, nRR = nRR},
					new {index = 2, nLL = nLL, nLR = nRR, nRL = nRL, nRR = nLR},
				}.ToList()
					.Where(_ => (_.nLL is not null && _.nLR is not null && _.nRL is not null && _.nRR is not null))
					.Where(_ => isBalanced(_.nLL.height, _.nLR.height))
					.Where(_ => isBalanced(_.nRL.height, _.nRR.height))
					.Where(_ => isBalanced(group(_.nLL.height, _.nLR.height), group(_.nRL.height, _.nRR.height)))
					.Select(_ => new {aabb1 = (_.nLL.aabb | _.nLR.aabb), aabb2 = (_.nRL.aabb | _.nRR.aabb), index = _.index})
					.OrderBy(_ => (_.aabb1 & _.aabb2).Area())
					.ThenBy(_ => _.aabb1.Area() + _.aabb2.Area());
					// .ThenBy(_ => Math.Max(_.aabb1.Area(), _.aabb2.Area()))
					// .ThenBy(_ => Math.Abs(_.aabb1.Area() - _.aabb2.Area()))
				
				var index = list.Any()? list.First().index: 0;
				switch (index)
				{
					case 1: Swap(nLR, nRL); break;
					case 2: Swap(nLR, nRR); break;
					default: break;
				}
				// Console.WriteLine("4 nodes : " + index.ToString());
			}

			Balance(node.parent);
		}

		public override string ToString()
		{
			string s = "";
			// foreach (Node node in nodes)
			// {
			// 	s += node.aabb + "\n";
			// }
			return s;
		}
	}
}
