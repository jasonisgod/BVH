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
		
		public Tree(int dim_ = 3, double gap_ = 0.05)
		{
			dim = dim_;
			gap = gap_;
			root = null;
		}

		public List<Node> Query(AABB aabb, bool touch = false)
		{
			Stack<Node> stack = new Stack<Node>();
			stack.Push(root);

			List<Node> ans = new List<Node>();

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
			Node newParent = new Node(leaf.aabb | sibling.aabb, oldParent);

			if (oldParent is not null)
			{
				if (oldParent.left == sibling)
				{
					oldParent.left = newParent;
				}
				else
				{
					oldParent.right = newParent;
				}

				newParent.left = sibling;
				newParent.right = leaf;
				sibling.parent = newParent;
				leaf.parent = newParent;
			}
			else
			{
				newParent.left = sibling;
				newParent.right = leaf;
				sibling.parent = newParent;
				leaf.parent = newParent;
				root = newParent;
			}
			Update(newParent);
			
			node = leaf.parent;
			while (node is not null)
			{
				Balance(node);
				node = node.parent;
			}
			
		}

		public void Remove(Node leaf)
		{
			if (leaf == root)
			{
				root = null;
				return;
			}

			Node parent = leaf.parent;
			Node grandParent = parent.parent;
			Node sibling;

			if (parent.left == leaf)
			{
				sibling = parent.right;
			}
			else
			{
				sibling = parent.left;
			}

			if (grandParent is not null)
			{
				if (grandParent.left == parent)
				{
					grandParent.left = sibling;
				}
				else
				{
					grandParent.right = sibling;
				}

				sibling.parent = grandParent;

				Node index = grandParent;
				while (index is not null)
				{
					Balance(index);
					index = index.parent;
				}
			}
			else
			{
				root = sibling;
				sibling.parent = null;
			}
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
				// node.aabb.fatten(gap);
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
			
			// L-RL-RR
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

		public void Draw(String file, double size, int scale)
		{
			int width = Convert.ToInt32(size * scale);
			int height = Convert.ToInt32(size * scale);
			Bitmap bitmap = new Bitmap(width, height); //, System.Drawing.Imaging.PixelFormat.Format32bppArgb
			Graphics g = Graphics.FromImage(bitmap);
			g.Clear(Color.White);
			_Draw(root, 0);
			bitmap.Save(file, System.Drawing.Imaging.ImageFormat.Png);
			
			void _Draw(Node node, int depth)
			{
				if (node is null) 
				{
					return;
				}

				AABB aabb = node.aabb;
				int x = Convert.ToInt32(scale * aabb.lb[0]);
				int y = Convert.ToInt32(scale * aabb.lb[1]);
				int w = Convert.ToInt32(scale * (aabb.ub[0] - aabb.lb[0]));
				int h = Convert.ToInt32(scale * (aabb.ub[1] - aabb.lb[1]));

				if (node.isLeaf())
				{
					Pen pen = new Pen(Color.Blue, 1);
					g.DrawEllipse(pen, x, y, w, h);
				} else {
					Pen pen = new Pen(Color.Black, 1);
					g.DrawRectangle(pen, x, y, w, h);
					// int rgb = Math.Min(255, 128 + depth * 20);
					// SolidBrush brush = new SolidBrush(Color.FromArgb(rgb, rgb, rgb));
					//g.FillRectangle(brush, x, y, w, h);
				}

				_Draw(node.right, depth + 1);
				_Draw(node.left, depth + 1);
			}
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
