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

				if (node == null)
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
			if (root == null)
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

			if (oldParent != null)
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
			while (node != null)
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

			if (grandParent != null)
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
				while (index != null)
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
			if (node == null) 
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
			if (node1 == null || node2 == null) 
			{
				return;
			}

			Node par1 = node1.parent;
			Node par2 = node2.parent;
			if (par1 == null || par2 == null)
			{
				return;
			}

			node1.parent = par2;
			node2.parent = par1;
			if (par1.left == node1) par1.left = node2;
				else par1.right = node2;
			if (par2.left == node2) par2.left = node1;
				else par2.right = node1;

			Update(node1);
			Update(node2);
		}

		private void Balance(Node node)
		{
			if (node == null)
			{
				return;
			}
			
			if (node.height < 2)
			{
				Balance(node.parent);
				return;
			}

			bool b(Node node1, Node node2)
			{
				return Math.Abs(node1.height - node2.height) < 2;
			}

			Node u(Node node1, Node node2)
			{
				return new Node() {height = Math.Max(node1.height, node2.height) + 1};
			};
			
			// Console.WriteLine("balance()");
			// 3-nodes
			{
				Node? nL = node.left;
				Node? nR = node.right;
				Node? n1 = nL.left;
				Node? n2 = nL.right;
				Node? n3 = nR.left;
				Node? n4 = nR.right;

				var list = new[]
				{
					new {index = 0, n0 = nR, n1 = n1, n2 = n2},
					new {index = 1, n0 = n1, n1 = nR, n2 = n2},
					new {index = 2, n0 = n2, n1 = nR, n2 = n1},
					new {index = 0, n0 = nL, n1 = n3, n2 = n4},
					new {index = 3, n0 = n3, n1 = nL, n2 = n4},
					new {index = 4, n0 = n4, n1 = nL, n2 = n3},
				}.ToList()
					.Where(_ => (_.n0 != null && _.n1 != null && _.n2 != null))
					.Where(_ => b(_.n1, _.n2) && b(u(_.n1, _.n2), _.n0))
					.Select(_ => new {aabb1 = (_.n1.aabb | _.n2.aabb), aabb2 = _.n0.aabb, index = _.index})
					.OrderBy(_ => (_.aabb1 & _.aabb2).Area())
					// .ThenBy(_ => _.aabb1.Area() + _.aabb2.Area());
					.ThenBy(_ => Math.Max(_.aabb1.Area(), _.aabb2.Area()));
					// .ThenBy(_ => Math.Abs(_.aabb1.Area() - _.aabb2.Area()))

				var index = list.Any()? list.First().index: 0;
				// Console.WriteLine("3 nodes : " + index.ToString());
				switch (index)
				{
					case 1: Swap(nR, n1); break;
					case 2: Swap(nR, n2); break;
					case 3: Swap(nL, n3); break;
					case 4: Swap(nL, n4); break;
					default: break;
				}
			}

			// 4-nodes
			{
				Node nL = node.left;
				Node nR = node.right;
				Node n1 = nL.left;
				Node n2 = nL.right;
				Node n3 = nR.left;
				Node n4 = nR.right;

				var list = new[]
				{
					new {index = 0, n1 = n1, n2 = n2, n3 = n3, n4 = n4},
					new {index = 1, n1 = n1, n2 = n3, n3 = n2, n4 = n4},
					new {index = 2, n1 = n1, n2 = n4, n3 = n3, n4 = n2},
				}.ToList()
					.Where(_ => (_.n1 != null && _.n2 != null && _.n3 != null && _.n4 != null))
					.Where(_ => b(_.n1, _.n2) && b(_.n3, _.n4) && b(u(_.n1, _.n2), u(_.n3, _.n4)))
					.Select(_ => new {aabb1 = (_.n1.aabb | _.n2.aabb), aabb2 = (_.n3.aabb | _.n4.aabb), index = _.index})
					.OrderBy(_ => (_.aabb1 & _.aabb2).Area())
					.ThenBy(_ => _.aabb1.Area() + _.aabb2.Area());
					// .ThenBy(_ => Math.Max(_.aabb1.Area(), _.aabb2.Area()))
					// .ThenBy(_ => Math.Abs(_.aabb1.Area() - _.aabb2.Area()))
				
				var index = list.Any()? list.First().index: 0;
				// Console.WriteLine("4 nodes : " + index.ToString());
				switch (index)
				{
					case 1: Swap(n2, n3); break;
					case 2: Swap(n2, n4); break;
					default: break;
				}
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

		void _Draw(Graphics g, Node node, int depth, int scale)
		{
			if (node == null) 
			{
				return;
			}

			AABB aabb = node.aabb;
			bool isLeaf = node.isLeaf();
			int height = node.height;

			int x = Convert.ToInt32(scale * aabb.lb[0]);
			int y = Convert.ToInt32(scale * aabb.lb[1]);
			int w = Convert.ToInt32(scale * (aabb.ub[0] - aabb.lb[0]));
			int h = Convert.ToInt32(scale * (aabb.ub[1] - aabb.lb[1]));

			Color color = isLeaf? Color.Blue: Color.Black;
			int penWidth = 1; //isLeaf? 1: Math.Max(1, 5 - depth);
			Pen pen = new Pen(color, penWidth);
			int rgb = Math.Min(255, 128 + depth * 20);
			SolidBrush brush = new SolidBrush(Color.FromArgb(rgb, rgb, rgb));

			if (isLeaf)
			{
				// Console.WriteLine(String.Format("{0} {1} {2} {3}", x, y, w, h));
				g.DrawEllipse(pen, x, y, w, h);
			} else {
				// Console.WriteLine("isLeaf false");
				//g.FillRectangle(brush, x, y, w, h);
				g.DrawRectangle(pen, x, y, w, h);
			}

			_Draw(g, node.right, depth + 1, scale);
			_Draw(g, node.left, depth + 1, scale);
		}

		public void Draw(String file, double baseLength, int scale)
		{
			Bitmap bitmap = new Bitmap(
				Convert.ToInt32(baseLength*scale), 
				Convert.ToInt32(baseLength*scale), 
				System.Drawing.Imaging.PixelFormat.Format32bppArgb);
			
			Graphics g = Graphics.FromImage(bitmap);
			g.Clear(Color.White);

			_Draw(g, root, 0, scale);
			bitmap.Save(file, System.Drawing.Imaging.ImageFormat.Png);
		}

	}

}
