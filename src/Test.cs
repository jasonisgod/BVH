using System; 
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Threading;

namespace BVH
{
	public class Test
	{
		public static void Main(string[] args)
		{
			int n = args.Length > 0? Convert.ToInt32(args[0]): 100; // number of circles
			int seed = args.Length > 1? Convert.ToInt32(args[1]): System.DateTime.Now.GetHashCode(); // random seed

			Random random = new Random(seed);
			Console.WriteLine("Seed: " + seed.ToString());
			
			int dim = 2; // dimension
			double gap = 2; // margin gap
			double radius = 3; // radius of circles
			double size = 1024; // size of world

			BVH.Tree tree = new BVH.Tree(dim, gap);
			var nodes = new List<Node>();

			for (int i = 0; i < n; i++)
			{
				List<double> center = _.List<double>(2);
				center[0] = size * (random.NextDouble() * 0.8 + 0.1);
				center[1] = size * (random.NextDouble() * 0.8 + 0.1);
				Node node = new Node(new AABB(2, gap, center, radius), null);
				nodes.Add(node);
				tree.Insert(node);
				//Thread.Sleep(50);
			}

			DrawTree(tree, @"img/test.png", size, 1);

			while (true)
			{
				DrawTree(tree, @"img/test.png", size, 1);
				var line = Console.ReadLine().Split(new char[] {' '}, StringSplitOptions.RemoveEmptyEntries);
				if (line.Length == 2)
				{
					if (line[0][0] == 'i')
					{
						tree.Insert(nodes[Convert.ToInt32(line[1])]);
					}
					if (line[0][0] == 'r')
					{
						tree.Remove(nodes[Convert.ToInt32(line[1])]);
					}
				}
			}

		}
		
		public static void DrawTree(Tree tree, String file, double size, int scale)
		{
			int width = Convert.ToInt32(size * scale);
			int height = Convert.ToInt32(size * scale);
			Bitmap bitmap = new Bitmap(width, height); //, System.Drawing.Imaging.PixelFormat.Format32bppArgb
			Graphics g = Graphics.FromImage(bitmap);
			
			g.Clear(Color.White);
			_Draw(tree.root, 0);
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
	}
}
