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
			double gap = 0.1; // box gap
			double radius = 5; // radius of circles
			double size = 1024; // size of world

			BVH.Tree tree = new BVH.Tree(dim, gap);

			for (int i = 0; i < n; i++)
			{
				List<double> center = _.List<double>(2);
				center[0] = size * (random.NextDouble() * 0.8 + 0.1);
				center[1] = size * (random.NextDouble() * 0.8 + 0.1);

				// Console.Write("Insert\n");
				tree.Insert(new Node(new AABB(2, gap, center, radius), null));
				//Thread.Sleep(50);
			}
			// Console.Write("Tree generated!\n");

			tree.Draw(@"img/test.png", size, 1);
			//Console.WriteLine(tree);
		}

	}

}
