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
			
			double maxDisp = 0.1; // box gap
			double radius = 0.5; // radius of circles
			double baseLength = 1024; // max box
			List<bool> periodicity = new List<bool>() {false, false};
			List<double> boxSize = new List<double>() {baseLength, baseLength};

			BVH.Tree tree = new BVH.Tree(2, maxDisp, periodicity, boxSize, n);

			for (int i = 0; i < n; i++)
			{
				List<double> position = _.List<double>(2);
				position[0] = boxSize[0] * (random.NextDouble() * 0.8 + 0.1);
				position[1] = boxSize[1] * (random.NextDouble() * 0.8 + 0.1);

				// Insert particle into tree.
				// Console.Write("Insert\n");
				tree.insertParticle(i, position, radius * (1 + 5 * random.NextDouble()));
				//Thread.Sleep(50);
			}
			// Console.Write("Tree generated!\n");

			tree.Draw(@"img/test.png", baseLength, 1);
			//Console.WriteLine(tree);
		}

	}

}
