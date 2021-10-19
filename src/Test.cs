using System; 
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Threading;

namespace aabb
{
    public class Test
    {
        public static void Main(string[] args)
        {
            /*****************************************************************/
            /*      Set parameters, initialise variables and objects.        */
            /*****************************************************************/

            int nSweeps = 100000; // The number of Monte Carlo sweeps.
            int sampleInterval = 100; // The number of sweeps per sample.
            int nSmall = Convert.ToInt32(args[0]); // The number of small particles.
            int nLarge = 0; // The number of large particles.
            double diameterSmall = 5; // The diameter of the small particles.
            double diameterLarge = 50; // The diameter of the large particles.
            //double density = 0.1; // The system density.
            double maxDisp = 0.1; // Maximum trial displacement (in units of diameter).

            // Total particles.
            int nParticles = nSmall + nLarge;

            // Number of samples.
            int nSamples = nSweeps / sampleInterval;

            // Particle radii.
            double radiusSmall = 0.5 * diameterSmall;
            double radiusLarge = 0.5 * diameterLarge;

            // Output formatting flag.
            int format = Convert.ToInt32(Math.Floor(Math.Log10(nSamples)));

            // Set the periodicity of the simulation box.
            List<bool> periodicity = new List<bool>() {false, false};

            // Work out base length of simulation box.
            double baseLength = 1024; //Math.Pow((Math.PI * (nSmall * diameterSmall + nLarge * diameterLarge)) / (4.0 * density), 1.0 / 2.0);
            List<double> boxSize = new List<double>() {baseLength, baseLength};

            // Initialise the random number generator.
            Random random = new Random();

            // Initialise the AABB trees.
            aabb.Tree treeSmall = new aabb.Tree(2, maxDisp, periodicity, boxSize, nSmall);

            // Initialise particle position vectors.
            List<List<double>> positionsSmall = _.List<List<double>>(nSmall, _.List<double>(boxSize.Count));

            /*****************************************************************/
            /*             Generate the initial AABB trees.                  */
            /*****************************************************************/

			Console.Write("\nInserting small particles into AABB tree ...\n");
			for (int i = 0;i < nSmall;i++)
            {
                // Initialise the particle position vector.
                List<double> position = _.List<double>(2);
				position[0] = boxSize[0] * (random.NextDouble() * 0.8 + 0.1);
				position[1] = boxSize[1] * (random.NextDouble() * 0.8 + 0.1);

                // Insert particle into tree.
                treeSmall.insertParticle(i, position, radiusSmall * (1 + 5 * random.NextDouble()));
				//Thread.Sleep(50);

                // Store the position.
                positionsSmall[i] = _.List<double>(position);
            }
            Console.Write("Tree generated!\n");

			treeSmall.Draw(@"test.png", baseLength, 1);
			//Console.WriteLine(treeSmall);
        }

    }

}
