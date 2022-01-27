using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Threading;

namespace BVH
{
	public class AABB
	{
		// public const double OFFSET = 2;

		public int dim;
		public double gap;
		public List<double> lb;
		public List<double> ub;

		public AABB(int dim_, double gap_)
		{
			Reset(dim_, gap_);
			for (int i = 0; i < dim; i++)
			{
				lb[i] = 0;
				ub[i] = 0;
			}
			Fatten();
		}

		public AABB(int dim_, double gap_, List<double> lb_, List<double> ub_)
		{
			Reset(dim_, gap_);
			for (int i = 0; i < dim; i++)
			{
				lb[i] = lb_[i];
				ub[i] = ub_[i];
			}
			Fatten();
		}

		public AABB(int dim_, double gap_, List<double> center, double radius)
		{
			Reset(dim_, gap_);
			for (int i = 0; i < dim; i++)
			{
				lb[i] = center[i] - radius;
				ub[i] = center[i] + radius;
			}
			Fatten();
		}

		private void Reset(int dim_, double gap_)
		{
			dim = dim_;
			gap = gap_;
			lb = _.List<double>(dim);
			ub = _.List<double>(dim);
		}

		private void Fatten()
		{
			for (int i = 0; i < dim; i++)
			{
				// double offset = gap * (ub[i] - lb[i]);
				lb[i] -= gap;
				ub[i] += gap;
			}
		}

		public double Area()
		{
			double sum = 0;

			for (int d1 = 0; d1 < dim; d1++)
			{
				double product = 1;
				for (int d2 = 0; d2 < dim; d2++)
				{
					if (d1 == d2)
					{
						continue;
					}
					product *= ub[d2] - lb[d2];
				}
				sum += product;
			}

			return 2.0 * sum;
		}

		public static AABB operator| (AABB x, AABB y)
		{
			AABB aabb = new AABB(x.dim, x.gap);

			for (int i = 0; i < x.dim; i++)
			{
				aabb.lb[i] = Math.Min(x.lb[i], y.lb[i]);
				aabb.ub[i] = Math.Max(x.ub[i], y.ub[i]);
			}

			aabb.Fatten();
			return aabb;
		}

		public static AABB operator& (AABB x, AABB y)
		{
			AABB aabb = new AABB(x.dim, x.gap);

			for (int i = 0; i < x.dim; i++)
			{
				aabb.lb[i] = Math.Max(x.lb[i], y.lb[i]);
				aabb.ub[i] = Math.Min(x.ub[i], y.ub[i]);

				if (aabb.lb[i] > aabb.ub[i])
				{
					return new AABB(x.dim, x.gap);
				}
			}
			
			aabb.Fatten();
			return aabb;
		}

		public bool Contains(AABB aabb)
		{
			for (int i = 0; i < dim; i++)
			{
				if (aabb.lb[i] < lb[i])
				{
					return false;
				}
				if (aabb.ub[i] > ub[i])
				{
					return false;
				}
			}
			return true;
		}

		public bool Overlaps(AABB aabb, bool touch)
		{
			bool ov = true;

			for (int i = 0; i < dim && ov; ++i)
			{
				if (aabb.ub[i] < lb[i] || aabb.lb[i] > ub[i])
				{
					ov = false;
				}
				if (touch && (aabb.ub[i] == lb[i] || aabb.lb[i] == ub[i]))
				{
					ov = false;
				}
			}
			
			return ov;
		}

		public List<double> Center()
		{
			List<double> center = _.List<double>(dim);

			for (int i = 0; i < dim; i++)
			{
				center[i] = 0.5 * (lb[i] + ub[i]);
			}

			return center;
		}
		
		public override string ToString()
		{
			if (dim != 2)
			{
				return "AABB";
			}
			return ("AABB " + 
				lb[0].ToString("0.00").PadLeft(8) + " " + 
				lb[1].ToString("0.00").PadLeft(8) + " | " + 
				ub[0].ToString("0.00").PadLeft(8) + " " + 
				ub[1].ToString("0.00").PadLeft(8));
		}
	}
}
