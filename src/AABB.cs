using System; 
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Threading;

namespace aabb
{
	/*! \brief The axis-aligned bounding box object.

	    Axis-aligned bounding boxes (AABBs) store information for the minimum
	    orthorhombic bounding-box for an object. Support is provided for
	    dimensions >= 2. (In 2D the bounding box is either a rectangle,
	    in 3D it is a rectangular prism.)

	    Class member functions provide functionality for merging AABB objects
	    and testing overlap with other AABBs.
	 */
	public class AABB
	{
		/// Lower bound of AABB in each dimension.
		public List<double> lowerBound = _.List<double>();

		/// Upper bound of AABB in each dimension.
		public List<double> upperBound = _.List<double>();

		/// The position of the AABB centre.
		public List<double> centre = _.List<double>();

		/// The AABB's surface area.
		public double surfaceArea;

		public static AABB zero = new AABB(0);

		/// Constructor.
		public AABB()
		{
		}

		//! Constructor.
		/*! \param dimension
		        The dimensionality of the system.
		 */
		public AABB(int dimension)
		{
			//Debug.Assert(dimension >= 2);
			lowerBound.Resize(dimension);
			upperBound.Resize(dimension);
		}

		//! Constructor.
		/*! \param lowerBound_
		        The lower bound in each dimension.

		    \param upperBound_
		        The upper bound in each dimension.
		 */
		public AABB(List<double> lowerBound_, List<double> upperBound_)
		{
			this.lowerBound = _.List<double>(lowerBound_);
			this.upperBound = _.List<double>(upperBound_);
			// Validate the dimensionality of the bounds vectors.
			if (lowerBound.Count != upperBound.Count)
			{
				throw new System.ArgumentException("[ERROR]: Dimensionality mismatch!");
			}

			// Validate that the upper bounds exceed the lower bounds.
			for (int i = 0;i < lowerBound.Count;i++)
			{
				// Validate the bound.
				if (lowerBound[i] > upperBound[i])
				{
					throw new System.ArgumentException("[ERROR]: AABB lower bound is greater than the upper bound!");
				}
			}

			surfaceArea = computeSurfaceArea();
			centre = _.List<double>(computeCentre());
		}

		public void fatten(double k)
		{
			for (int i = 0;i < lowerBound.Count;i++)
			{
				// double offset = k * (upperBound[i] - lowerBound[i]);
				double offset = 2;
				lowerBound[i] -= offset;
				upperBound[i] += offset;
			}
		}

		/// Compute the surface area of the box.
		public double computeSurfaceArea()
		{
			// Sum of "area" of all the sides.
			double sum = 0;

			// General formula for one side: hold one dimension constant
			// and multiply by all the other ones.
			for (int d1 = 0; d1 < lowerBound.Count; d1++)
			{
				// "Area" of current side.
				double product = 1;

				for (int d2 = 0; d2 < lowerBound.Count; d2++)
				{
					if (d1 == d2)
					{
						continue;
					}

					double dx = upperBound[d2] - lowerBound[d2];
					product *= dx;
				}

				// Update the sum.
				sum += product;
			}

			return 2.0 * sum;
		}

		public double Area()
		{
			return computeSurfaceArea();
		}
		
		/// Get the surface area of the box.
		public double getSurfaceArea()
		{
			return surfaceArea;
		}

		public static AABB operator| (AABB aabb1, AABB aabb2)
		{
			return Union(aabb1, aabb2);
		}

		public static AABB operator& (AABB aabb1, AABB aabb2)
		{
			return Intersect(aabb1, aabb2);
		}

		//! Union two AABBs into this one.
		/*! \param aabb1
		        A reference to the first AABB.

		    \param aabb2
		        A reference to the second AABB.
		 */
		public static AABB Union(AABB aabb1, AABB aabb2)
		{
			if (aabb1.lowerBound.Count != aabb2.lowerBound.Count) return AABB.zero;
			if (aabb1.upperBound.Count != aabb2.upperBound.Count) return AABB.zero;

			int dimension = aabb1.lowerBound.Count;
			AABB aabb = new AABB(dimension);

			for (int i = 0;i < dimension;i++)
			{
				aabb.lowerBound[i] = Math.Min(aabb1.lowerBound[i], aabb2.lowerBound[i]);
				aabb.upperBound[i] = Math.Max(aabb1.upperBound[i], aabb2.upperBound[i]);
			}

			aabb.surfaceArea = aabb.computeSurfaceArea();
			aabb.centre = _.List<double>(aabb.computeCentre());

			return aabb;
		}

		public static AABB Intersect(AABB aabb1, AABB aabb2)
		{
			if (aabb1.lowerBound.Count != aabb2.lowerBound.Count) return AABB.zero;
			if (aabb1.upperBound.Count != aabb2.upperBound.Count) return AABB.zero;

			int dimension = aabb1.lowerBound.Count;
			AABB aabb = new AABB(dimension);

			for (int i = 0;i < dimension;i++)
			{
				aabb.lowerBound[i] = Math.Max(aabb1.lowerBound[i], aabb2.lowerBound[i]);
				aabb.upperBound[i] = Math.Min(aabb1.upperBound[i], aabb2.upperBound[i]);
				if (aabb.lowerBound[i] > aabb.upperBound[i])
					return AABB.zero;
			}

			aabb.surfaceArea = aabb.computeSurfaceArea();
			aabb.centre = _.List<double>(aabb.computeCentre());

			return aabb;
		}

		//! Test whether the AABB is contained within this one.
		/*! \param aabb
		        A reference to the AABB.

		    \return
		        Whether the AABB is fully contained.
		 */
		public bool contains(AABB aabb)
		{
			Debug.Assert(aabb.lowerBound.Count == lowerBound.Count);

			for (int i = 0;i < lowerBound.Count;i++)
			{
				if (aabb.lowerBound[i] < lowerBound[i])
				{
					return false;
				}
				if (aabb.upperBound[i] > upperBound[i])
				{
					return false;
				}
			}

			return true;
		}

		//! Test whether the AABB overlaps this one.
		/*! \param aabb
		        A reference to the AABB.

		    \param touchIsOverlap
		        Does touching constitute an overlap?

		    \return
		        Whether the AABB overlaps.
		 */
		public bool overlaps(AABB aabb, bool touchIsOverlap)
		{
			Debug.Assert(aabb.lowerBound.Count == lowerBound.Count);

			bool rv = true;

			if (touchIsOverlap)
			{
				for (int i = 0; i < lowerBound.Count; ++i)
				{
					if (aabb.upperBound[i] < lowerBound[i] || aabb.lowerBound[i] > upperBound[i])
					{
						rv = false;
						break;
					}
				}
			}
			else
			{
				for (int i = 0; i < lowerBound.Count; ++i)
				{
					if (aabb.upperBound[i] <= lowerBound[i] || aabb.lowerBound[i] >= upperBound[i])
					{
						rv = false;
						break;
					}
				}
			}

			return rv;
		}

		//! Compute the centre of the AABB.
		/*! \returns
		        The position vector of the AABB centre.
		 */
		public List<double> computeCentre()
		{
			List<double> position = _.List<double>(lowerBound.Count);

			for (int i = 0;i < position.Count;i++)
			{
				position[i] = 0.5 * (lowerBound[i] + upperBound[i]);
			}

			return _.List<double>(position);
		}

		//! Set the dimensionality of the AABB.
		/*! \param dimension
		        The dimensionality of the system.
		 */
		public void setDimension(int dimension)
		{
			Debug.Assert(dimension >= 2);

			lowerBound.Resize(dimension);
			upperBound.Resize(dimension);
		}

		public override string ToString()
		{
			if (lowerBound.Count != 2)
			{
				return "AABB";
			}
			return ("AABB " + 
				lowerBound[0].ToString("0.00").PadLeft(8) + " " + 
				lowerBound[1].ToString("0.00").PadLeft(8) + " | " + 
				upperBound[0].ToString("0.00").PadLeft(8) + " " + 
				upperBound[1].ToString("0.00").PadLeft(8));
		}
	}

}
