using System; 
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Threading;

namespace aabb
{
	/*! \brief A node of the AABB tree.

		Each node of the tree contains an AABB object which corresponds to a
		particle, or a group of particles, in the simulation box. The AABB
		objects of individual particles are "fattened" before they are stored
		to avoid having to continually update and rebalance the tree when
		displacements are small.

		Nodes are aware of their position within in the tree. The isLeaf member
		function allows the tree to query whether the node is a leaf, i.e. to
		determine whether it holds a single particle.
	 */
	public class Node
	{
		/// Constructor.
		public Node()
		{
		}

		/// The fattened axis-aligned bounding box.
		public AABB aabb = new AABB();

		/// Index of the parent node.
		public int parent;

		/// Index of the next node.
		public int next;

		/// Index of the left-hand child.
		public int left;

		/// Index of the right-hand child.
		public int right;

		/// Height of the node. This is 0 for a leaf and -1 for a free node.
		public int height;

		/// The index of the particle that the node contains (leaf nodes only).
		public int particle;

		//! Test whether the node is a leaf.
		/*! \return
				Whether the node is a leaf node.
		 */
		public bool isLeaf()
		{
			return (left == -1);
		}
	}

}
