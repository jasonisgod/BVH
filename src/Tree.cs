using System; 
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Threading;

namespace aabb
{
    public class Tree
	{
		/// The index of the root node.
		public int root;

		/// The dynamic tree.
		public List<Node> nodes = _.List<Node>();

		/// The current number of nodes in the tree.
		private int nodeCount;

		/// The current node capacity.
		private int nodeCapacity;

		/// The position of node at the top of the free list.
		private int freeList;

		/// The dimensionality of the system.
		private int dimension;

		/// Whether the system is periodic along at least one axis.
		private bool isPeriodic;

		/// The skin thickness of the fattened AABBs, as a fraction of the AABB base length.
		private double skinThickness;

		/// Whether the system is periodic along each axis.
		private List<bool> periodicity = _.List<bool>();

		/// The size of the system in each dimension.
		private List<double> boxSize = _.List<double>();

		/// The position of the negative minimum image.
		private List<double> negMinImage = _.List<double>();

		/// The position of the positive minimum image.
		private List<double> posMinImage = _.List<double>();

		/// A map between particle and node indices.
		private Dictionary<int, int> particleMap = new Dictionary<int, int>();

		/// Does touching count as overlapping in tree queries?
		private bool touchIsOverlap;

		//! Constructor (non-periodic).
		/*! \param dimension_
		        The dimensionality of the system.

		    \param skinThickness_
		        The skin thickness for fattened AABBs, as a fraction
		        of the AABB base length.

		    \param nParticles
		        The number of particles (for fixed particle number systems).

		    \param touchIsOverlap
		        Does touching count as overlapping in query operations?
		 */
		public Tree(int dimension_ = 3, double skinThickness_ = 0.05, int nParticles = 16, bool touchIsOverlap_ = true)
		{
			this.dimension = dimension_;
			this.isPeriodic = false;
			this.skinThickness = skinThickness_;
			this.touchIsOverlap = touchIsOverlap_;
			// Validate the dimensionality.
			if ((dimension < 2))
			{
				throw new System.ArgumentException("[ERROR]: Invalid dimensionality!");
			}

			// Initialise the periodicity vector.
			periodicity.Resize(dimension, false);

			// Initialise the tree.
			root = -1;
			nodeCount = 0;
			nodeCapacity = nParticles * 2;
			nodes.Resize(nodeCapacity);

			// Build a linked list for the list of free nodes.
			for (int i = 0;i < nodeCapacity - 1;i++)
			{
				nodes[i].next = i + 1;
				nodes[i].height = -1;
			}
			nodes[nodeCapacity - 1].next = -1;
			nodes[nodeCapacity - 1].height = -1;

			// Assign the index of the first free node.
			freeList = 0;
		}

		//! Constructor (custom periodicity).
		/*! \param dimension_
		        The dimensionality of the system.

		    \param skinThickness_
		        The skin thickness for fattened AABBs, as a fraction
		        of the AABB base length.

		    \param periodicity_
		        Whether the system is periodic in each dimension.

		    \param boxSize_
		        The size of the simulation box in each dimension.

		    \param nParticles
		        The number of particles (for fixed particle number systems).

		    \param touchIsOverlap
		        Does touching count as overlapping in query operations?
		 */
		public Tree(int dimension_, double skinThickness_, List<bool> periodicity_, List<double> boxSize_, int nParticles = 16, bool touchIsOverlap_ = true)
		{
			this.dimension = dimension_;
			this.skinThickness = skinThickness_;
			this.periodicity = _.List<bool>(periodicity_);
			this.boxSize = _.List<double>(boxSize_);
			this.touchIsOverlap = touchIsOverlap_;
			// Validate the dimensionality.
			if (dimension < 2)
			{
				throw new System.ArgumentException("[ERROR]: Invalid dimensionality!");
			}

			// Validate the dimensionality of the vectors.
			if ((periodicity.Count != dimension) || (boxSize.Count != dimension))
			{
				throw new System.ArgumentException("[ERROR]: Dimensionality mismatch!");
			}

			// Initialise the tree.
			root = -1;
			nodeCount = 0;
			nodeCapacity = nParticles * 2;
			nodes.Resize(nodeCapacity);

			// Build a linked list for the list of free nodes.
			for (int i = 0;i < nodeCapacity - 1;i++)
			{
				nodes[i].next = i + 1;
				nodes[i].height = -1;
			}
			nodes[nodeCapacity - 1].next = -1;
			nodes[nodeCapacity - 1].height = -1;

			// Assign the index of the first free node.
			freeList = 0;

			// Check periodicity.
			isPeriodic = false;
			posMinImage.Resize(dimension);
			negMinImage.Resize(dimension);
			for (int i = 0;i < dimension;i++)
			{
				posMinImage[i] = 0.5 * boxSize[i];
				negMinImage[i] = -0.5 * boxSize[i];

				if (periodicity[i])
				{
					isPeriodic = true;
				}
			}
		}

		//! Set the periodicity of the simulation box.
		/*! \param periodicity_
		        Whether the system is periodic in each dimension.
		 */
		public void setPeriodicity(List<bool> periodicity_)
		{
			periodicity = _.List<bool>(periodicity_);
		}

		//! Set the size of the simulation box.
		/*! \param boxSize_
		        The size of the simulation box in each dimension.
		 */
		public void setBoxSize(List<double> boxSize_)
		{
			boxSize = _.List<double>(boxSize_);
		}

		//! Allocate a new node.
		/*! \return
		        The index of the allocated node.
		 */
		private int allocateNode()
		{
			// Exand the node pool as needed.
			if (freeList == -1)
			{
				Debug.Assert(nodeCount == nodeCapacity);

				// The free list is empty. Rebuild a bigger pool.
				nodeCapacity *= 2;
				nodes.Resize(nodeCapacity);

				// Build a linked list for the list of free nodes.
				for (int i = nodeCount;i < nodeCapacity - 1;i++)
				{
					nodes[i].next = i + 1;
					nodes[i].height = -1;
				}
				nodes[nodeCapacity - 1].next = -1;
				nodes[nodeCapacity - 1].height = -1;

				// Assign the index of the first free node.
				freeList = nodeCount;
			}

			// Peel a node off the free list.
			int node = freeList;
			freeList = nodes[node].next;
			nodes[node].parent = -1;
			nodes[node].left = -1;
			nodes[node].right = -1;
			nodes[node].height = 0;
			nodes[node].aabb.setDimension(dimension);
			nodeCount++;

			return node;
		}

		//! Free an existing node.
		/*! \param node
		        The index of the node to be freed.
		 */
		private void freeNode(int node)
		{
			Debug.Assert(node < nodeCapacity);
			Debug.Assert(0 < nodeCount);

			nodes[node].next = freeList;
			nodes[node].height = -1;
			freeList = node;
			nodeCount--;
		}

		//! Insert a particle into the tree (point particle).
		/*! \param index
		        The index of the particle.

		    \param position
		        The position vector of the particle.

		    \param radius
		        The radius of the particle.
		 */
		public void insertParticle(int particle, List<double> position, double radius)
		{
			// Make sure the particle doesn't already exist.
			if (particleMap.Count(particle) != 0)
			{
				throw new System.ArgumentException("[ERROR]: Particle already exists in tree!");
			}

			// Validate the dimensionality of the position vector.
			if (position.Count != dimension)
			{
				throw new System.ArgumentException("[ERROR]: Dimensionality mismatch!");
			}

			// Allocate a new node for the particle.
			int node = allocateNode();

			// AABB size in each dimension.
			List<double> size = _.List<double>(dimension);

			// Compute the AABB limits.
			for (int i = 0;i < dimension;i++)
			{
				nodes[node].aabb.lowerBound[i] = position[i] - radius;
				nodes[node].aabb.upperBound[i] = position[i] + radius;
				size[i] = nodes[node].aabb.upperBound[i] - nodes[node].aabb.lowerBound[i];
			}

			// Fatten the AABB.
			nodes[node].aabb.fatten(skinThickness);

			nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
			nodes[node].aabb.centre = nodes[node].aabb.computeCentre();

			// Zero the height.
			nodes[node].height = 0;

			// Insert a new leaf into the tree.
			insertLeaf(node);

			// Add the new particle to the map.
			particleMap.Add(particle, node);

			// Store the particle index.
			nodes[node].particle = particle;
		}

		//! Insert a particle into the tree (arbitrary shape with bounding box).
		/*! \param index
		        The index of the particle.

		    \param lowerBound
		        The lower bound in each dimension.

		    \param upperBound
		        The upper bound in each dimension.
		 */
		public void insertParticle(int particle, List<double> lowerBound, List<double> upperBound)
		{
			// Make sure the particle doesn't already exist.
			if (particleMap.Count(particle) != 0)
			{
				throw new System.ArgumentException("[ERROR]: Particle already exists in tree!");
			}

			// Validate the dimensionality of the bounds vectors.
			if ((lowerBound.Count != dimension) || (upperBound.Count != dimension))
			{
				throw new System.ArgumentException("[ERROR]: Dimensionality mismatch!");
			}

			// Allocate a new node for the particle.
			int node = allocateNode();

			// AABB size in each dimension.
			List<double> size = _.List<double>(dimension);

			// Compute the AABB limits.
			for (int i = 0;i < dimension;i++)
			{
				// Validate the bound.
				if (lowerBound[i] > upperBound[i])
				{
					throw new System.ArgumentException("[ERROR]: AABB lower bound is greater than the upper bound!");
				}

				nodes[node].aabb.lowerBound[i] = lowerBound[i];
				nodes[node].aabb.upperBound[i] = upperBound[i];
				size[i] = upperBound[i] - lowerBound[i];
			}

			// Fatten the AABB.
			nodes[node].aabb.fatten(skinThickness);

			nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
			nodes[node].aabb.centre = nodes[node].aabb.computeCentre();

			// Zero the height.
			nodes[node].height = 0;

			// Insert a new leaf into the tree.
			insertLeaf(node);

			// Add the new particle to the map.
			particleMap.Add(particle, node);

			// Store the particle index.
			nodes[node].particle = particle;
		}

        /// Return the number of particles in the tree.
		public int nParticles()
		{
			return (int)particleMap.Count;
		}

		//! Remove a particle from the tree.
		/*! \param particle
		        The particle index (particleMap will be used to map the node).
		 */
		public void removeParticle(int particle)
		{
			// The particle doesn't exist.
			if (!particleMap.ContainsKey(particle))
			{
				throw new System.ArgumentException("[ERROR]: Invalid particle index!");
			}

            // Extract the node index.
            int node = particleMap[particle];

			// Erase the particle from the map.
			particleMap.Remove(particle);

			Debug.Assert(node < nodeCapacity);
			Debug.Assert(nodes[node].isLeaf());

			removeLeaf(node);
			freeNode(node);
		}

		/// Remove all particles from the tree.
		public void removeAll()
		{
			// Iterator pointing to the start of the particle map.
			Dictionary<int, int>.Enumerator it = particleMap.GetEnumerator();

			// Iterate over the map.
			while (it.MoveNext())
			{
				// Extract the node index.
				int node = it.Current.Value;

				Debug.Assert(node < nodeCapacity);
				Debug.Assert(nodes[node].isLeaf());

				removeLeaf(node);
				freeNode(node);

			}

			// Clear the particle map.
			particleMap.Clear();
		}

		//! Update the tree if a particle moves outside its fattened AABB.
		/*! \param particle
		        The particle index (particleMap will be used to map the node).

		    \param position
		        The position vector of the particle.

		    \param radius
		        The radius of the particle.

		    \param alwaysReinsert
		        Always reinsert the particle, even if it's within its old AABB (default:false)

		    \return
		        Whether the particle was reinserted.
		 */
		public bool updateParticle(int particle, List<double> position, double radius, bool alwaysReinsert = false)
		{
			// Validate the dimensionality of the position vector.
			if (position.Count != dimension)
			{
				throw new System.ArgumentException("[ERROR]: Dimensionality mismatch!");
			}

			// AABB bounds vectors.
			List<double> lowerBound = _.List<double>(dimension);
			List<double> upperBound = _.List<double>(dimension);

			// Compute the AABB limits.
			for (int i = 0;i < dimension;i++)
			{
				lowerBound[i] = position[i] - radius;
				upperBound[i] = position[i] + radius;
			}

			// Update the particle.
			return updateParticle(particle, lowerBound, upperBound, alwaysReinsert);
		}

		//! Update the tree if a particle moves outside its fattened AABB.
		/*! \param particle
		        The particle index (particleMap will be used to map the node).

		    \param lowerBound
		        The lower bound in each dimension.

		    \param upperBound
		        The upper bound in each dimension.

		    \param alwaysReinsert
		        Always reinsert the particle, even if it's within its old AABB (default: false)
		 */
		public bool updateParticle(int particle, List<double> lowerBound, List<double> upperBound, bool alwaysReinsert = false)
		{
			// Validate the dimensionality of the bounds vectors.
			if ((lowerBound.Count != dimension) && (upperBound.Count != dimension))
			{
				throw new System.ArgumentException("[ERROR]: Dimensionality mismatch!");
			}

			// The particle doesn't exist.
			if (!particleMap.ContainsKey(particle))
			{
				throw new System.ArgumentException("[ERROR]: Invalid particle index!");
			}

            // Extract the node index.
            int node = particleMap[particle];

			// Erase the particle from the map.
			particleMap.Remove(particle);

			Debug.Assert(node < nodeCapacity);
			Debug.Assert(nodes[node].isLeaf());

			// AABB size in each dimension.
			List<double> size = _.List<double>(dimension);

			// Compute the AABB limits.
			for (int i = 0;i < dimension;i++)
			{
				// Validate the bound.
				if (lowerBound[i] > upperBound[i])
				{
					throw new System.ArgumentException("[ERROR]: AABB lower bound is greater than the upper bound!");
				}

				size[i] = upperBound[i] - lowerBound[i];
			}

			// Create the new AABB.
			AABB aabb = new AABB(lowerBound, upperBound);

			// No need to update if the particle is still within its fattened AABB.
			if (!alwaysReinsert && nodes[node].aabb.contains(aabb))
			{
				return false;
			}

			// Remove the current leaf.
			removeLeaf(node);

			// Fatten the new AABB.
			aabb.fatten(skinThickness);

			// Assign the new AABB.
			nodes[node].aabb = aabb;

			// Update the surface area and centroid.
			nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
			nodes[node].aabb.centre = nodes[node].aabb.computeCentre();

			// Insert a new leaf node.
			insertLeaf(node);

			return true;
		}

		//! Query the tree to find candidate interactions for a particle.
		/*! \param particle
		        The particle index.

		    \return particles
		        A vector of particle indices.
		 */
		public List<int> query(int particle)
		{
			// Make sure that this is a valid particle.
			if (particleMap.Count(particle) == 0)
			{
				throw new System.ArgumentException("[ERROR]: Invalid particle index!");
			}

			// Test overlap of particle AABB against all other particles.
			return _.List<int>(query(particle, nodes[particleMap[particle]].aabb));
		}

		//! Query the tree to find candidate interactions for an AABB.
		/*! \param particle
		        The particle index.

		    \param aabb
		        The AABB.

		    \return particles
		        A vector of particle indices.
		 */
		public List<int> query(int particle, AABB aabb)
		{
			List<int> stack = _.List<int>();
			stack.Capacity = 256;
			stack.Add(root);

			List<int> particles = _.List<int>();

			while (stack.Count > 0)
			{
				int node = stack[stack.Count - 1];
				stack.RemoveAt(stack.Count - 1);

				// Copy the AABB.
				AABB nodeAABB = nodes[node].aabb;

				if (node == -1)
				{
					continue;
				}

				if (isPeriodic)
				{
					List<double> separation = _.List<double>(dimension);
					List<double> shift = _.List<double>(dimension);
					for (int i = 0;i < dimension;i++)
					{
						separation[i] = nodeAABB.centre[i] - aabb.centre[i];
					}

					bool isShifted = minimumImage(separation, shift);

					// Shift the AABB.
					if (isShifted)
					{
						for (int i = 0;i < dimension;i++)
						{
							nodeAABB.lowerBound[i] += shift[i];
							nodeAABB.upperBound[i] += shift[i];
						}
					}
				}

				// Test for overlap between the AABBs.
				if (aabb.overlaps(nodeAABB, touchIsOverlap))
				{
					// Check that we're at a leaf node.
					if (nodes[node].isLeaf())
					{
						// Can't interact with itself.
						if (nodes[node].particle != particle)
						{
							particles.Add(nodes[node].particle);
						}
					}
					else
					{
						stack.Add(nodes[node].left);
						stack.Add(nodes[node].right);
					}
				}
			}

			return _.List<int>(particles);
		}

		//! Query the tree to find candidate interactions for an AABB.
		/*! \param aabb
		        The AABB.

		    \return particles
		        A vector of particle indices.
		 */
		public List<int> query(AABB aabb)
		{
			// Make sure the tree isn't empty.
			if (particleMap.Count == 0)
			{
				return _.List<int>();
			}

			// Test overlap of AABB against all particles.
			return _.List<int>(query(int.MaxValue, aabb));
		}

		//! Get a particle AABB.
		/*! \param particle
		        The particle index.
		 */
		public AABB getAABB(int particle)
		{
			return nodes[particleMap[particle]].aabb;
		}
        
        //! Insert a leaf into the tree.
        /*! \param leaf
                The index of the leaf node.
         */
		private void insertLeaf(int leaf)
		{
			if (root == -1)
			{
				root = leaf;
				nodes[root].parent = -1;
				return;
			}

			// Find the best sibling for the node.

			AABB leafAABB = nodes[leaf].aabb;
			int index = root;

			while (!nodes[index].isLeaf())
			{
				// Extract the children of the node.
				int left = nodes[index].left;
				int right = nodes[index].right;

				double surfaceArea = nodes[index].aabb.getSurfaceArea();

				AABB combinedAABB = AABB.Union(nodes[index].aabb, leafAABB);
				double combinedSurfaceArea = combinedAABB.getSurfaceArea();

				// Cost of creating a new parent for this node and the new leaf.
				double cost = 2.0 * combinedSurfaceArea;

				// Minimum cost of pushing the leaf further down the tree.
				double inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

				// Cost of descending to the left.
				double costLeft;
				if (nodes[left].isLeaf())
				{
					AABB aabb = AABB.Union(leafAABB, nodes[left].aabb);
					costLeft = aabb.getSurfaceArea() + inheritanceCost;
				}
				else
				{
					AABB aabb = AABB.Union(leafAABB, nodes[left].aabb);
					double oldArea = nodes[left].aabb.getSurfaceArea();
					double newArea = aabb.getSurfaceArea();
					costLeft = (newArea - oldArea) + inheritanceCost;
				}

				// Cost of descending to the right.
				double costRight;
				if (nodes[right].isLeaf())
				{
					AABB aabb = AABB.Union(leafAABB, nodes[right].aabb);
					costRight = aabb.getSurfaceArea() + inheritanceCost;
				}
				else
				{
					AABB aabb = AABB.Union(leafAABB, nodes[right].aabb);
					double oldArea = nodes[right].aabb.getSurfaceArea();
					double newArea = aabb.getSurfaceArea();
					costRight = (newArea - oldArea) + inheritanceCost;
				}

				// Descend according to the minimum cost.
				if ((cost < costLeft) && (cost < costRight))
				{
					break;
				}

				// Descend.
				if (costLeft < costRight)
				{
					index = left;
				}
				else
				{
					index = right;
				}
			}

			int sibling = index;

			// Create a new parent.
			int oldParent = nodes[sibling].parent;
			int newParent = allocateNode();
			nodes[newParent].parent = oldParent;
			nodes[newParent].aabb = AABB.Union(leafAABB, nodes[sibling].aabb);
			nodes[newParent].aabb.fatten(skinThickness);
			nodes[newParent].height = nodes[sibling].height + 1;

			// The sibling was not the root.
			if (oldParent != -1)
			{
				if (nodes[oldParent].left == sibling)
				{
					nodes[oldParent].left = newParent;
				}
				else
				{
					nodes[oldParent].right = newParent;
				}

				nodes[newParent].left = sibling;
				nodes[newParent].right = leaf;
				nodes[sibling].parent = newParent;
				nodes[leaf].parent = newParent;
			}
			// The sibling was the root.
			else
			{
				nodes[newParent].left = sibling;
				nodes[newParent].right = leaf;
				nodes[sibling].parent = newParent;
				nodes[leaf].parent = newParent;
				root = newParent;
			}
			
			// Walk back up the tree fixing heights and AABBs.
			index = nodes[leaf].parent;
			while (index != -1)
			{
				index = balance(index);
				update(index);
				optimize(index);
				// optimize2(index);
				index = nodes[index].parent;
			}
			
        }

		//! Remove a leaf from the tree.
		/*! \param leaf
		        The index of the leaf node.
		 */
		private void removeLeaf(int leaf)
		{
			if (leaf == root)
			{
				root = -1;
				return;
			}

			int parent = nodes[leaf].parent;
			int grandParent = nodes[parent].parent;
			int sibling;

			if (nodes[parent].left == leaf)
			{
				sibling = nodes[parent].right;
			}
			else
			{
				sibling = nodes[parent].left;
			}

			// Destroy the parent and connect the sibling to the grandparent.
			if (grandParent != -1)
			{
				if (nodes[grandParent].left == parent)
				{
					nodes[grandParent].left = sibling;
				}
				else
				{
					nodes[grandParent].right = sibling;
				}

				nodes[sibling].parent = grandParent;
				freeNode(parent);

				// Adjust ancestor bounds.
				int index = grandParent;
				while (index != -1)
				{
					index = balance(index);
					update(index);
					optimize(index);
					//optimize2(index);
					index = nodes[index].parent;
				}
			}
			else
			{
				root = sibling;
				nodes[sibling].parent = -1;
				freeNode(parent);
			}
		}

		private void update(int node)
		{
			if (node == -1) return;
			if (nodes[node].isLeaf()) return;

			int left = nodes[node].left;
			int right = nodes[node].right;

			nodes[node].height = 1 + Math.Max(nodes[left].height, nodes[right].height);
			nodes[node].aabb = AABB.Union(nodes[left].aabb, nodes[right].aabb);
			nodes[node].aabb.fatten(skinThickness);
		}

		private void swap(int node1, int node2)
		{
			if (node1 == -1 || node2 == -1) return;

			int par1 = nodes[node1].parent;
			int par2 = nodes[node2].parent;
			if (par1 == -1 || par2 == -1) return;

			nodes[node1].parent = par2;
			nodes[node2].parent = par1;
			if (nodes[par1].left == node1) nodes[par1].left = node2;
				else nodes[par1].right = node2;
			if (nodes[par2].left == node2) nodes[par2].left = node1;
				else nodes[par2].right = node1;
		}

		private int balance(int node)
		{
			//return node;
			Debug.Assert(node != -1);

			if (nodes[node].isLeaf() || (nodes[node].height < 2))
			{
				return node;
			}

			int left = nodes[node].left;
			int right = nodes[node].right;

			Debug.Assert(left < nodeCapacity);
			Debug.Assert(right < nodeCapacity);

			int currentBalance = nodes[right].height - nodes[left].height;

			// Rotate right branch up.
			if (currentBalance > 1)
			{
				int rightLeft = nodes[right].left;
				int rightRight = nodes[right].right;

				Debug.Assert(rightLeft < nodeCapacity);
				Debug.Assert(rightRight < nodeCapacity);

				// Swap node and its right-hand child.
				nodes[right].left = node;
				nodes[right].parent = nodes[node].parent;
				nodes[node].parent = right;

				// The node's old parent should now point to its right-hand child.
				if (nodes[right].parent != -1)
				{
					if (nodes[nodes[right].parent].left == node)
					{
						nodes[nodes[right].parent].left = right;
					}
					else
					{
						Debug.Assert(nodes[nodes[right].parent].right == node);
						nodes[nodes[right].parent].right = right;
					}
				}
				else
				{
					root = right;
				}

				// Rotate.
				if (nodes[rightLeft].height > nodes[rightRight].height)
				{
					nodes[right].right = rightLeft;
					nodes[node].right = rightRight;
					nodes[rightRight].parent = node;
					nodes[node].aabb = AABB.Union(nodes[left].aabb, nodes[rightRight].aabb);
					nodes[right].aabb = AABB.Union(nodes[node].aabb, nodes[rightLeft].aabb);

					nodes[node].height = 1 + Math.Max(nodes[left].height, nodes[rightRight].height);
					nodes[right].height = 1 + Math.Max(nodes[node].height, nodes[rightLeft].height);
				}
				else
				{
					nodes[right].right = rightRight;
					nodes[node].right = rightLeft;
					nodes[rightLeft].parent = node;
					nodes[node].aabb = AABB.Union(nodes[left].aabb, nodes[rightLeft].aabb);
					nodes[right].aabb = AABB.Union(nodes[node].aabb, nodes[rightRight].aabb);

					nodes[node].height = 1 + Math.Max(nodes[left].height, nodes[rightLeft].height);
					nodes[right].height = 1 + Math.Max(nodes[node].height, nodes[rightRight].height);
				}

				return right;
			}

			// Rotate left branch up.
			if (currentBalance < -1)
			{
				int leftLeft = nodes[left].left;
				int leftRight = nodes[left].right;

				Debug.Assert(leftLeft < nodeCapacity);
				Debug.Assert(leftRight < nodeCapacity);

				// Swap node and its left-hand child.
				nodes[left].left = node;
				nodes[left].parent = nodes[node].parent;
				nodes[node].parent = left;

				// The node's old parent should now point to its left-hand child.
				if (nodes[left].parent != -1)
				{
					if (nodes[nodes[left].parent].left == node)
					{
						nodes[nodes[left].parent].left = left;
					}
					else
					{
						Debug.Assert(nodes[nodes[left].parent].right == node);
						nodes[nodes[left].parent].right = left;
					}
				}
				else
				{
					root = left;
				}

				// Rotate.
				if (nodes[leftLeft].height > nodes[leftRight].height)
				{
					nodes[left].right = leftLeft;
					nodes[node].left = leftRight;
					nodes[leftRight].parent = node;
					nodes[node].aabb = AABB.Union(nodes[right].aabb, nodes[leftRight].aabb);
					nodes[left].aabb = AABB.Union(nodes[node].aabb, nodes[leftLeft].aabb);

					nodes[node].height = 1 + Math.Max(nodes[right].height, nodes[leftRight].height);
					nodes[left].height = 1 + Math.Max(nodes[node].height, nodes[leftLeft].height);
				}
				else
				{
					nodes[left].right = leftRight;
					nodes[node].left = leftLeft;
					nodes[leftLeft].parent = node;
					nodes[node].aabb = AABB.Union(nodes[right].aabb, nodes[leftLeft].aabb);
					nodes[left].aabb = AABB.Union(nodes[node].aabb, nodes[leftRight].aabb);

					nodes[node].height = 1 + Math.Max(nodes[right].height, nodes[leftLeft].height);
					nodes[left].height = 1 + Math.Max(nodes[node].height, nodes[leftRight].height);
				}

				return left;
			}
			
			return node;
		}
		
		private void optimize(int node)
		{
			if (node == -1) return;
			if (nodes[node].isLeaf()) return;

			int left = nodes[node].left;
			int right = nodes[node].right;

			if (nodes[left].isLeaf() || nodes[right].isLeaf()) return;
			
			int n1 = nodes[left].left;
			int n2 = nodes[left].right;
			int n3 = nodes[right].left;
			int n4 = nodes[right].right;

			AABB u12 = AABB.Union(nodes[n1].aabb, nodes[n2].aabb); double a12 = u12.Area();
			AABB u13 = AABB.Union(nodes[n1].aabb, nodes[n3].aabb); double a13 = u13.Area();
			AABB u14 = AABB.Union(nodes[n1].aabb, nodes[n4].aabb); double a14 = u14.Area();
			AABB u23 = AABB.Union(nodes[n2].aabb, nodes[n3].aabb); double a23 = u23.Area();
			AABB u24 = AABB.Union(nodes[n2].aabb, nodes[n4].aabb); double a24 = u24.Area();
			AABB u34 = AABB.Union(nodes[n3].aabb, nodes[n4].aabb); double a34 = u34.Area();

			double x1234 = AABB.Intersect(u12, u34).Area();
			double x1324 = AABB.Intersect(u13, u24).Area();
			double x1423 = AABB.Intersect(u14, u23).Area();
			double xMin = Math.Min(Math.Min(x1234, x1324), x1423);
			bool b1234 = (xMin == x1234);
			bool b1324 = (xMin == x1324);
			bool b1423 = (xMin == x1423);

			if (Convert.ToInt32(b1234) + Convert.ToInt32(b1324) + Convert.ToInt32(b1423) > 1)
			{
				double a1234 = b1234? u12.Area() + u34.Area(): double.MaxValue;
				double a1324 = b1324? u13.Area() + u24.Area(): double.MaxValue;
				double a1423 = b1423? u14.Area() + u23.Area(): double.MaxValue;
				double aMin = Math.Min(Math.Min(a1234, a1324), a1423);
				b1234 = (aMin == a1234);
				b1324 = (aMin == a1324);
				b1423 = (aMin == a1423);

				if (Convert.ToInt32(b1234) + Convert.ToInt32(b1324) + Convert.ToInt32(b1423) > 1)
				{
					double r1234 = b1234? Math.Max(a12, a34) / Math.Min(a12, a34): double.MaxValue;
					double r1324 = b1324? Math.Max(a13, a24) / Math.Min(a13, a24): double.MaxValue;
					double r1423 = b1423? Math.Max(a14, a23) / Math.Min(a14, a23): double.MaxValue;
					double rMin = Math.Min(Math.Min(r1234, r1324), r1423);
					b1234 = (rMin == r1234);
					b1324 = (rMin == r1324);
					b1423 = (rMin == r1423);
				} 
			}
			
			if (b1234) {
				// pass
			} else if (b1324) {
				swap(n2, n3);
			} else if (b1423) {
				swap(n2, n4);
			}
			
			update(left);
			update(right);
		}

		private void optimize2(int node)
		{
			if (node == -1) return;
			if (nodes[node].isLeaf()) return;

			for (int i = 0; i < 2; i++)
			{
				int left = nodes[node].left;
				int right = nodes[node].right;
				if (i == 1) 
				{
					int tmp = left; left = right; right = tmp;
				}

				if (nodes[left].isLeaf()) continue;
				if (nodes[left].height < nodes[right].height) continue;

				int n1 = nodes[left].left;
				int n2 = nodes[left].right;
				double aL = nodes[left].aabb.Area();
				double aR = nodes[right].aabb.Area();
				double a1 = nodes[n1].aabb.Area();
				double a2 = nodes[n2].aabb.Area();
				double aR1 = AABB.Union(nodes[right].aabb, nodes[n1].aabb).Area();
				double aR2 = AABB.Union(nodes[right].aabb, nodes[n2].aabb).Area();
				double aMin = Math.Min(Math.Min(a1 + aR2, a2 + aR1), aL + aR);
				if (aMin == a1 + aR2) swap(right, n1);
				if (aMin == a2 + aR1) swap(right, n2);
				update(left);
				update(right);
				// Console.Write(node);
				// Console.Write("\n");
				// Console.Write(aL + aR); Console.Write(" ");
				// Console.Write(a1 + aR2); Console.Write(" ");
				// Console.Write(a2 + aR1); Console.Write(" ");
				// Console.Write("\n");

			}
		}

		//! Compute the height of the tree.
		/*! \return
		        The height of the entire tree.
		 */
		private int computeHeight()
		{
			return computeHeight(root);
		}

		//! Compute the height of a sub-tree.
		/*! \param node
		        The index of the root node.

		    \return
		        The height of the sub-tree.
		 */
		private int computeHeight(int node)
		{
			Debug.Assert(node < nodeCapacity);

			if (nodes[node].isLeaf())
			{
				return 0;
			}

			int height1 = computeHeight(nodes[node].left);
			int height2 = computeHeight(nodes[node].right);

			return (int)(1 + Math.Max(height1, height2));
		}

		//! Get the height of the tree.
		/*! \return
		        The height of the binary tree.
		 */
		public int getHeight()
		{
			if (root == -1)
			{
				return 0;
			}
			return nodes[root].height;
		}

		//! Get the number of nodes in the tree.
		/*! \return
		        The number of nodes in the tree.
		 */
		public int getNodeCount()
		{
			return nodeCount;
		}

		//! Compute the maximum balancance of the tree.
		/*! \return
		        The maximum difference between the height of two
		        children of a node.
		 */
		public int computeMaximumBalance()
		{
			int maxBalance = 0;
			for (int i = 0; i < nodeCapacity; i++)
			{
				if (nodes[i].height <= 1)
				{
					continue;
				}

				Debug.Assert(nodes[i].isLeaf() == false);

				int balance = (int)Math.Abs(nodes[nodes[i].left].height - nodes[nodes[i].right].height);
				maxBalance = Math.Max(maxBalance, balance);
			}

			return maxBalance;
		}

		//! Compute the surface area ratio of the tree.
		/*! \return
		        The ratio of the sum of the node surface area to the surface
		        area of the root node.
		 */
		public double computeSurfaceAreaRatio()
		{
			if (root == -1)
			{
				return 0.0;
			}

			double rootArea = nodes[root].aabb.computeSurfaceArea();
			double totalArea = 0.0;

			for (int i = 0; i < nodeCapacity;i++)
			{
				if (nodes[i].height < 0)
				{
					continue;
				}

				totalArea += nodes[i].aabb.computeSurfaceArea();
			}

			return totalArea / rootArea;
		}

		/// Validate the tree.
		public void validate()
		{
            //#if ! NDEBUG
			validateStructure(root);
			validateMetrics(root);

			int freeCount = 0;
			int freeIndex = freeList;

			while (freeIndex != -1)
			{
				Debug.Assert(freeIndex < nodeCapacity);
				freeIndex = nodes[freeIndex].next;
				freeCount++;
			}

			Debug.Assert(getHeight() == computeHeight());
			Debug.Assert((nodeCount + freeCount) == nodeCapacity);
            //#endif
		}

		/// Rebuild an optimal tree.
		public void rebuild()
		{
			List<int> nodeIndices = _.List<int>(nodeCount);
			int count = 0;

			for (int i = 0;i < nodeCapacity;i++)
			{
				// Free node.
				if (nodes[i].height < 0)
				{
					continue;
				}

				if (nodes[i].isLeaf())
				{
					nodes[i].parent = -1;
					nodeIndices[count] = i;
					count++;
				}
				else
				{
					freeNode(i);
				}
			}

			while (count > 1)
			{
				double minCost = double.MaxValue;
				int iMin = -1;
				int jMin = -1;

				for (int i = 0;i < count;i++)
				{
					AABB aabbi = nodes[nodeIndices[i]].aabb;

					for (int j = i + 1;j < count;j++)
					{
						AABB aabbj = nodes[nodeIndices[j]].aabb;
						AABB aabb = AABB.Union(aabbi, aabbj);
						double cost = aabb.getSurfaceArea();

						if (cost < minCost)
						{
							iMin = (int)i;
							jMin = (int)j;
							minCost = cost;
						}
					}
				}

				int index1 = nodeIndices[iMin];
				int index2 = nodeIndices[jMin];

				int parent = allocateNode();
				nodes[parent].left = index1;
				nodes[parent].right = index2;
				nodes[parent].height = 1 + Math.Max(nodes[index1].height, nodes[index2].height);
				nodes[parent].aabb = AABB.Union(nodes[index1].aabb, nodes[index2].aabb);
				nodes[parent].parent = -1;

				nodes[index1].parent = parent;
				nodes[index2].parent = parent;

				nodeIndices[jMin] = nodeIndices[count - 1];
				nodeIndices[iMin] = parent;
				count--;
			}

			root = nodeIndices[0];

			validate();
		}


		//! Assert that the sub-tree has a valid structure.
		/*! \param node
		        The index of the root node.
		 */
		private void validateStructure(int node)
		{
			if (node == -1)
			{
				return;
			}

			if (node == root)
			{
				Debug.Assert(nodes[node].parent == -1);
			}

			int left = nodes[node].left;
			int right = nodes[node].right;

			if (nodes[node].isLeaf())
			{
				Debug.Assert(left == -1);
				Debug.Assert(right == -1);
				Debug.Assert(nodes[node].height == 0);
				return;
			}

			Debug.Assert(left < nodeCapacity);
			Debug.Assert(right < nodeCapacity);

			Debug.Assert(nodes[left].parent == node);
			Debug.Assert(nodes[right].parent == node);

			validateStructure(left);
			validateStructure(right);
		}

		//! Assert that the sub-tree has valid metrics.
		/*! \param node
		        The index of the root node.
		 */
		private void validateMetrics(int node)
		{
			if (node == -1)
			{
				return;
			}

			int left = nodes[node].left;
			int right = nodes[node].right;

			if (nodes[node].isLeaf())
			{
				Debug.Assert(left == -1);
				Debug.Assert(right == -1);
				Debug.Assert(nodes[node].height == 0);
				return;
			}

			Debug.Assert(left < nodeCapacity);
			Debug.Assert(right < nodeCapacity);

			int height1 = nodes[left].height;
			int height2 = nodes[right].height;
			int height = (int)(1 + Math.Max(height1, height2));
			//height; // Unused variable in Release build
			Debug.Assert(nodes[node].height == height);

			AABB aabb = AABB.Union(nodes[left].aabb, nodes[right].aabb);

			for (int i = 0;i < dimension;i++)
			{
				Debug.Assert(aabb.lowerBound[i] == nodes[node].aabb.lowerBound[i]);
				Debug.Assert(aabb.upperBound[i] == nodes[node].aabb.upperBound[i]);
			}

			validateMetrics(left);
			validateMetrics(right);
		}

		//! Apply periodic boundary conditions.
		/* \param position
		        The position vector.
		 */
		private void periodicBoundaries(List<double> position)
		{
			for (int i = 0;i < dimension;i++)
			{
				if (position[i] < 0)
				{
					position[i] += boxSize[i];
				}
				else
				{
					if (position[i] >= boxSize[i])
					{
						position[i] -= boxSize[i];
					}
				}
			}
		}

		//! Compute minimum image separation.
		/*! \param separation
		        The separation vector.

		    \param shift
		        The shift vector.

		    \return
		        Whether a periodic shift has been applied.
		 */
		private bool minimumImage(List<double> separation, List<double> shift)
		{
			bool isShifted = false;

			for (int i = 0;i < dimension;i++)
			{
				if (separation[i] < negMinImage[i])
				{
					separation[i] += Convert.ToDouble(periodicity[i]) * boxSize[i];
					shift[i] = Convert.ToDouble(periodicity[i]) * boxSize[i];
					isShifted = true;
				}
				else
				{
					if (separation[i] >= posMinImage[i])
					{
						separation[i] -= Convert.ToDouble(periodicity[i]) * boxSize[i];
						shift[i] = -Convert.ToDouble(periodicity[i]) * boxSize[i];
						isShifted = true;
					}
				}
			}

			return isShifted;
		}

		public override string ToString()
		{
			string s = "";
			foreach (Node node in nodes)
			{
				s += node.aabb + "\n";
			}
			return s;
		}

		void _Draw(Graphics g, int node, int depth, int scale)
		{
			if (node == -1) return;

			AABB aabb = nodes[node].aabb;
			bool isLeaf = nodes[node].isLeaf();
			int height = nodes[node].height;
			if (aabb.lowerBound.Count != 2) return;

			int x = Convert.ToInt32(scale * aabb.lowerBound[0]);
			int y = Convert.ToInt32(scale * aabb.lowerBound[1]);
			int w = Convert.ToInt32(scale * (aabb.upperBound[0] - aabb.lowerBound[0]));
			int h = Convert.ToInt32(scale * (aabb.upperBound[1] - aabb.lowerBound[1]));

			Color color = isLeaf? Color.Blue: Color.Black;
			int penWidth = 1;//isLeaf? 1: Math.Max(1, 5 - depth);
			Pen pen = new Pen(color, penWidth);
			int rgb = Math.Min(255, 128+depth*20);
			SolidBrush brush = new SolidBrush(Color.FromArgb(rgb, rgb, rgb));

			if (isLeaf)
			{
				g.DrawEllipse(pen, x, y, w, h);
			} else {
				//g.FillRectangle(brush, x, y, w, h);
				g.DrawRectangle(pen, x, y, w, h);
			}

			_Draw(g, nodes[node].right, depth + 1, scale);
			_Draw(g, nodes[node].left, depth + 1, scale);
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
