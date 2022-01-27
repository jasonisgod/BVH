# Implementation of BVH using AABBs

![GitHub last commit](https://img.shields.io/github/last-commit/jasonisgod/BVH) 
![GitHub commit checks state](https://img.shields.io/github/checks-status/jasonisgod/BVH/master) 
![GitHub](https://img.shields.io/github/license/jasonisgod/BVH) 

## About
A C# implementation of a dynamic bounding volume hierarchy
([BVH](https://en.wikipedia.org/wiki/Bounding_volume_hierarchy)) using
axis-aligned bounding boxes ([AABBs](https://en.wikipedia.org/wiki/Minimum_bounding_box)).
The data structure provides an efficient way of detecting potential overlap
between objects of arbitrary shape and size and is commonly used in
computer game engines for collision detection and ray tracing.

Because of their speed and flexibility, AABB trees are also well suited
to overlap detection in physics applications, such as molecular simulation.
They are particularly helpful for systems where there is a large size disparity
between particle species, or whenever the particle density is extremely
inhomogeneous.

![Alt text](README.png?raw=true "Title")

## Optimization


### Priority (high to low)
- Tree Balance
	- Minimize the difference height (at most 1) between left and right for every nodes
- Intersect area
	- Minimize the intersect area of left and right for every nodes
- Union area
	- Minimize the union area of left and right for every nodes
- Difference of area
	- Minimize the difference of area between left and right for every nodes

### Operation

```
           0
    L              R
LL     LR      RL      RR
```
- L-RL-RR
	- Swap(L, RL)
	- Swap(L, RR)
- R-LL-LR
	- Swap(R, LL)
	- Swap(R, LR)
- LL-LR-RL-RR
	- Swap(LL, RL)
	- Swap(LL, RR)

