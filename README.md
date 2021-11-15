# Implementation of BVH using AABBs

![GitHub last commit](https://img.shields.io/github/last-commit/jasonisgod/BVH) 
![GitHub commit checks state](https://img.shields.io/github/checks-status/jasonisgod/BVH/master) 

![CRAN/METACRAN](https://img.shields.io/cran/l/devtools) 
![NodePing uptime](https://img.shields.io/nodeping/uptime/jkiwn052-ntpp-4lbb-8d45-ihew6d9ucoei)

![PingPong status](https://img.shields.io/pingpong/status/sp_2e80bc00b6054faeb2b87e2464be337e) 
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

![Alt text](img/test.png?raw=true "Title")

## Optimization
```
			0
	L				R
1		2		3		4
```

Priority
1. Tree Balance
2. Min(Overlap area)
3. Min(Total area)
4. Min(Diff of pairs)

Options
- Swap(L, 3)
- Swap(L, 4)
- Swap(R, 1)
- Swap(R, 2)
- Swap(1, 3)
- Swap(1, 4)

## TODO
- node-aabb-polygon
	- capitalize method name
	- edit insertLeaf(): find nearest as sibling
	- object <-> nodes pointer mapping
	- create aabb by polygon
- closest point search
	- dolfin BoundingBoxTree
- unity GUI
	- insert/find/update/remove
