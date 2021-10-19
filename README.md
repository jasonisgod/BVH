# AABB Tree



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
- one class one file
- tree optimize 
	- (L3, L4, R1, R2, 13, 14)
	- list.Where(...).OrderBy(...).First();
- node-aabb-polygon
	- object <-> nodes pointer mapping
	- create aabb by polygon
- closest point search
	- dolfin BoundingBoxTree
- improve main run
- unity GUI
	- mouse event
	- insert/find/update/remove
