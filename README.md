# Conforming Weighted Delaunay Triangulations

This repo is a rearrangement of the code for the work ["**Conforming Weighted Delaunay Triangulations**"](https://dl.acm.org/doi/abs/10.1145/3414685.3417776) by [Marc Alexa](https://cg.tu-berlin.de/people/marc-alexa). [[Source Code]](https://cragl.cs.gmu.edu/iheartla/evaluation/static/cases_res/Conforming%20Weighted%20Delaunay%20Triangulations/cwdt3.cc)

I am attempting to make it more user-friendly and easier to use.

The code is now executable, but I have not yet finalized the CMakeLists. However, configuring this code is easy. You only need to properly configure the CGAL, Eigen, libigl, Gurobi, and CLI11 libraries in order to run this code.

If you use this code, you need to cite:

```
@article{10.1145/3414685.3417776,
	author = {Alexa, Marc},
	title = {Conforming weighted delaunay triangulations},
	year = {2020},
	issue_date = {December 2020},
	publisher = {Association for Computing Machinery},
	address = {New York, NY, USA},
	volume = {39},
	number = {6},
	issn = {0730-0301},
	url = {https://doi.org/10.1145/3414685.3417776},
	doi = {10.1145/3414685.3417776},
	abstract = {Given a set of points together with a set of simplices we show how to compute weights associated with the points such that the weighted Delaunay triangulation of the point set contains the simplices, if possible. For a given triangulated surface, this process provides a tetrahedral mesh conforming to the triangulation, i.e. solves the problem of meshing the triangulated surface without inserting additional vertices. The restriction to weighted Delaunay triangulations ensures that the orthogonal dual mesh is embedded, facilitating common geometry processing tasks.We show that the existence of a single simplex in a weighted Delaunay triangulation for given vertices amounts to a set of linear inequalities, one for each vertex. This means that the number of inequalities for a given triangle mesh is quadratic in the number of mesh elements, making the naive approach impractical. We devise an algorithm that incrementally selects a small subset of inequalities, repeatedly updating the weights, until the weighted Delaunay triangulation contains all constrained simplices or the problem becomes infeasible. Applying this algorithm to a range of triangle meshes commonly used graphics demonstrates that many of them admit a conforming weighted Delaunay triangulation, in contrast to conforming or constrained Delaunay that require additional vertices to split the input primitives.},
	journal = {ACM Trans. Graph.},
	month = nov,
	articleno = {248},
	numpages = {16},
	keywords = {simplicial meshes, conforming meshes}
}
```



