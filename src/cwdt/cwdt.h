#ifndef CWDT_H
#define CWDT_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <array>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include "gurobi_c++.h"
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT Weight;
typedef K::Point_3 Point;
typedef CGAL::Regular_triangulation_vertex_base_3<K> Vb0;
typedef CGAL::Regular_triangulation_cell_base_3<K> Cb0;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, Vb0> Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<int, K, Cb0> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<K, Tds> Rt;
typedef Rt::Vertex_iterator Vertex_iterator;
typedef Rt::Finite_vertices_iterator Finite_vertices_iterator;
typedef Rt::Finite_cells_iterator Finite_cells_iterator;
typedef Rt::Vertex_handle VH;
typedef Rt::Cell_handle CH;

namespace CWDT {
	class edge {
	public:
		edge(const int v1, const int v2) {
			ev_ = { v1, v2 };
		}

		bool operator==(edge const& e) const {
			return (ev_[0] == e.ev_[0] || ev_[0] == e.ev_[1])
				&& (ev_[1] == e.ev_[0] && ev_[1] == e.ev_[1]);
		}

		int& operator[](const int i) {
			return ev_[i];
		}

		const int& operator[](const int i) const {
			return ev_[i];
		}

	private:
		std::array<int, 2> ev_;
	};

	class triangle {
	public:
		triangle() {};
		triangle(const int v1, const int v2, const int v3) {
			tv_ = { v1, v2, v3 };
		}

		bool operator==(triangle const& t) const {
			return ((tv_[0] == t.tv_[0]) || (tv_[0] == t.tv_[1]) || (tv_[0] == t.tv_[2]))
				&& ((tv_[1] == t.tv_[0]) || (tv_[1] == t.tv_[1]) || (tv_[1] == t.tv_[2]))
				&& ((tv_[2] == t.tv_[0]) || (tv_[2] == t.tv_[1]) || (tv_[2] == t.tv_[2]));
		}

		int& operator[](int i) {
			return tv_[i];
		}

		const int& operator[](int i) const {
			return tv_[i];
		}
	private:
		std::array<int, 3> tv_;
	};

	class poly {
	public:
		poly() {};

		void insert_boundary_edge(
			const edge& e
			) {
			be_.push_back(e);
		}

		std::list<edge>& be() { return be_; }
		std::vector<int>& iv() { return iv_; }
		std::vector<triangle>& tri() { return tri_; }
		K::Plane_3& plane() { return plane_; }
	private:
		std::list<edge> be_;
		std::vector<int> iv_;
		std::vector<triangle> tri_;
		K::Plane_3 plane_;
	};
}

template<>
struct std::hash<CWDT::edge> {
	size_t operator()(CWDT::edge const& key) const {
		return static_cast<size_t>(key[0]) * static_cast<size_t>(key[1]);
	}
};

template<>
struct std::hash<CWDT::triangle> {
	size_t operator()(CWDT::triangle const& key) const {
		return static_cast<size_t>(key[0]) * static_cast<size_t>(key[1]) * static_cast<size_t>(key[2]);
	}
};

namespace CWDT {
	class processor {
	public:
		processor() {
			env_ = new GRBEnv();
			gabriel_ = false;
			minimize_heights_ = true;
		}

		int nv() const { return nv_; }
		int np() const { return np_; }
		GRBEnv* env() const { return env_; }
		Rt& T() { return T_; }

		int ntri();

		/** \brief Read in constrained mesh information.
		 * \param[in] file_path
		 * \param[in] in: ifstream of input file
		 * \return success or not */
		bool read_off(const std::string& file_path);
		bool read_igl(const std::string& file_path);
		bool write_off(const std::string& file_path);
		bool write_woff(const std::string& file_path);
		bool write_tet_off(const std::string& file_path,
			const double& shrink_factor = 1.);

		/** \brief Pre-compute the plane for each polygon. */
		void compute_planes();

		/** \brief Triangulate the specified polygon.
		 * \param[in, out] p: polygon
		 * \return the number of triangles after triangulation */
		int triangulate_poly(poly& p);
		int triangulate_cpoly_w(poly& p);
		int triangulate();

		/** \brief Get the adjacency relationship between tets, and store in tetNeighbors. */
		void buildTetNeighbors();

		int solve(int nr,
			double mintol, double maxtol,
			double height_factor, double weight_factor);

		int split_missing_edges();
		int split_missing_faces();

		void missing_faces(std::vector<int>& cmt);

		int peel_by_winding_number(double thr);

	private:
		int nv_; // vertex nb
		int np_; // polygon nb

		std::vector<Point> vertices_; /* constrained vertices */
		std::vector<poly> polygons_; /* constrained polygons */
		std::vector<double> weights_; /* vertex weights */

		std::unordered_map<edge, std::vector<int>> E2P_; /* map[edge] to the vertex indices that make up this edge */

		std::vector<std::array<int, 4>> tets_;
		std::vector<std::array<int, 4>> tetNeighbors_; /* the index of the adjacent tet corresponding to each vertex in tet,
														-1 means there is no adjacent tet */

		std::vector<VH> id2vh_;
		std::vector<bool> ve_;

		bool gabriel_;

		GRBEnv* env_;

		bool minimize_heights_, write_missing_, write_max_constraints_, write_witness_;

		Rt T_;
	};
}
#endif
