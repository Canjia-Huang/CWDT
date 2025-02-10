#include "cwdt.h"
#include "macro.h"

#include <string>
// libigl
#include <igl/read_triangle_mesh.h>
#include <igl/bfs_orient.h>
#include <igl/winding_number.h>
// CGAL
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb2;
#include <CGAL/Triangulation_face_base_with_info_2.h>
struct FaceInfo2 {
    FaceInfo2() {}
    int nesting_level;
    bool in_domain() {
        return nesting_level % 2 == 1;
    }
};
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Fbb;
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb> Fb2;
typedef CGAL::Triangulation_data_structure_2<Vb2, Fb2> Tds2;
typedef CGAL::No_constraint_intersection_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds2, Itag> CDT;

namespace CWDT {
    void mark_domains(CDT& ct,
        CDT::Face_handle start, int index,
        std::list<CDT::Edge>& border
    ) {
        if (start->info().nesting_level != -1) {
            return;
        }
        std::list<CDT::Face_handle> queue;
        queue.push_back(start);
        while (!queue.empty()) {
            CDT::Face_handle fh = queue.front();
            queue.pop_front();
            if (fh->info().nesting_level == -1) {
                fh->info().nesting_level = index;
                for (int i = 0; i < 3; i++) {
                    CDT::Edge e(fh, i);
                    CDT::Face_handle n = fh->neighbor(i);
                    if (n->info().nesting_level == -1) {
                        if (ct.is_constrained(e)) {
                            border.push_back(e);
                        }
                        else {
                            queue.push_back(n);
                        }
                    }
                }
            }
        }
    }
    void mark_domains(CDT& cdt
    ) {
        for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) {
            it->info().nesting_level = -1;
        }
        std::list<CDT::Edge> border;
        mark_domains(cdt, cdt.infinite_face(), 0, border);
        while (!border.empty()) {
            CDT::Edge e = border.front();
            border.pop_front();
            CDT::Face_handle n = e.first->neighbor(e.second);
            if (n->info().nesting_level == -1) {
                mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
            }
        }
    }


    bool processor::read_OFF(std::ifstream& in
    ) {
        if (!in.good()) {
            WARNING("Input file is not good!");
            return false;
        }

        std::string line;
        std::getline(in, line);

        std::istringstream isf(line);

        /** read the first line **/
        std::string ft;
        isf >> ft;

        if (ft.size() < 3 || (ft.substr(ft.size() - 3)) != "OFF") return false;
        bool woff = (ft == "WOFF"); // with weight?

        /** read the second line **/
        int ne;
        std::getline(in, line);
        std::istringstream iss(line);
        iss >> nv_ >> np_ >> ne;

        /** read vertices (and weights) **/
        vertices_.resize(nv_);
        weights_.resize(nv_);

        for (int i = 0; i < nv_; ++i) {
            std::getline(in, line);
            std::istringstream iss(line);
            double x, y, z, w = 0.0;
            if (woff) iss >> x >> y >> z >> w;
            else iss >> x >> y >> z;
            vertices_[i] = Point(x, y, z);
            weights_[i] = w;
        }

        /** read polygons **/
        polygons_.resize(np_);

        double mina = 180.0;
        int nf = 0;
        for (int i = 0; i < np_; i++) {
            std::getline(in, line);
            std::istringstream iss(line);
            int d;
            iss >> d;
            if (d < 3) WARNING("Degree < 3");

            int* id = (int*)malloc(sizeof(int) * d);

            for (int v = 0; v < d; ++v) iss >> id[v];
            int prev_v = d - 1;
            for (int v = 0; v < d; ++v) {
                edge e(id[prev_v], id[v]);

                polygons_[i].insert_boundary_edge(e);
                E2P_[e].push_back(i);

                prev_v = v;
            }
            for (int v = 0; v < d; ++v) {
                double a = CGAL::approximate_angle(vertices_[id[v]], vertices_[id[(v + 1) % d]], vertices_[id[(v + 2) % d]]);
                if (a < mina) mina = a;
            }

            free(id);
        }

        /*std::vector<int> degree(nv_, 0);
        for (const auto& ep : E2P_) {
            ++degree[ep.first[0]];
            ++degree[ep.first[1]];
        }
        for (int i = 0; i < nv_; i++) {
            if (degree[i] < 3) WARNING("Degree of vertex " << i << " is " << degree[i]);
        }*/

        return true;
    }

    void processor::compute_planes(
    ) {
        for (int i = 0; i < np_; ++i) {
            std::vector<Point> v;
            for (const auto& e : polygons_[i].be()) {
                v.push_back(vertices_[e[0]]);
            }
                
            if (v.size() == 3) {
                polygons_[i].plane() = K::Plane_3(v[0], v[1], v[2]);
            }
            else {
                CGAL::linear_least_squares_fitting_3(v.begin(), v.end(), polygons_[i].plane(), CGAL::Dimension_tag<0>());
            }
        }
    }

    int processor::triangulate_poly(poly& p
    ) {
        int nf = p.be().size() + 2 * p.iv().size() - 2;
        p.tri().resize(nf);

        if (nf == 1) {
            auto eit = p.be().begin();
            for (int v = 0; v < 3; v++, ++eit) {
                p.tri()[0][v] = (*eit)[0];
            }
            return 1;
        }

        CDT cdt;
        std::vector<std::pair<CDT::Point, int>> points;
        for (const auto& e : p.be()) {
            points.push_back(std::make_pair(p.plane().to_2d(vertices_[e[0]]), e[0]));
        }
        if (!p.iv().empty()) {
            for (const auto& v : p.iv()) {
                points.push_back(std::make_pair(p.plane().to_2d(vertices_[v]), v));
            }
        }

        cdt.insert(points.begin(), points.end());
        std::unordered_map<int, CDT::Vertex_handle> vhm;
        CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
        for (; vit != cdt.finite_vertices_end(); ++vit) {
            vhm[vit->info()] = vit;
        }
        for (const auto& e : p.be()) {
            cdt.insert_constraint(vhm[e[0]], vhm[e[1]]);
        }
        mark_domains(cdt);
        CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
        int fc = 0;
        for (; fit != cdt.finite_faces_end(); ++fit) {
            if (fit->info().in_domain()) {
                for (int v = 0; v < 3; v++) {
                    p.tri()[fc][v] = fit->vertex(v)->info();
                }
                fc++;
            }
        }

        return fc;
    }

	void processor::buildTetNeighbors(
	) {
		VERBOSE("Building tet neighbors");

		const int NT = tets_.size();

        /** finding vertex-adjacent tets **/
        std::vector<std::vector<int>> vertex_adjacent_tets(nv_); // tet indices adjacent to each vertex

        for (int i = 0; i < NT; ++i) {
            for (int j = 0; j < 4; ++j) vertex_adjacent_tets[tets_[i][j]].push_back(i);
        }

        /** finding tet-adjacent tets **/
        tetNeighbors_.resize(NT);

        for (int i = 0; i < NT; ++i) {
            for (int j = 0; j < 4; ++j) {
                // find tets opposite to vertex j

                triangle cft(tets_[i][j % 4], tets_[i][(j + 1) % 4], tets_[i][(j + 2) % 4]); // a facet of tet[i]

                int bc; // the vertex index with the least number of adjacent tets
                bc = cft[0];
                if (vertex_adjacent_tets[cft[1]].size() < vertex_adjacent_tets[bc].size()) {
                    bc = cft[1];
                }
                if (vertex_adjacent_tets[cft[2]].size() < vertex_adjacent_tets[bc].size()) {
                    bc = cft[2];
                }

                int curr = -1;
                for (const auto& at : vertex_adjacent_tets[bc]) { // vertex[bc]'s neighboring tet index
                    if (at != i) {
                        for (int k = 0; k < 4; k++) {
                            triangle aft(tets_[at][k % 4], tets_[at][(k + 1) % 4], tets_[at][(k + 2) % 4]);
                            if (aft == cft) {
                                curr = at;
                                break;
                            }
                        }
                    }
                }

                tetNeighbors_[i][j] = curr;
            }
        }

        VERBOSE("Done!");
	}
}
