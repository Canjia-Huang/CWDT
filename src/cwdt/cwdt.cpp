#include "cwdt/cwdt.h"

#include <fstream>
#include <string>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <igl/bfs_orient.h>
#include <igl/read_triangle_mesh.h>
#include <igl/winding_number.h>

#ifdef CWDT_VERBOSE
#	define VERBOSE(x) std::cout << "\033[33m" << "[" << __FUNCTION__ << "]" << "\033[0m" << " " << x << std::endl // [yellow] + white cout
#	define WARNING(x) std::cerr << "\033[33m" << "[" << __FILE__ << " " << __LINE__ << "]" << "\033[0m" << " " << "\033[31m" << x << "\033[0m" << std::endl // [yellow] + red cout
#else
#	define VERBOSE(x)
#	define WARNING(x)
#endif

#ifdef CWDT_DEBUG
#endif

// CGAL
typedef K::Vector_3 Vector;
typedef K::Weighted_point_3 Weighted_point;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb2;
struct FaceInfo2 {
    FaceInfo2(): nesting_level(0) {}

    int nesting_level;
    bool in_domain() const {
        return nesting_level % 2 == 1;
    }
};
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb> Fb2;
typedef CGAL::Triangulation_data_structure_2<Vb2, Fb2> Tds2;
typedef CGAL::No_constraint_intersection_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds2, Itag> CDT;
typedef CGAL::Regular_triangulation_vertex_base_2<K> Vbr0;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K, Vbr0> Vb2r;
typedef CGAL::Regular_triangulation_face_base_2<K> Cbr0;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K, Cbr0> Cb2r;
typedef CGAL::Triangulation_data_structure_2<Vb2r, Cb2r> Tds2r;

typedef CGAL::Regular_triangulation_2<K, Tds2r> Rt2;

namespace CWDT {
    void mark_domains(
        const CDT& ct,
        const CDT::Face_handle start, const int index,
        std::list<CDT::Edge>& border
        ) {
        if (start->info().nesting_level != -1)
            return;

        std::list<CDT::Face_handle> queue;
        queue.push_back(start);
        while (!queue.empty()) {
            CDT::Face_handle fh = queue.front();
            queue.pop_front();
            if (fh->info().nesting_level == -1) {
                fh->info().nesting_level = index;
                for (int i = 0; i < 3; i++) {
                    CDT::Edge e(fh, i);
                    if (CDT::Face_handle n = fh->neighbor(i);
                        n->info().nesting_level == -1) {
                        if (ct.is_constrained(e))
                            border.push_back(e);
                        else
                            queue.push_back(n);
                    }
                }
            }
        }
    }

    void mark_domains(
        const CDT& cdt
        ) {
            for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) {
                it->info().nesting_level = -1;
            }
            std::list<CDT::Edge> border;
            mark_domains(cdt, cdt.infinite_face(), 0, border);
            while (!border.empty()) {
                const CDT::Edge e = border.front();
                border.pop_front();
                if (const CDT::Face_handle n = e.first->neighbor(e.second);
                    n->info().nesting_level == -1)
                    mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
            }
    }

    int processor::ntri(
        ) {
        int nt = 0;
        for (auto& p : polygons_)
            nt += p.tri().size();
        return nt;
    }

    void processor::compute_planes(
        ) {
        for (int i = 0; i < np_; ++i) {
            std::vector<Point> v;
            for (const auto& e : polygons_[i].be())
                v.push_back(vertices_[e[0]]);

            if (v.size() == 3)
                polygons_[i].plane() = K::Plane_3(v[0], v[1], v[2]);
            else
                CGAL::linear_least_squares_fitting_3(v.begin(), v.end(), polygons_[i].plane(), CGAL::Dimension_tag<0>());
        }
    }

    int processor::triangulate_poly(
        poly& p
        ) {
        const int nf = p.be().size() + 2 * p.iv().size() - 2;
        p.tri().resize(nf);

        if (nf == 1) {
            auto eit = p.be().begin();
            for (int v = 0; v < 3; v++, ++eit)
                p.tri()[0][v] = (*eit)[0];
            return 1;
        }

        CDT cdt;
        std::vector<std::pair<CDT::Point, int>> points;
        for (const auto& e : p.be())
            points.push_back(std::make_pair(p.plane().to_2d(vertices_[e[0]]), e[0]));
        if (!p.iv().empty()) {
            for (const auto& v : p.iv())
                points.push_back(std::make_pair(p.plane().to_2d(vertices_[v]), v));
        }

        cdt.insert(points.begin(), points.end());
        std::unordered_map<int, CDT::Vertex_handle> vhm;
        CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
        for (; vit != cdt.finite_vertices_end(); ++vit)
            vhm[vit->info()] = vit;
        for (const auto& e : p.be())
            cdt.insert_constraint(vhm[e[0]], vhm[e[1]]);
        mark_domains(cdt);
        CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
        int fc = 0;
        for (; fit != cdt.finite_faces_end(); ++fit) {
            if (fit->info().in_domain()) {
                for (int v = 0; v < 3; v++)
                    p.tri()[fc][v] = fit->vertex(v)->info();
                fc++;
            }
        }

        return fc;
    }

    int processor::triangulate_cpoly_w(
        poly& p
        ) {
        const int nf = p.be().size() + 2 * p.iv().size() - 2;
        p.tri().resize(nf);
        if (nf == 1) { // it's a triangle, so just create the single triangle from the edges
            auto eit = p.be().begin();
            for (int v = 0; v < 3; v++, ++eit)
                p.tri()[0][v] = (*eit)[0];
            return 1;
        }

        Rt2 rt2;
        std::vector<std::pair<Rt2::Point, int>> points;
        for (const auto& e : p.be()) {
            points.push_back(std::make_pair(
                Rt2::Point(p.plane().to_2d(vertices_[e[0]]), weights_[e[0]]),
                e[0]));
        }
        if (!p.iv().empty()) {
            for (const auto& v : p.iv())
                points.push_back(std::make_pair(
                    Rt2::Point(p.plane().to_2d(vertices_[v]), weights_[v]),
                    v));
        }
        rt2.insert(points.begin(), points.end());
        std::unordered_map<int, Rt2::Vertex_handle> vhm;

        for (Rt2::Finite_vertices_iterator vit = rt2.finite_vertices_begin(); vit != rt2.finite_vertices_end(); ++vit)
            vhm[vit->info()] = vit;

        int fc = 0;
        for (Rt2::Finite_faces_iterator fit = rt2.finite_faces_begin(); fit != rt2.finite_faces_end(); ++fit) {
            for (int v = 0; v < 3; ++v)
                p.tri()[fc][v] = fit->vertex(v)->info();
            fc++;
        }
        if (fc != nf) {
            WARNING("Something very wrong?" << fc << " = " << nf);

            exit(1);
        }

        return fc;
    }

    int processor::triangulate(
        ) {
        int nf = 0;
        for (auto& p : polygons_)
            nf += triangulate_poly(p);
        return nf;
    }

	void processor::buildTetNeighbors(
	    ) {
		VERBOSE("Building tet neighbors");

		const int NT = tets_.size();

        /** finding vertex-adjacent tets **/
        std::vector<std::vector<int>> vertex_adjacent_tets(nv_); // tet indices adjacent to each vertex

        for (int i = 0; i < NT; ++i)
            for (int j = 0; j < 4; ++j)
                vertex_adjacent_tets[tets_[i][j]].push_back(i);

        /** finding tet-adjacent tets **/
        tetNeighbors_.resize(NT);

        for (int i = 0; i < NT; ++i) {
            for (int j = 0; j < 4; ++j) {
                // find tets opposite to vertex j

                triangle cft(tets_[i][j % 4], tets_[i][(j + 1) % 4], tets_[i][(j + 2) % 4]); // a facet of tet[i]

                int bc = cft[0]; // the vertex index with the least number of adjacent tets
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
                            if (triangle aft(tets_[at][k % 4], tets_[at][(k + 1) % 4], tets_[at][(k + 2) % 4]);
                                aft == cft) {
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

    int processor::solve(
        int nr,
        double mintol, double maxtol,
        double height_factor, double weight_factor
        ) {
        int nic = 0;

        std::vector<VH>(nv_).swap(id2vh_);

        std::vector<std::pair<Weighted_point, int>> Pi(nv_);

        int nf = ntri();
        VERBOSE("Number of constrained triangles: " << nf);

        std::vector<bool> dve(nf, false), ignore(nf, false);
        std::vector<std::unordered_map<int, double>> constr_tol(nf);
        std::vector<std::unordered_set<int>> constr(nf);

        Eigen::VectorXd H;
        Eigen::MatrixXd G;
        if (!gabriel_) {
            H = Eigen::VectorXd::Zero(nf);
            G = Eigen::MatrixXd::Zero(nv_, 3);
        }

        try {
            GRBModel model = GRBModel(*env_);
            model.set(GRB_IntParam_OutputFlag, 0);

            /** Create variables **/
            GRBVar* w = model.addVars(NULL, NULL, NULL, NULL, NULL, nv_);

            for (int i = 0; i < nv_; ++i) {
                if (weights_[i] > 0.0)
                    w[i].set(GRB_DoubleAttr_Start, weights_[i]);
            }

            GRBVar* h, * g;
            if (!gabriel_) {
                std::vector<double> lbf(nf, -GRB_INFINITY);
                std::vector<double> lbv(3 * nv_, -GRB_INFINITY);
                h = model.addVars(&lbf[0], NULL, NULL, NULL, NULL, nf);
                g = model.addVars(&lbv[0], NULL, NULL, NULL, NULL, 3 * nv_);
            }

            /** Set objective **/
            GRBQuadExpr qe;
            GRBLinExpr le;
            std::vector<double> onev(nv_, weight_factor); //,b(nv);
            std::vector<double> onep(nf, height_factor); //,b(nv);
            // squared norm of height vector
            if (!gabriel_)
                qe.addTerms(&onep[0], h, h, nf);
            // squared norm of weight vector
            qe.addTerms(&onev[0], w, w, nv_);
            le.addTerms(&onev[0], w, nv_);

            if (minimize_heights_)
                model.setObjective(qe, GRB_MINIMIZE);
            else
                model.setObjective(le, GRB_MINIMIZE);
            int ntc = 0;
            int nutc = 0;

            Eigen::MatrixXd Z(nf, 3);
            for (int r = 0; r < nr; ++r) {
                VERBOSE("Round " << r);

                for (int i = 0; i < nv_; i++)
                    Pi[i] = std::make_pair(Weighted_point(vertices_[i], weights_[i]), i);

                T_.clear();
                T_.insert(Pi.begin(), Pi.end());

                // set up integer ids for vertices and mark vertices in the triangulation
                std::vector<bool>(nv_, false).swap(ve_);
                for (Finite_vertices_iterator vit = T_.finite_vertices_begin(); vit != T_.finite_vertices_end(); ++vit) {
                    id2vh_[vit->info()] = vit;
                    ve_[vit->info()] = true;
                }

                int nacv = 0;
                int nmv = 0;
                int nmt = 0;
                int ti = 0;

                nmv = nv_ - T_.number_of_vertices();
                VERBOSE("Number of missing vertices: " << nmv << " / " << nv_);
                if (nmv < 0) {
                    for (int i = 0; i < nv_; ++i)
                        if (!ve_[i]) {
                            Eigen::RowVector3d p(
                                vertices_[i][0],
                                vertices_[i][1],
                                vertices_[i][2]);

                            Point z(G(i, 0) + vertices_[i][0], G(i, 1) + vertices_[i][1], G(i, 2) + vertices_[i][2]);
                            VH vh = T_.nearest_power_vertex(z);
                            int vi = vh->info();
                            if (vi == i)
                                VERBOSE("Weird!");
                            VERBOSE("Adding " << i << " - " << vi);

                            Eigen::RowVector3d pi(
                                vertices_[vi][0],
                                vertices_[vi][1],
                                vertices_[vi][2]);

                            Eigen::Vector3d a = 2.0 * (pi - p);
                            double rhs = pi.squaredNorm() - p.squaredNorm() - p * a - 1e-6;

                            GRBLinExpr l = 0;
                            if (!gabriel_) {
                                l += a[0] * g[i];
                                l += a[1] * g[i + nv_];
                                l += a[2] * g[i + nv_ + nv_];
                            }
                            l -= w[i];
                            l += w[vi];
                            model.addConstr(l, GRB_LESS_EQUAL, rhs, std::to_string(-i));
                            ntc++;
                        }
                }

                std::vector<int> ac(np_, 0);
                for (int i = 0; i < np_; ++i) {
                    for (int ii = 0; ii < polygons_[i].tri().size(); ii++, ti++) {
                        if (!ignore[ti]) {
                            int v0 = polygons_[i].tri()[ii][0];
                            int v1 = polygons_[i].tri()[ii][1];
                            int v2 = polygons_[i].tri()[ii][2];

                            CH ch;
                            int c0, c1, c2, c3 = -1;
                            if (ve_[v0] && ve_[v1] && ve_[v2] &&
                                T_.is_facet(
                                    id2vh_[v0], id2vh_[v1], id2vh_[v2],
                                    ch, c0, c1, c2)) {
                                c3 = 6 - c0 - c1 - c2;
                                if (!gabriel_ && !minimize_heights_)
                                    continue;
                                if (T_.is_Gabriel(ch, c3))
                                    continue;
                            }

                            Eigen::RowVector3d p0(
                                vertices_[v0][0],
                                vertices_[v0][1],
                                vertices_[v0][2]);
                            Eigen::RowVector3d p1(
                                vertices_[v1][0],
                                vertices_[v1][1],
                                vertices_[v1][2]);
                            Eigen::RowVector3d p2(
                                vertices_[v2][0],
                                vertices_[v2][1],
                                vertices_[v2][2]);

                            Eigen::Vector3d n(
                                polygons_[i].plane().orthogonal_vector()[0],
                                polygons_[i].plane().orthogonal_vector()[1],
                                polygons_[i].plane().orthogonal_vector()[2]);

                            double x0s = p0.squaredNorm();
                            double x1s = p1.squaredNorm();
                            double x2s = p2.squaredNorm();

                            Eigen::Matrix3d M;
                            M.row(0) = 2.0 * (p1 - p0);
                            M.row(1) = 2.0 * (p2 - p0);
                            M.row(2) = n;

                            Eigen::Matrix3d Mi = M.inverse();

                            Eigen::Vector3d r(
                                x1s - x0s,
                                x2s - x0s,
                                p0 * n);

                            Eigen::Vector3d m = Mi * r;

                            int vi;
                            double tol = mintol;
                            if (c3 < 0 || gabriel_) {
                                Eigen::Vector3d b(
                                    weights_[v0] - weights_[v1],
                                    weights_[v0] - weights_[v2],
                                    0);

                                if (!gabriel_)
                                    b[2] = H[ti];

                                Z.row(ti) = m + Mi * b;
                                Point p(Z(ti, 0), Z(ti, 1), Z(ti, 2));

                                // because of what's going on behind the scenes in CGAL
                                // the following is highly inefficient and should be optimized

                                vi = T_.nearest_power_vertex(p)->info();

                                auto vit = constr_tol[ti].find(vi);
                                if (vit != constr_tol[ti].end()) {
                                    tol = vit->second;
                                    tol *= 2.0;
                                    if (tol > maxtol) {
                                        ignore[ti] = true;
                                        nic++;
                                        VERBOSE("Ignoring " << ti << " because it has reach max tolerance");
                                        GRBConstr* c = model.getConstrs();
                                        for (int ci = 0; ci < model.get(GRB_IntAttr_NumConstrs); ++ci) {
                                            int ctid = atoi(c[ci].get(GRB_StringAttr_ConstrName).c_str());
                                            if (ctid == ti)
                                                model.remove(c[ci]);
                                        }
                                        model.update();
                                        continue;
                                    }
                                    vit->second = tol;
                                }
                                else
                                    constr_tol[ti][vi] = mintol;
                            }
                            else {
                                VH vh = ch->vertex(c3);
                                if (T_.is_infinite(vh) || constr[ti].count(vh->info()) > 0)
                                    vh = T_.mirror_vertex(ch, c3);
                                if (T_.is_infinite(vh) || constr[ti].count(vh->info()) > 0)
                                    continue;
                                vi = vh->info();
                                if (constr_tol[ti].find(vi) != constr_tol[ti].end())
                                    continue;
                                constr_tol[ti][vi] = mintol;
                            }

                            Eigen::RowVector3d pi(
                                vertices_[vi][0],
                                vertices_[vi][1],
                                vertices_[vi][2]);

                            Eigen::RowVector3d a = 2.0 * (pi - p0);
                            double rhs = -a * m + pi.squaredNorm() - x0s - tol;

                            a *= Mi;
                            GRBLinExpr l = 0;
                            l += (a[0] + a[1] - 1.0) * w[v0];
                            l -= a[0] * w[v1];
                            l -= a[1] * w[v2];
                            if (!gabriel_)
                                l += a[2] * h[ti];
                            l += w[vi];
                            model.addConstr(l, GRB_LESS_EQUAL, rhs, std::to_string(ti));
                            constr[ti].insert(vi);
                            ntc++;
                            nmt++;
                            ac[i] = 1;
                            if (c3 < 0) ac[i] = 2;
                        }
                    }
                }

                if (nmt == 0) {
                    VERBOSE("No constraint violations anymore, breaking in round: " << r);
                    int mnc = 0, mnci = -1, tncv = 0;
                    std::unordered_map<int, int> chist;
                    for (int csi = 0; csi < ti; csi++) {
                        int cts = constr_tol[csi].size();
                        if (cts > mnc) {
                            mnci = csi;
                            mnc = constr_tol[csi].size();
                        }
                        tncv += cts;
                        chist[cts]++;
                    }
                    VERBOSE( "Constraint triangle " << mnci << " has max # of constraints: " << mnc);
                    VERBOSE("Total constraints: " << tncv << ", per constrained triangle: " << ((double)tncv / (double)ntri()));
                    VERBOSE("Histogram");
                    for (int csi = 0; csi <= mnc; csi++)
                        VERBOSE(csi << ": " << chist[csi]);
                    return nic;
                }
                else
                    VERBOSE("Number of missing triangles: " << nmt << " / " << ti << std::endl
                            << "Total number of constraints added: " << ntc);

                model.optimize();

                int status = model.get(GRB_IntAttr_Status);
                VERBOSE("Status: ");

                if (status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD) {
                    model.computeIIS();
                    VERBOSE("\nThe following constraint cannot be satisfied: ");
                    GRBConstr* c = model.getConstrs();
                    int tid = -1;
                    for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i) {
                        if (c[i].get(GRB_IntAttr_IISConstr) == 1) {
                            tid = atoi(c[i].get(GRB_StringAttr_ConstrName).c_str());
                            std::cerr << tid << std::endl;
                            if (tid >= 0)
                                break;
                        }
                    }
                    if (ignore[tid])
                        std::cerr << "Already ignoring tid: " << tid << std::endl;
                    ignore[tid] = true;
                    nic++;
                    for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i) {
                        int ctid = atoi(c[i].get(GRB_StringAttr_ConstrName).c_str());
                        if (ctid == tid)
                            model.remove(c[i]);
                    }
                    model.update();
                }
                else {
                    for (int i = 0; i < nv_; i++)
                        weights_[i] = w[i].get(GRB_DoubleAttr_X);
                    if (!gabriel_) {
                        for (int i = 0; i < nf; ++i)
                            H[i] = h[i].get(GRB_DoubleAttr_X);
                        for (int i = 0; i < nv_; ++i) {
                            G(i, 0) = g[i].get(GRB_DoubleAttr_X);
                            G(i, 1) = g[i + nv_].get(GRB_DoubleAttr_X);
                            G(i, 2) = g[i + nv_ + nv_].get(GRB_DoubleAttr_X);
                        }
                    }
                }
            }
        }
        catch (GRBException e) {
            WARNING("Error code = " << e.getErrorCode());
            WARNING(e.getMessage());
        }

        std::vector<double>(nv_, 0.).swap(weights_);

        return -1;
    }

    int processor::split_missing_edges(
        ) {
        VERBOSE("Checking edges... ");
        std::vector<std::pair<edge, std::vector<int>>> me;
        for (const auto& it : E2P_) {
            edge e = it.first;
            CH ch;
            int c0, c1;
            if (!ve_[e[0]] || !ve_[e[1]] ||
                !T_.is_edge(id2vh_[e[0]], id2vh_[e[1]], ch, c0, c1))
                me.emplace_back(it);
        }
        const int nme = me.size();
        VERBOSE("Edges missing (encroached): " << nme);
        if (nme == 0)
            return 0;

        std::set<int> mp;
        for (int i = 0; i < nme; ++i) {
            edge e = me[i].first;
            Point em = CGAL::midpoint(vertices_[e[0]], vertices_[e[1]]);
            vertices_.push_back(em);
            weights_.push_back(0.0);

            for (const auto& pid : me[i].second) {
                auto eit = std::find(polygons_[pid].be().begin(), polygons_[pid].be().end(), e);
                if (eit == polygons_[pid].be().end())
                    VERBOSE("Edge " << e[0] << ", " << e[1] << " not found in poly " << pid);

                const int v = (*eit)[0];
                (*eit)[0] = nv_;
                polygons_[pid].be().insert(eit, edge{v, nv_});

                E2P_[edge{ e[0], nv_ }].push_back(pid);
                E2P_[edge{ e[1], nv_ }].push_back(pid);
                mp.insert(pid);
            }
            E2P_.erase(e);
            ++nv_;
        }
        VERBOSE("Modified " << mp.size());
        std::vector<bool> af(np_, true);

        for (const auto& i : mp) {
            triangulate_poly(polygons_[i]);
            af[i] = false;
        }

        return nme;
    }

    int processor::split_missing_faces(
        ) {
        std::unordered_set<int> mp;
        std::unordered_set<edge> ebs;
        std::vector<std::pair<Point, int>> ccp;
        int ni = 0;

        for (int i = 0; i < np_; ++i) {
            std::unordered_set<edge> eis;

            for (int ii = 0; ii < polygons_[i].tri().size(); ++ii) {
                triangle t = polygons_[i].tri()[ii];
                VERBOSE("Checking triangle " << t[0] << ", " << t[1] << ", " << t[2]);
                // we can assume all vertices are there, because all edges are there
                CH ch;
                int c0, c1, c2;

                if (!T_.is_facet(
                    id2vh_[t[0]], id2vh_[t[1]], id2vh_[t[2]],
                    ch, c0, c1, c2)) {
                    bool icc = true;
                    Point cc = CGAL::circumcenter(vertices_[t[0]], vertices_[t[1]], vertices_[t[2]]);
                    VERBOSE("cc: " << cc);

                    for (int j = 0; j < 3; ++j) {
                        // consider edge j, j+2
                        edge e = edge{ t[j], t[(j + 2) % 3] };
                        // is it part of the boundary of P[i]?
                        auto eit = std::find(polygons_[i].be().begin(), polygons_[i].be().end(), e);
                        bool ob = (eit != polygons_[i].be().end());
                        // is the angle opposite e obtuse?
                        bool wa = (CGAL::angle(vertices_[t[j]], vertices_[t[(j + 1) % 3]], vertices_[t[(j + 2) % 3]]) != CGAL::ACUTE);
                        // if angle is obtuse split edge in any case
                        // if edge is on boundary and encroached by cc also split
                        if (wa) {
                            if (ob) {
                                ebs.insert(e);
                            }
                            else {
                               VERBOSE("Should insert on interior edge " << e[0] << " - " << e[1] << " in poly " << i);
                                const auto& ieit = eis.find(e);
                                if (ieit == eis.end()) {
                                   VERBOSE("Not found for polygon " << i << "  -> Inserting");
                                    Point em = CGAL::midpoint(vertices_[e[0]], vertices_[e[1]]);
                                    vertices_.push_back(em);
                                    weights_.push_back(0.0);
                                    polygons_[i].iv().push_back(nv_);
                                    mp.insert(i);
                                    ++nv_;
                                    ++ni;
                                    eis.insert(e);
                                }
                            }
                            icc = false;
                            break;
                        }

                        if (ob) {
                            Point em = CGAL::midpoint(vertices_[e[0]], vertices_[e[1]]);
                            if ((cc - em).squared_length() <= (vertices_[e[0]] - em).squared_length()) {
                                ebs.insert(e);
                                icc = false;
                                break;
                            }
                        }
                    }

                    if (icc) {
                       VERBOSE("Potentially inserting cc");
                       ccp.emplace_back(cc, i);
                    }
                }
            }
        }

       VERBOSE("Interior insertions: " << ni);
       VERBOSE("Encroached edges by ccs: " << (int)ebs.size());

       for (const auto& e : ebs) {
           Point em = CGAL::midpoint(vertices_[e[0]], vertices_[e[1]]);
           vertices_.push_back(em);
           weights_.push_back(0.0);

           for (auto pid : E2P_[e])  {
               auto eit = std::find(polygons_[pid].be().begin(), polygons_[pid].be().end(), e);
               if (eit == polygons_[pid].be().end())
                   VERBOSE("Edge " << e[0] << ", " << e[1] << " not found in poly " << pid);

               int v = (*eit)[0];
               (*eit)[0] = nv_;
               polygons_[pid].be().insert(eit, edge{v, nv_});

               E2P_[edge{ e[0], nv_ }].push_back(pid);
               E2P_[edge{ e[1], nv_ }].push_back(pid);
               mp.insert(pid);
           }
           E2P_.erase(e);
           ++nv_;
       }

       VERBOSE("Modified " << mp.size() << "by edge splitting");
       std::vector<bool> af(np_, true);

       if (mp.size() == 0) {
           for (const auto& cp : ccp) {
               vertices_.push_back(cp.first);
               weights_.push_back(0.0);
               polygons_[cp.second].iv().push_back(nv_);
               mp.insert(cp.second);
               ++nv_;
           }
       }

       for (const auto& i : mp) {
           triangulate_poly(polygons_[i]);
           af[i] = false;
       }

       return ni + ebs.size() + ccp.size();
    }

    void processor::missing_faces(
        std::vector<int>& cmt
        ) {
        std::vector<int>().swap(cmt);
        for (int i = 0; i < np_; ++i) {
            for (int ii = 0; ii < polygons_[i].tri().size(); ++ii) {
                triangle t = polygons_[i].tri()[ii];
                CH ch;
                if (!(ve_[t[0]] && ve_[t[1]] && ve_[t[2]])) {
                    cmt.push_back(i);
                    continue;
                }
                int c0, c1, c2;
                if (!T_.is_facet(
                    id2vh_[t[0]], id2vh_[t[1]], id2vh_[t[2]],
                    ch, c0, c1, c2)) {
                    cmt.push_back(i);
                }
            }
        }
    }

    int processor::peel_by_winding_number(
        double thr
        ) {
        VERBOSE("Peeling... ");
        Eigen::MatrixXd BC(T_.number_of_finite_cells(), 3);

        Finite_cells_iterator cit = T_.finite_cells_begin();
        for (int ti = 0; cit != T_.finite_cells_end(); ++cit, ++ti) {
            Point c = CGAL::centroid(cit->vertex(0)->point().point(),
                cit->vertex(1)->point().point(),
                cit->vertex(2)->point().point(),
                cit->vertex(3)->point().point());
            BC(ti, 0) = c[0];
            BC(ti, 1) = c[1];
            BC(ti, 2) = c[2];
        }

        Eigen::VectorXd WN;
        Eigen::VectorXi CO;
        int nf = 0;
        for (int i = 0; i < np_; ++i)
            nf += polygons_[i].tri().size();

        Eigen::MatrixXi F(nf, 3);
        nf = 0;
        for (int i = 0; i < np_; ++i) {
            for (int j = 0; j < polygons_[i].tri().size(); j++)
                F.row(nf++) << polygons_[i].tri()[j][0], polygons_[i].tri()[j][1], polygons_[i].tri()[j][2];
        }

        igl::bfs_orient(F, F, CO);
        Eigen::MatrixXd EV(nv_, 3);
        for (int i = 0; i < nv_; ++i) {
            for (int j = 0; j < 3; ++j)
                EV(i, j) = vertices_[i][j];
        }
        igl::winding_number(EV, F, BC, WN);

        int npt = 0;
        cit = T_.finite_cells_begin();
        for (int ti = 0; cit != T_.finite_cells_end(); ++cit, ++ti)
            if (fabs(WN(ti)) > thr)
                cit->info() = 0;
            else {
                cit->info() = -1;
                npt++;
            }
        VERBOSE("done - " << npt << " tets removed");

        return npt;
    }

    bool processor::read_off(
        const std::string& file_path
        ) {
        std::ifstream in(file_path);

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
            std::istringstream iss1(line);
            double x, y, z, w = 0.0;
            if (woff) iss1 >> x >> y >> z >> w;
            else iss1 >> x >> y >> z;
            vertices_[i] = Point(x, y, z);
            weights_[i] = w;
        }

        /** read polygons **/
        polygons_.resize(np_);

        double mina = 180.0;
        for (int i = 0; i < np_; i++) {
            std::getline(in, line);
            std::istringstream iss2(line);
            int d;
            iss2 >> d;
            if (d < 3) WARNING("Degree < 3");

            auto id = static_cast<int*>(malloc(sizeof(int) * d));

            for (int v = 0; v < d; ++v) iss2 >> id[v];
            int prev_v = d - 1;
            for (int v = 0; v < d; ++v) {
                edge e(id[prev_v], id[v]);

                polygons_[i].insert_boundary_edge(e);
                E2P_[e].push_back(i);

                prev_v = v;
            }
            for (int v = 0; v < d; ++v) {
                if (double a = CGAL::approximate_angle(
                    vertices_[id[v]],
                    vertices_[id[(v + 1) % d]],
                    vertices_[id[(v + 2) % d]]);
                    a < mina
                    ) mina = a;
            }

            free(id);
        }

        return true;
    }

    bool processor::read_igl(
        const std::string& file_path
        ) {
        Eigen::MatrixXd EV;
        Eigen::MatrixXi F;
        igl::read_triangle_mesh(file_path, EV, F);
        nv_ = EV.rows();
        vertices_.resize(nv_);
        weights_.resize(nv_, 0.0);
        for (int i = 0; i < nv_; i++)
            vertices_[i] = Point(EV(i, 0), EV(i, 1), EV(i, 2));

        np_ = F.rows();
        polygons_ = std::vector<poly>(np_);
        for (int i = 0; i < np_; i++) {
            constexpr int d = 3;
            for (int v = 0; v < d; v++) {
                edge e{ F(i, v), F(i, (v + 1) % d) };
                polygons_[i].be().push_back(e);
                E2P_[e].push_back(i);
            }
        }

        return true;
    }

    bool processor::write_off(
        const std::string& file_path
        ) {
        std::ofstream out(file_path);

        if (!out.good()) {
            WARNING("Output file is not good!");
            return false;
        }

        out << "OFF" << std::endl;
        out << nv_ << " " << ntri() << " 0" << std::endl;
        for (int i = 0; i < nv_; i++)
            out << vertices_[i].x() << " " << vertices_[i].y() << " " << vertices_[i].z() << std::endl;
        for (int i = 0; i < np_; i++) {
            for (int ii = 0; ii < polygons_[i].tri().size(); ii++) {
                const int v0 = polygons_[i].tri()[ii][0];
                const int v1 = polygons_[i].tri()[ii][1];
                const int v2 = polygons_[i].tri()[ii][2];

                out << "3 " << v0 << " " << v1 << " " << v2 << std::endl;
            }
        }
        out.close();

        return true;
    }

    bool processor::write_woff(
        const std::string& file_path
        ) {
        std::ofstream out(file_path);

        if (!out.good()) {
            WARNING("Output file is not good!");
            return false;
        }

        out << "WOFF" << std::endl;
        out << nv_ << " " << np_ << " 0" << std::endl;
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (int i = 0; i < nv_; i++)
            out << vertices_[i].x() << " " << vertices_[i].y() << " " << vertices_[i].z() << " " << weights_[i] << std::endl;

        for (int i = 0; i < np_; i++) {
            out << polygons_[i].be().size();
            for (const auto& e : polygons_[i].be())
                out << " " << e[0];
            out << std::endl;
        }
        out.close();

        return true;
    }

    bool processor::write_tet_off(
        const std::string& file_path,
        const double& shrink_factor
        ) {
        std::ofstream out(file_path);

        if (!out.good()) {
            WARNING("Output file is not good!");
            return false;
        }

        out << "OFF" << std::endl;
        out << nv_ << " " << ntri() << " 0" << std::endl;
        for (int i = 0; i < nv_; i++)
            out << vertices_[i][0] << " " << vertices_[i][1] << " " << vertices_[i][2] << std::endl;

        out.close();

        return true;
    }
}
