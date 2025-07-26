#include "cwdt/cwdt.h"

int main(){
	using namespace CWDT;

    int nr = 100;
    int ns = 100;
    bool w = false;
    bool peel = false;
    double thr = 0.5;
    double height_factor = 1.0;
    double weight_factor = 1.0;
    double tolerance = 1e-9;

	processor P;
    P.read_igl("../data/spot_quadrangulated.obj");
	P.compute_planes();
	P.triangulate();
	int onv = P.nv();
	int onf = P.ntri();

	P.env()->start();

    int splits = 0;
    for (; splits < ns; splits++)
    {
        std::cout << "Loop " << splits << std::endl;
        int nif = P.solve(nr, tolerance, 1e2 * tolerance, height_factor, weight_factor);
        std::cout << "Number of ignored triangles: " << nif << std::endl;
        if (nif == 0)
            break;
        int nse = P.split_missing_edges();
        if (nse == 0)
            P.split_missing_faces();
        P.env()->resetParams();
    }
    std::vector<int> mf;
    P.missing_faces(mf);
    if (mf.size() == 0)
    {
        std::cout << "---\nFound regular triangulation of constraints after " << splits << " splits" << std::endl;
    }

    int npt = 0;
    if (peel)
        npt = P.peel_by_winding_number(thr);

    std::cerr << "Vertices: " << P.nv() << " ( " << onv << " + " << (P.nv() - onv) << " )" << std::endl;
    std::cout << "Surface Triangles: " << P.ntri() << " ( " << onf << " + " << (P.ntri() - onf) << " )" << std::endl;
    std::cout << "Tets: " << P.T().number_of_finite_cells() << " - " << npt << " (peeled) = " << (P.T().number_of_finite_cells() - npt) << std::endl;



    P.write_off("D://C_Project//CWDT//data//result.off");
    P.write_tet_off("D://C_Project//CWDT//data//result_tets.off");

	return 1;
}
