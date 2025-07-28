#include "cwdt/cwdt.h"
#include "CLI11/CLI11.hpp"

using namespace CWDT;

int main(int argc, char** argv) {
    /* Settings */
    std::string input_mesh_path;
    bool not_peel = false;
    int ns = 100;
    double thr = 0.5;
    int nr = 100;
    double height_factor = 1.0;
    double weight_factor = 1.0;
    double tolerance = 1e-9;
    double shrink_factor = 0.3;

    /* App */
    CLI::App app{"CWDT"};
    argv = app.ensure_utf8(argv);

    app.add_option("--ns", ns,
        "Maximum number of splits (default: 100).");
    app.add_option("--nr", nr,
        "Maximum number of rounds (default: 100).");
    app.add_option("--hf", height_factor,
        "Height factor (default: 1.0).");
    app.add_option("--wf", weight_factor,
        "Weight factor (default: 1.0).");
    app.add_option("--tol", tolerance,
        "Minimum tolerance (default: 1e-9).");
    app.add_flag("--np", not_peel,
        "Not utilize the winding number to eliminate external tetrahedra (default: false).");
    app.add_option("--pthr", thr,
        "Threshold during the peel operation (default: 0.5), where tetrahedra with a winding number of "
        "their centroids less than this thr are identified as external.");
    app.add_option("--shrink", shrink_factor,
        "The shrink factor of each tetrahedron in the result used for visualization (default: 0.3).");
    app.add_option(
        "input_mesh_path",
        input_mesh_path,
        "Constrained mesh path."
        )->check(CLI::ExistingFile)->required();

    CLI11_PARSE(app, argc, argv);

    // std::cout << "ns: " << ns << std::endl;
    // std::cout << "nr: " << nr << std::endl;
    // std::cout << "hf: " << height_factor << std::endl;
    // std::cout << "wf: " << weight_factor << std::endl;
    // std::cout << "tol: " << tolerance << std::endl;
    // std::cout << "np: " << not_peel << std::endl;
    // std::cout << "pthr: " << thr << std::endl;
    // std::cout << "s: " << shrink_factor << std::endl;
    // std::cout << "input_mesh_path: " << input_mesh_path << std::endl;
    // return 1;

    /* Process */
    const std::filesystem::path absolute_path = std::filesystem::absolute(input_mesh_path);
    const std::filesystem::path parent_dir = absolute_path.parent_path();
    const std::string output_path = parent_dir.string() + "/";

    processor P;
    P.read_igl(input_mesh_path);
	P.compute_planes();
	P.triangulate();
	const int onv = P.nv();
	const int onf = P.ntri();

	P.env()->start();

    int splits = 0;
    for (; splits < ns; splits++) {
        std::cout << "Loop " << splits << std::endl;

        const int nif = P.solve(nr, tolerance, 1e2 * tolerance, height_factor, weight_factor);
        std::cout << "Number of ignored triangles: " << nif << std::endl;
        if (nif == 0)
            break;

        if (const int nse = P.split_missing_edges();
            nse == 0)
            P.split_missing_faces();
        P.env()->resetParams();
    }

    std::vector<int> mf;
    P.missing_faces(mf);
    if (mf.size() == 0)
        std::cout << "---\nFound regular triangulation of constraints after " << splits << " splits" << std::endl;

    int npt = 0;
    if (!not_peel)
        npt = P.peel_by_winding_number(thr);

    std::cout << "Vertices: " << P.nv() << " ( " << onv << " + " << (P.nv() - onv) << " )" << std::endl;
    std::cout << "Surface Triangles: " << P.ntri() << " ( " << onf << " + " << (P.ntri() - onf) << " )" << std::endl;
    std::cout << "Tets: " << P.T().number_of_finite_cells() << " - " << npt << " (peeled) = " << (P.T().number_of_finite_cells() - npt) << std::endl;

    P.write_off(output_path + "result_surface_mesh.off");
    P.write_woff(output_path + "result_weighted_surface_mesh.off");
    P.write_tets_obj(output_path + "result_tets.obj", shrink_factor, !not_peel);
    P.write_tets_mesh(output_path + "result_tets.mesh");

	return 1;
}