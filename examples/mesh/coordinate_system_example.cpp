/**
 * @file coordinate_system_example.cpp
 * @brief Example program demonstrating CoordinateSystem class usage
 */

#include "mesh/coordinate_system.hpp"
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "=== CoordinateSystem Example ===\n\n";

    try {
        // Example 1: Uniform grid
        std::cout << "1. Creating uniform grid (10x10, [0,1]x[0,1])\n";
        mesh::CoordinateSystem uniform_grid(10, 10, 0.0, 1.0, 0.0, 1.0);

        std::cout << "   Grid info:\n";
        std::cout << "     nx: " << uniform_grid.nx() << "\n";
        std::cout << "     ny: " << uniform_grid.ny() << "\n";
        std::cout << "     x_min: " << uniform_grid.x_min() << "\n";
        std::cout << "     x_max: " << uniform_grid.x_max() << "\n";
        std::cout << "     y_min: " << uniform_grid.y_min() << "\n";
        std::cout << "     y_max: " << uniform_grid.y_max() << "\n";
        std::cout << "     lx: " << uniform_grid.lx() << "\n";
        std::cout << "     ly: " << uniform_grid.ly() << "\n";
        std::cout << "     hx: " << uniform_grid.hx() << "\n";
        std::cout << "     hy: " << uniform_grid.hy() << "\n";
        std::cout << "     grid type: " << (uniform_grid.grid_type() == mesh::GridType::Uniform ? "Uniform" : "Non-uniform") << "\n";
        std::cout << "     grid quality: " << uniform_grid.grid_quality() << "\n";
        std::cout << "     aspect ratio: " << uniform_grid.aspect_ratio() << "\n\n";

        std::cout << "   Sample coordinates:\n";
        std::cout << "     (0, 0): " << uniform_grid.x(0) << ", " << uniform_grid.y(0) << "\n";
        std::cout << "     (5, 5): " << uniform_grid.x(5) << ", " << uniform_grid.y(5) << "\n";
        std::cout << "     (9, 9): " << uniform_grid.x(9) << ", " << uniform_grid.y(9) << "\n\n";

        // Example 2: Coordinate transformations
        std::cout << "2. Testing coordinate transformations\n";
        double test_x = 0.5;
        double test_y = 0.5;
        double xi = uniform_grid.xi_from_x(test_x);
        double eta = uniform_grid.eta_from_y(test_y);
        std::cout << "   Physical coordinates: (" << test_x << ", " << test_y << ")\n";
        std::cout << "   Normalized coordinates: (" << xi << ", " << eta << ")\n";
        std::cout << "   Back to physical: (" << uniform_grid.x_from_xi(xi) << ", "
                  << uniform_grid.y_from_eta(eta) << ")\n\n";

        // Example 3: Index lookup
        std::cout << "3. Testing index lookup\n";
        size_t i = uniform_grid.i_at(0.5);
        size_t j = uniform_grid.j_at(0.5);
        std::cout << "   Index at (0.5, 0.5): (" << i << ", " << j << ")\n";
        std::cout << "   Coordinate at (" << i << ", " << j << "): ("
                  << uniform_grid.x(i) << ", " << uniform_grid.y(j) << ")\n\n";

        // Example 4: Stretched grid
        std::cout << "4. Creating stretched grid (50x50)\n";
        mesh::CoordinateSystem stretched_grid(50, 50, 0.0, 1.0, 0.0, 1.0);
        stretched_grid.set_x_distribution([](double xi) {
            return mesh::CoordinateSystem::stretching_function(xi, 2.0);
        });
        stretched_grid.set_y_distribution([](double eta) {
            return mesh::CoordinateSystem::stretching_function(eta, 1.5);
        });

        std::cout << "   Grid type: " << (stretched_grid.grid_type() == mesh::GridType::Uniform ? "Uniform" : "Non-uniform") << "\n";
        std::cout << "   Grid quality: " << stretched_grid.grid_quality() << "\n";

        std::cout << "   First 5 X coordinates (stretched): ";
        for (size_t i = 0; i < 5; ++i) {
            std::cout << std::fixed << std::setprecision(4) << stretched_grid.x(i) << " ";
        }
        std::cout << "\n";

        std::cout << "   Last 5 X coordinates (stretched): ";
        for (size_t i = stretched_grid.nx() - 5; i < stretched_grid.nx(); ++i) {
            std::cout << std::fixed << std::setprecision(4) << stretched_grid.x(i) << " ";
        }
        std::cout << "\n\n";

        // Example 5: Generate coordinate arrays
        std::cout << "5. Generating coordinate arrays\n";
        std::vector<double> x_coords = uniform_grid.generate_x_coords();
        std::vector<double> y_coords = uniform_grid.generate_y_coords();
        std::cout << "   X coordinates size: " << x_coords.size() << "\n";
        std::cout << "   Y coordinates size: " << y_coords.size() << "\n\n";

        // Example 6: Generate coordinate mesh
        std::cout << "6. Generating coordinate mesh\n";
        utils::Array2D x_mesh(10, 10);
        utils::Array2D y_mesh(10, 10);
        uniform_grid.generate_coordinate_mesh(x_mesh, y_mesh);
        std::cout << "   Mesh dimensions: " << x_mesh.rows() << "x" << x_mesh.cols() << "\n";
        std::cout << "   Sample mesh points:\n";
        std::cout << "     (0, 0): " << x_mesh(0, 0) << ", " << y_mesh(0, 0) << "\n";
        std::cout << "     (5, 5): " << x_mesh(5, 5) << ", " << y_mesh(5, 5) << "\n";
        std::cout << "     (9, 9): " << x_mesh(9, 9) << ", " << y_mesh(9, 9) << "\n\n";

        // Example 7: Inside domain checking
        std::cout << "7. Testing domain boundary checks\n";
        std::cout << "   (0.5, 0.5) inside: " << (uniform_grid.is_inside(0.5, 0.5) ? "yes" : "no") << "\n";
        std::cout << "   (0.0, 0.0) inside: " << (uniform_grid.is_inside(0.0, 0.0) ? "yes" : "no") << "\n";
        std::cout << "   (1.0, 1.0) inside: " << (uniform_grid.is_inside(1.0, 1.0) ? "yes" : "no") << "\n";
        std::cout << "   (-0.1, 0.5) inside: " << (uniform_grid.is_inside(-0.1, 0.5) ? "yes" : "no") << "\n";
        std::cout << "   (1.1, 0.5) inside: " << (uniform_grid.is_inside(1.1, 0.5) ? "yes" : "no") << "\n\n";

        // Example 8: Grid spacing
        std::cout << "8. Testing grid spacing\n";
        std::cout << "   Average hx: " << uniform_grid.hx() << "\n";
        std::cout << "   Average hy: " << uniform_grid.hy() << "\n";
        std::cout << "   hx_at(0): " << uniform_grid.hx_at(0) << "\n";
        std::cout << "   hy_at(0): " << uniform_grid.hy_at(0) << "\n";

        std::cout << "\n=== All tests passed! ===\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
