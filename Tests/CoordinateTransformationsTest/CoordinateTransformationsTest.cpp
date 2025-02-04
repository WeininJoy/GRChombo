/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "IntVect.H"

// Other includes
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"
#include <limits>

// Chombo namespace
#include "UsingNamespace.H"

constexpr int ulp = 15; /* units in the last place */

bool almost_equal(double value, double correct_value,
                  int a_ulp /* units in the last place */)
{
    double diff = fabs(value - correct_value);
    double epsilon = std::numeric_limits<double>::epsilon();
    return (diff <
            std::max(epsilon * std::abs(correct_value + value), epsilon) *
                a_ulp);
}

bool check_tensor(const Tensor<2, double> &tensor,
                  const Tensor<2, double> &correct_tensor,
                  const std::string &test_name)
{
    bool failed = false;
    FOR(i, j)
    {
        if (!almost_equal(tensor[i][j], correct_tensor[i][j], ulp))
        {
            std::cout << "Failed " << test_name << " in component [" << i
                      << "][" << j << "]\n";
            std::cout << "value: " << tensor[i][j] << "\n";
            std::cout << "correct_value: " << correct_tensor[i][j] << std::endl;
            failed = true;
        }
    }
    return failed;
}

bool check_vector(const Tensor<1, double> &vector,
                  const Tensor<1, double> &correct_vector,
                  const std::string &test_name)
{
    bool failed = false;
    FOR(i)
    {
        if (!almost_equal(vector[i], correct_vector[i], ulp))
        {
            std::cout << "Failed " << test_name << " in component [" << i
                      << "]\n";
            std::cout << "value: " << vector[i] << "\n";
            std::cout << "correct_value: " << correct_vector[i] << std::endl;
            failed = true;
        }
    }
    return failed;
}

int main()
{
    bool failed = false;
    const double dx = 0.1;
    IntVect iv;
    iv[0] = 0;
    iv[1] = 0;
    iv[2] = 0;

    Coordinates<double> coords(iv, dx);
    const double x = coords.x;
    const double y = coords.y;
    const double z = coords.z;
    const double r = coords.get_radius();
    double rho2 = simd_max(x * x + y * y, 1e-12);
    double rho = sqrt(rho2);
    double r2sin2theta = rho2;

    /* for debugging
    std::cout << "x " << x << std::endl;
    std::cout << "y " << y << std::endl;
    std::cout << "z " << z << std::endl;
    std::cout << "r " << r << std::endl;
    std::cout << "rho " << rho << std::endl;
    */

    using namespace TensorAlgebra;
    using namespace CoordinateTransformations;

    // 1. Test spherical <-> cartesian coordinate
    
    // Test if inv_jac is really the inverse of the jacobian
    Tensor<2, double> jac_spher = spherical_jacobian(x, y, z);
    Tensor<2, double> inv_jac_spher = inverse_spherical_jacobian(x, y, z);
    Tensor<2, double> inv_jac_check_spher = compute_inverse(jac_spher);
    failed |= check_tensor(inv_jac_spher, inv_jac_check_spher, "inverse_jacobian (spherical)");

    // Test tensor transformations
    Tensor<2, double> Mij_cart;
    FOR(i, j) { Mij_cart[i][j] = 0.; }
    Mij_cart[0][0] = 1.;
    Mij_cart[1][1] = 1.;
    Mij_cart[2][2] = 1.;

    Tensor<2, double> Mij_spher;
    FOR(i, j) { Mij_spher[i][j] = 0.; }
    Mij_spher[0][0] = 1.;
    Mij_spher[1][1] = r * r;
    Mij_spher[2][2] = r2sin2theta;

    // Test cartesian_to_spherical_LL
    Tensor<2, double> Mij_spher_check;
    Mij_spher_check = cartesian_to_spherical_LL(Mij_cart, x, y, z);
    failed |=
        check_tensor(Mij_spher_check, Mij_spher, "cartesian_to_spherical_LL");

    // Test spherical_to_cartesian_LL
    Tensor<2, double> Mij_cart_check;
    Mij_cart_check = spherical_to_cartesian_LL(Mij_spher, x, y, z);
    failed |=
        check_tensor(Mij_cart_check, Mij_cart, "spherical_to_cartesian_LL");

    // Test cartesian_to_spherical_UU
    Tensor<2, double> Mij_spher_UU;
    Tensor<2, double> Mij_spher_UU_check;
    Mij_spher_UU_check =
        cartesian_to_spherical_UU(compute_inverse_sym(Mij_cart), x, y, z);
    Mij_spher_UU = compute_inverse_sym(Mij_spher);
    failed |= check_tensor(Mij_spher_UU_check, Mij_spher_UU,
                           "cartesian_to_spherical_UU");

    // Test spherical_to_cartesian_UU
    Tensor<2, double> Mij_cart_UU;
    Tensor<2, double> Mij_cart_UU_check;
    Mij_cart_UU_check =
        spherical_to_cartesian_UU(compute_inverse_sym(Mij_spher), x, y, z);
    Mij_cart_UU = compute_inverse_sym(Mij_cart);
    failed |= check_tensor(Mij_cart_UU_check, Mij_cart_UU,
                           "spherical_to_cartesian_UU");

    // Test vector transformations
    Tensor<1, double> si_cart;
    si_cart[0] = x / r;
    si_cart[1] = y / r;
    si_cart[2] = z / r;

    Tensor<1, double> si_spher;
    si_spher[0] = 1.0;
    si_spher[1] = 0.0;
    si_spher[2] = 0.0;

    // Test cartesian_to_spherical_U
    Tensor<1, double> si_spher_U_check;
    si_spher_U_check = cartesian_to_spherical_U(si_cart, x, y, z);
    failed |=
        check_vector(si_spher_U_check, si_spher, "cartesian_to_spherical_U");

    // Test spherical_to_cartesian_U
    Tensor<1, double> si_cart_U_check;
    si_cart_U_check = spherical_to_cartesian_U(si_spher, x, y, z);
    failed |=
        check_vector(si_cart_U_check, si_cart, "spherical_to_cartesian_U");

    // Test cartesian_to_spherical_L
    Tensor<1, double> si_spher_L_check;
    si_spher_L_check = cartesian_to_spherical_L(si_cart, x, y, z);
    failed |=
        check_vector(si_spher_L_check, si_spher, "cartesian_to_spherical_L");

    // Test spherical_to_cartesian_L
    Tensor<1, double> si_cart_L_check;
    si_cart_L_check = spherical_to_cartesian_L(si_spher, x, y, z);
    failed |=
        check_vector(si_cart_L_check, si_cart, "spherical_to_cartesian_L");

    // Test area_element_sphere
    double area_element_spher = r * sqrt(rho2);
    double area_element_check_spher = area_element_sphere(Mij_spher);
    if (!almost_equal(area_element_spher, area_element_check_spher, ulp))
    {
        std::cout << "Failed area_element (sphere)\n";
        std::cout << "value: " << area_element_spher << "\n";
        std::cout << "correct_value: " << area_element_check_spher << std::endl;
        failed = true;
    }


    // 2. Test cylindrical <-> cartesian coordinate

    // Test if inv_jac is really the inverse of the jacobian
    Tensor<2, double> jac_cylin = cylindrical_jacobian(x, y, z);
    Tensor<2, double> inv_jac_cylin = inverse_cylindrical_jacobian(x, y, z);
    Tensor<2, double> inv_jac_check_cylin = compute_inverse(jac_cylin);
    failed |= check_tensor(inv_jac_cylin, inv_jac_check_cylin, "inverse_jacobian (cylindrical)");

    // Test tensor transformations

    Tensor<2, double> Mij_cylin;
    FOR(i, j) { Mij_cylin[i][j] = 0.; }
    Mij_cylin[0][0] = 1.;
    Mij_cylin[1][1] = rho2;
    Mij_cylin[2][2] = 1.;

    // Test cartesian_to_cylindrical_LL
    Tensor<2, double> Mij_cylin_check;
    Mij_cylin_check = cartesian_to_cylindrical_LL(Mij_cart, x, y, z);
    failed |=
        check_tensor(Mij_cylin_check, Mij_cylin, "cartesian_to_cylindrical_LL");

    // Test cylindrical_to_cartesian_LL
    Mij_cart_check = cylindrical_to_cartesian_LL(Mij_cylin, x, y, z);
    failed |=
        check_tensor(Mij_cart_check, Mij_cart, "cylindrical_to_cartesian_LL");

    // Test cartesian_to_cylindrical_UU
    Tensor<2, double> Mij_cylin_UU;
    Tensor<2, double> Mij_cylin_UU_check;
    Mij_cylin_UU_check =
        cartesian_to_cylindrical_UU(compute_inverse_sym(Mij_cart), x, y, z);
    Mij_cylin_UU = compute_inverse_sym(Mij_cylin);
    failed |= check_tensor(Mij_cylin_UU_check, Mij_cylin_UU,
                           "cartesian_to_cylindrical_UU");

    // Test cylindrical_to_cartesian_UU
    Mij_cart_UU_check =
        cylindrical_to_cartesian_UU(compute_inverse_sym(Mij_cylin), x, y, z);
    Mij_cart_UU = compute_inverse_sym(Mij_cart);
    failed |= check_tensor(Mij_cart_UU_check, Mij_cart_UU,
                           "cylindrical_to_cartesian_UU");

    // Test vector transformations
    si_cart[0] = x / rho;
    si_cart[1] = y / rho;
    si_cart[2] = 0.0;

    Tensor<1, double> si_cylin;
    si_cylin[0] = 1.0;
    si_cylin[1] = 0.0;
    si_cylin[2] = 0.0;

    // Test cartesian_to_cylindrical_U
    Tensor<1, double> si_cylin_U_check;
    si_cylin_U_check = cartesian_to_cylindrical_U(si_cart, x, y, z);
    failed |=
        check_vector(si_cylin_U_check, si_cylin, "cartesian_to_cylindrical_U");

    // Test cylindrical_to_cartesian_U
    si_cart_U_check = cylindrical_to_cartesian_U(si_cylin, x, y, z);
    failed |=
        check_vector(si_cart_U_check, si_cart, "cylindrical_to_cartesian_U");

    // Test cartesian_to_cylindrical_L
    Tensor<1, double> si_cylin_L_check;
    si_cylin_L_check = cartesian_to_cylindrical_L(si_cart, x, y, z);
    failed |=
        check_vector(si_cylin_L_check, si_cylin, "cartesian_to_cylindrical_L");

    // Test cylindrical_to_cartesian_L
    si_cart_L_check = cylindrical_to_cartesian_L(si_cylin, x, y, z);
    failed |=
        check_vector(si_cart_L_check, si_cart, "cylindrical_to_cartesian_L");

    // Test area_element_cylinder
    
    double area_element_cylin = sqrt(rho2);
    double area_element_check_cylin = area_element_cylinder(Mij_cylin);
    if (!almost_equal(area_element_cylin, area_element_check_cylin, ulp))
    {
        std::cout << "Failed area_element (cylinder)\n";
        std::cout << "value: " << area_element_cylin << "\n";
        std::cout << "correct_value: " << area_element_check_cylin << std::endl;
        failed = true;
    }
    

    if (failed)
    {
        std::cout << "Coordinate transformations test failed..." << std::endl;
        return failed;
    }
    else
    {
        std::cout << "Coordinate transformations test passed..." << std::endl;
        return 0;
    }
}
