#ifndef OIL2D_HPP
#define OIL2D_HPP

#include <cmath>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>

#define MATRIX_CSV "matrix.csv"

std::string find_file(const std::string &aFileName, const std::string &aLookupStartDir = "/home")
{
    std::string aFilePath {""};
    for(const auto& entry : std::filesystem::recursive_directory_iterator(aLookupStartDir))
    {
        if(entry.is_regular_file() && entry.path().filename() == aFileName)
        {
            aFilePath = entry.path().string();
            break;
        }
    }
    return aFilePath;
}

std::vector<std::vector<double>> permeability_matrix;

std::vector<std::vector<double>> load_permeability_matrix(const std::string& file_path) {
    std::ifstream file(file_path);
    std::vector<std::vector<double>> permeability_matrix;
    std::string line;

    while (getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row;
        std::string val;
        
        while (getline(iss, val, ',')) {
        	double scaled_value = std::stod(val) * 100;
                row.push_back(scaled_value);  // Convert string to double
        }
        
        permeability_matrix.push_back(row);
    }

    return permeability_matrix;
}



namespace ads {

struct vec2d {
    double x;
    double y;
};

inline vec2d operator-(const vec2d& a, const vec2d& b) {
    return {a.x - b.x, a.y - b.y};
}

inline vec2d operator+(const vec2d& b, const vec2d& a) {
    return {a.x + b.x, a.y + b.y};
}

inline vec2d operator*(double c, const vec2d& a) {
    return {c * a.x, c * a.y};
}

inline double dot(const vec2d& a, const vec2d& b) {
    return a.x * b.x + a.y * b.y;
}

inline double len_sq(const vec2d& a) {
    return dot(a, a);
}

inline double len(const vec2d& a) {
    return std::sqrt(len_sq(a));
}

inline double falloff(double r, double R, double t) {
    if (t < r)
        return 1.0;
    if (t > R)
        return 0.0;
    double h = (t - r) / (R - r);
    return std::pow((h - 1) * (h + 1), 2);
}

// r < R in [0, 1]
inline double bump(double r, double R, double x, double y) {
    double dx = x - 0.5;
    double dy = y - 0.5;
    double t = std::sqrt(dx * dx + dy * dy);
    return falloff(r / 2, R / 2, t);
}

struct pumps {
    std::vector<ads::vec2d> sources;
    std::vector<ads::vec2d> sinks;

    static constexpr double radius = 0.05;
    static constexpr double pumping_strength = 2;
    static constexpr double draining_strength = 1;

    double pumping(double x, double y) const {
        ads::vec2d v{x, y};
        double p = 0;
        for (const auto& pos : sources) {
            double dist = len(v - pos);
            p += pumping_strength * ads::falloff(0, radius, dist);
        }
        return p;
    }

    double draining(double u, double x, double y) const {
        ads::vec2d v{x, y};
        double p = 0;
        for (const auto& pos : sinks) {
            double dist = len(v - pos);
            double s = draining_strength * ads::falloff(0, radius, dist);
            p += u * s;
        }
        return p;
    }

    auto pumping_fun() const {
        return [this](double x, double y) { return pumping(x, y); };
    }
};

class oil2d : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    galois_executor executor{4};

    pumps process = pumps{
    	{{0.2, 0.8}, {0.2, 0.65}, {0.25, 0.45}, {0.4, 0.45}, {0.45, 0.55}}, 
    	{{0.3, 0.65}, {0.3, 0.3}, {0.75, 0.8},}
    };
    lin::tensor<double, 4> kq;
    output_manager<2> output;

public:
    explicit oil2d(const config_2d& config)
    : Base{config}
    , u{shape()}
    , u_prev{shape()}
    , kq{{x.basis.elements, y.basis.elements, x.basis.quad_order + 1, y.basis.quad_order + 1}}
    , output{x.B, y.B, 100} { }

    double init_state(double x, double y) {
        double r = 0.1;
        double R = 0.5;
        return 1e-3 * ads::bump(r, R, x, y);
    };

private:
    void before() override {
        // Load the binary matrix representing of geological layers
       	std::string permeability_matrix_path = find_file(MATRIX_CSV);
        permeability_matrix = load_permeability_matrix(permeability_matrix_path);
       
        
        fill_permeability_map();
        prepare_matrices();
        
        auto init = [this](double x, double y) { return init_state(x, y); };
        projection(u, init);
        solve(u);
        output.to_file(u, "out_%d.data", 0);
    }

    
    void fill_permeability_map() {
        for (auto e : elements()) {
            for (auto q : quad_points()) {
            	auto[x, y] = point(e, q); 

            	// X, y to index into the matrix
            	int matrix_x = static_cast<int>(x * 19); 
            	int matrix_y = static_cast<int>(y * 19);  

            	// Boundary check for safety
            	matrix_x = std::max(0, std::min(matrix_x, 19));
            	matrix_y = std::max(0, std::min(matrix_y, 19));

            	kq(e[0], e[1], q[0], q[1]) = permeability_matrix[matrix_y][matrix_x];
            	
            	// Drukowanie dla sprawdzenia czy dobrze wczytuje wartosci
            	//std::cout << "Permeability at (" << x << ", " << y << ") -> Matrix(" << matrix_x << ", " << matrix_y 
                      //<< ") = " << permeability_matrix[matrix_y][matrix_x] << std::endl;
            	
            }
        }
    }
    

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
    }

    void step(int /*iter*/, double t) override {
        compute_rhs(t);
        solve(u);
    }

    void compute_rhs(double t) {
        auto& rhs = u;

        zero(rhs);
        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);

                double mi = 10;
                double k = permeability(e, q);
                value_type u = eval_fun(u_prev, e, q);
                double h = forcing(x, t, u.val);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    double val = -k * std::exp(mi * u.val) * grad_dot(u, v) + h * v.val;
                    U(aa[0], aa[1]) += (u.val * v.val + steps.dt * val) * w * J;
                }
            }
            executor.synchronized([&] { update_global_rhs(rhs, U, e); });
        });
    }

    double energy(const vector_type& u) const {
        double E = 0;
        for (auto e : elements()) {
            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                value_type a = eval_fun(u, e, q);
                E += a.val * a.val * w * J;
            }
        }
        return E;
    }

    void after_step(int iter, double /*t*/) override {
        if (iter % 10 == 0) {
            std::cout << "Step " << iter << ", energy: " << energy(u) << std::endl;
        }
        if ((iter + 1) % 100 == 0) {
            output.to_file(u, "out_%d.data", iter + 1);
        }
    }

    double permeability(index_type e, index_type q) const { 
    	auto x = point(e, q);
    	return permeability_matrix[int(x[0]*20)][int(x[1]*20)]; 
    }

    double forcing(point_type x, double /*t*/, double u) const {
        return process.pumping(x[0], x[1]) - process.draining(x[0], x[1], u);
    }
};

}  // namespace ads

#endif  // OIL2D_HPP
