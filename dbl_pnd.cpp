/*
 *
 * This software comes without any warranty to be correct or useful.
 * Licensed under GPLv3, creator: alcubierre-drive@github, May 2019.
 *
 */


#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include <cmath>

namespace ode = boost::numeric::odeint;

// 1. Ordnung benötigt, also l, \dot l, \phi, \dot \phi als eigene Variablen.
typedef boost::array<double,4> coords;

typedef struct {
    coords x;
    double t;
} step;

boost::array<double,2> to_caresian( coords& x ) {
    boost::array<double,2> result = { sin(x[0])*x[2], cos(x[0])*x[2] };
    return result;
}

double square( double x ) { return x*x; }

typedef struct {
    double g = 1.0;
    double D_M = 4.0;
    double l_0_const = 1.0;
} parameters;

std::vector<step> run_simulation( coords x0, double t_max, double dt,
        const parameters& param ) {
    double g = param.g;
    double D_M = param.D_M;
    double l_0_const = param.l_0_const;

    // lambda function for one step
    auto double_pendulum = [g,D_M,l_0_const] (const coords& x, coords& dx,
            double) {
        // phi
        dx[0] = x[1];
        // dphi
        dx[1] = -2.0 * x[3] * x[1] / x[2] - g / x[2] * sin(x[0]);
        // l
        dx[2] = x[3];
        // dl
        dx[3] = x[2] * square(x[1]) + g * cos(x[0]) - D_M * (x[2]-l_0_const);
    };

    // lambda function to save one step
    std::vector<step> steps;
    auto save = [&steps] (const coords& x, double t) {
        step s = {x,t};
        steps.push_back(s);
    };

    // integrate the ode
    ode::integrate( double_pendulum, x0, 0.0, t_max, dt, save);

    return steps;
}

int alternate( int i ) {
    // 0,1, 2,3, 4,5, 6,... --> 0,1,-1,2,-2,3,-3,...
    i += 1;
    return (1-2*(i%2))*i/2;
}

void flatten_angles( std::vector<step>& results ) {
    for (unsigned i=1; i<results.size(); ++i) {
        double delta = results[i].x[0] - results[i-1].x[0];
        int mod = 1;
#define maxmod 50
        for (mod=1; mod<maxmod; ++mod) {
            if (abs(delta+(double)alternate(mod)*2.*M_PI) < 2.*M_PI) {
                break;
            }
        }
        if (mod!=maxmod/2) {
            results[i].x[0] += (double)(alternate(mod))*2.*M_PI;
        }
    }
    double delta_01 = results[1].x[0] - results[0].x[0];
    int num_pi = round((delta_01 / 2.0 / M_PI));
    results[0].x[0] += num_pi*2.0*M_PI;
}

int main (int argc, char** argv) {

    if (argc>1 && !strcmp(argv[1],"-h")) {
        printf("usage: %s [nmax] [tmax] [dt]\n",argv[0]);
        exit(1);
    }
    int nmax = 8000; if (argc>1) { nmax = atoi(argv[1]); }
    double tmax = 200.0; if (argc>2) { tmax = atof(argv[2]); }
    double dt = 0.1; if (argc>3) { dt = atof(argv[3]); }
    double phimax = 2.0*M_PI;

    parameters p;
    coords x0 = { 0.0, 0.0, 1.25+0.2, 0.0 };

    std::vector<step> final_results(nmax);

// faster execution with omp using all available cores
#pragma omp parallel for schedule(dynamic)
    for (int i=0; i<nmax; ++i) {
        x0[0] = (double)i/((double)nmax-1)*phimax;
        std::vector<step> result = run_simulation( x0, tmax, dt, p );
        final_results[i] = result[result.size()-1];
    }

    // try to subtract 2pi for a smoother result
    flatten_angles(final_results);

    std::cout << "# phi0 phi_end" << std::endl;
    for (unsigned i=0; i<final_results.size(); ++i) {
        std::cout << (double)i/((double)nmax-1)*phimax << " " <<
            final_results[i].x[0] << std::endl;
    }

}
