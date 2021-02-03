// Abstract class for a regge trajectory. Any user-defined class for a specific
// for of the RT can be used in amplitudes (e.g. pomeron_exchange) as long as it
// is derived from this class
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _REGGE_TRAJ_
#define _REGGE_TRAJ_

#include <complex>
#include <string>

class regge_trajectory
{
    public:

    // constructor
    regge_trajectory(std::string name = "")
    : _parent(name)
    {};

    regge_trajectory(int sig, std::string name = "")
    : _signature(sig), _parent(name)
    {};

    regge_trajectory(int sig, int min_j, std::string name = "")
    : _signature(sig), _minJ(min_j), _parent(name)
    {};

    // copy constructor
    regge_trajectory(const regge_trajectory & old)
    : _parent(old._parent), _signature(old._signature)
    {};

    // Only need a function to evaluate the trajectory at some s
    virtual std::complex<double> eval(double s) = 0;

    virtual std::complex<double> slope(double s = 0.){return 0.;};

    inline void set_minimum_spin(int j)
    {
        _minJ = j;
    }

    // These parameters define the trajectory
    // name, spin, and mass of the lowest lying resonance on the parent trajectory
    std::string _parent;
    int _signature;
    int _minJ = 0;
};


// Basic linear regge_trajectory
class linear_trajectory : public regge_trajectory
{
    private:

    // Intercept and slope
    double _a0, _aprime;

    public:

    // Empty constructor
    linear_trajectory(){};

    // Parameterized constructor
    linear_trajectory(int sig, double inter, double slope, std::string name = "")
    : regge_trajectory(sig, name),
      _a0(inter), _aprime(slope)
    {};

    // copy Constructor
    linear_trajectory(const linear_trajectory & old)
    : regge_trajectory(old),
      _a0(old._a0), _aprime(old._aprime)
    {};

    // Setting utility
    void set_params(double inter, double slope)
    {
        _a0 = inter; _aprime = slope;
    };

    std::complex<double> eval(double s)
    {
        return _a0 + _aprime * s;
    };

    std::complex<double> slope(double s = 0.)
    {
        return _aprime;
    };
};

#endif
