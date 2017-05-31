#ifndef POINTMASSLINEAR_H_
#define POINTMASSLINEAR_H_

/*
 * N-dimensional 2nd order point robot tracking linear desired dynamics with quadratic cost.
 * State, x = [x_rob ; (x_d - x_rob)]
 * Cost_to_go = x^T Q x + u^T R u
 * Final_cost = x^T Q_track_f x
*/

#include "types.h"
#include <eigen3/Eigen/Dense>
#include <OCProblemFH.h>

namespace ileqg {

class PointMassLinear: public ileqg::OCProblemFH
{
private :
    Matrix A_, B_, A_d_, Q_, Q_f_, R_;                          // State and cost matrices of the full compound state  x = [x_rob ; (x_d - x_rob)]
    Vector x_d_center_;


public :

    //Robot linear dynamics A_rob, B_rob and desired trajectory dynamics A_d
    // Q_track and Q_track_f are the weighting for the tracking error (x_d - x_rob)^T Q_track (x_d - x_rob)
    PointMassLinear(int time_horizon, const Matrix & A_rob, const Matrix & B_rob, const Matrix & A_d, const Vector & x_d_center, const Matrix & Q_track, const Matrix & Q_track_f, const Matrix & R)
        : OCProblemFH(time_horizon), R_(R) {
        reset(  A_rob,  B_rob, A_d, x_d_center,  Q_track,  Q_track_f);
    }

    void reset(const Matrix & A_rob, const Matrix & B_rob, const Matrix & A_d, const Vector & x_d_center, const Matrix & Q_track, const Matrix & Q_track_f){

        size_t x_dim = A_rob.rows()*2;
        size_t u_dim = B_rob.cols();
        A_.resize(x_dim, x_dim);
        A_.setZero();
        B_.resize(x_dim, u_dim);
        B_.setZero();
        Q_.resize(x_dim, x_dim);
        Q_.setZero();
        Q_f_.resize(x_dim, x_dim);
        Q_f_.setZero();
        A_d_.resize(x_dim/2, x_dim/2);
        A_d_.setZero();
        x_d_center_.resize(x_dim/2);
        x_d_center_.setZero();

        A_.topLeftCorner(x_dim/2, x_dim/2) = A_rob;
        A_.bottomLeftCorner(x_dim/2, x_dim/2) = - A_d - A_rob;

        Q_.bottomRightCorner(x_dim/2, x_dim/2) = Q_track;

        Q_f_.bottomRightCorner(x_dim/2, x_dim/2) = Q_track_f;

        A_d_ = A_d;
        x_d_center_ = x_d_center;

        B_ << B_rob , -B_rob;
    }

    // Dynamics
    void dynamics(const Vector & x, const Vector & u, Vector & x_dot){
        x_dot = A_*x + B_*u;
        x_dot.tail(x.rows()/2) = x_dot.tail(x.rows()/2) + (A_d_ * x_d_center_);
    }

    // Dynamics with 1st order state and control derivatives
    void dynamics(const Vector & x, const Vector & u, Matrix & A, Matrix & B) {
        A = A_;
        B = B_;
    }

    // Cost at sample k
    double cost(const Vector & x, const Vector & u, const int & k) {
        return ((x.transpose() * (Q_ * x)) + (u.transpose() * (R_ * u)))(0,0);
    }

    // Cost at sample T
    double finalCost(const Vector & x, const Vector & u) {
        return x.transpose() * Q_f_ * x;
    }

    // Cost at sample k with 2nd order state and control derivatives
    double cost(const Vector & x, const Vector & u, const int & k, double & q_0, Vector & q, Vector & r, Matrix & Q, Matrix & R) {
        R = 2*R_;
        Q = 2*Q_;
        r = 2 * R_ * u;
        q = 2 * Q_ * x;
        q_0 = cost(x,u,k);
        return q_0;
    }

    // Cost at sample T with 2nd order state derivatives
    double finalCost(const Vector & x, const Vector & u, double & q_0, Vector & q, Matrix & Q) {
        Q = 2*Q_f_;
        q = 2 * Q_f_ * x;
        q_0 = finalCost(x,u);
        return q_0;
    }
};

}

#endif /* POINTMASSLINEAR_H_ */
