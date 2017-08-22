#ifndef POINTMASSLINEAROBSTACLE_H_
#define POINTMASSLINEAROBSTACLE_H_

/*
 * N-dimensional 2nd order point robot tracking linear desired dynamics with quadratic cost.
 * State, x = [x_rob ; (x_d - x_rob) ; (x_obstacle - r_rob)]
 * Cost_to_go = x^T Q x + u^T R u
 * Final_cost = x^T Q_track_f x
*/

#include "types.h"
#include <eigen3/Eigen/Dense>
#include <OCProblemFH.h>
#include <vector>

namespace ileqg {

class PointMassLinearObstacle: public ileqg::SOCProblemFH
{
private :
    Matrix A_, B_, A_d_, A_obst_, Q_track_, Q_track_f_, R_, W_;                          // State and cost matrices of the full compound state  x = [x_rob ;  (x_d - x_rob) ; (x_obstacle - r_rob)]
    Vector x_d_center_, x_obstacle_center_;
    std::vector < std::vector <Matrix> > Sigma_;
    double f_;


public :

    //Robot linear dynamics A_rob, B_rob and desired trajectory dynamics A_d
    // Q_track and Q_track_f are the weighting for the tracking error (x_d - x_rob)^T Q_track (x_d - x_rob)
    PointMassLinearObstacle(int time_horizon, const Matrix & A_rob, const Matrix & B_rob, const Matrix & A_d, const Matrix & A_obst,const Vector & x_d_center,
                            const Vector & x_obstacle_center, const Matrix & Q_track, const Matrix & Q_track_f, const Matrix & R, const std::vector <std::vector <Matrix> > & Sigma,
                            const Matrix & W, double f)
        : SOCProblemFH(time_horizon), R_(R)
    {
        reset(A_rob, B_rob, A_d, A_obst , x_d_center, x_obstacle_center, Q_track, Q_track_f, Sigma, W,f);
    }

    void reset(const Matrix & A_rob, const Matrix & B_rob, const Matrix & A_d, const Matrix & A_obst, const Vector & x_d_center, const Vector & x_obstacle_center,
               const Matrix & Q_track, const Matrix & Q_track_f, const std::vector <std::vector <Matrix> > & Sigma, const Matrix & W, double f){

        size_t x_dim = A_rob.rows()*3;
        size_t u_dim = B_rob.cols();
        A_.resize(x_dim, x_dim);
        A_.setZero();
        B_.resize(x_dim, u_dim);
        B_.setZero();
        Q_track_.resize(x_dim/3, x_dim/3);
        Q_track_.setZero();
        Q_track_f_.resize(x_dim/3, x_dim/3);
        Q_track_f_.setZero();
        A_d_.resize(x_dim/3, x_dim/3);
        A_d_.setZero();
        W_.resize(x_dim/3, x_dim/3);
        W_.setZero();
        A_obst_.resize(x_dim/3, x_dim/3);
        A_obst_.setZero();
        x_d_center_.resize(x_dim/3);
        x_d_center_.setZero();
        x_obstacle_center_.resize(x_dim/3);
        x_obstacle_center_.setZero();

        A_.topLeftCorner(x_dim/3, x_dim/3) = A_rob;

        A_.block (x_dim/3, 0 , x_dim/3, x_dim/3) = -A_d - A_rob;
        A_.block (x_dim/3, x_dim/3 , x_dim/3, x_dim/3) = -A_d;

        A_.bottomLeftCorner(x_dim/3, x_dim/3) = -A_obst - A_rob;
        A_.bottomRightCorner(x_dim/3, x_dim/3) = -A_obst;


        Q_track_ = Q_track;
        Q_track_f_ = Q_track_f;
        A_d_ = A_d;
        A_obst_ = A_obst;
        W_ = W;
        f_ = f;

        x_d_center_ = x_d_center;
        x_obstacle_center_ = x_obstacle_center;

        B_ << B_rob , -B_rob, -B_rob;

        Sigma_.resize(Sigma.size());
        for (size_t i=0; i<Sigma.size(); i++) {
            Sigma_[i].resize(Sigma[i].size());
             // x_dim,x_dim
             for (size_t j = 0 ; j < Sigma[i].size(); j++) {
                 Sigma_[i][j] = Matrix(x_dim,x_dim);
                 Sigma_[i][j] = Sigma[i][j];
              }
        }
    }

    // Dynamics
    void dynamics(const Vector & x, const Vector & u, Vector & x_dot){
        x_dot = A_*x + B_*u;
        x_dot.tail(x.rows()/3) = x_dot.tail(x.rows()/3) + (A_obst_ * x_obstacle_center_);
        x_dot.block(x.rows()/3,0,x.rows()/3,1) = x_dot.block(x.rows()/3, 0, x.rows()/3, 1) + (A_d_ * x_d_center_);
    }

    // Dynamics with 1st order state and control derivatives
    void dynamics(const Vector & x, const Vector & u, Matrix & A, Matrix & B) {
        A = A_;
        B = B_;
    }

    // Cost at sample k
    double cost(const Vector & x, const Vector & u, const int & k) {
        Vector e_track, e_obst;
        e_track.resize(x.rows()/3);
        e_obst.resize(x.rows()/3);
        e_track << x.block (x.rows()/3,0, x.rows()/3, 1);
        e_obst << x.tail(x.rows()/3);

        // std::cout << "Q_track: " << Q_track_ << "\n e_track: " << e_track  <<"\n expression: " << (e_track.transpose() * Q_track_ * e_track ) + (u.transpose() * (R_ * u)) << std::endl;
        return ((e_track.transpose() * Q_track_ * e_track ) + (u.transpose() * (R_ * u))) (0) + (f_ * exp(- 0.5 * e_obst.transpose() * W_ * e_obst));
    }

    // Cost at sample T
    double finalCost(const Vector & x, const Vector & u) {
        Vector e_track, e_obst;
        e_track.resize(x.rows()/3);
        e_obst.resize(x.rows()/3);
        e_track << x.block (x.rows()/3,0, x.rows()/3, 1);
        e_obst << x.tail(x.rows()/3);
        return ((e_track.transpose() * Q_track_f_ * e_track )) (0) + (f_ * exp(- 0.5 * e_obst.transpose() * W_ * e_obst));
    }


    // Cost at sample k with 2nd order state and control derivatives
    double cost(const Vector & x, const Vector & u, const int & k, double & q_0, Vector & q, Vector & r, Matrix & Q, Matrix & R) {
        Vector e_track, e_obst;
        int x_dim = x.rows();
        e_track.resize(x.rows()/3);
        e_obst.resize(x.rows()/3);
        e_track << x.block (x_dim/3,0, x_dim/3, 1);
        e_obst << x.tail(x_dim/3);
        R = 2.0*R_;
        Q.setZero();
        q.setZero();
        Q.block (x_dim/3,x_dim/3,x_dim/3,x_dim/3) = 2.0*Q_track_;
        q.block(x_dim/3,0,x_dim/3,1) << 2.0*Q_track_*e_track;
        double l_obs = (f_ * exp(- 0.5 * e_obst.transpose() * W_ * e_obst));
        q.tail(x_dim/3) << - W_ * l_obs * e_obst;
        Q.bottomRightCorner(x_dim/3,x_dim/3) << W_ * W_ * l_obs*(e_obst * e_obst.transpose()) - l_obs * W_ * Matrix::Identity(x_dim/3,x_dim/3);
        r = 2 * R_ * u;
        q_0 = cost(x,u,k);
        //std::cout << "Q: " << Q  <<"\n q: " << q << "\n q0: " << q_0 << "\n e_track: " << e_track << "\n e_obst: " << e_obst << std::endl;
        return q_0;
    }

    // Cost at sample T with 2nd order state derivatives
    double finalCost(const Vector & x, const Vector & u, double & q_0, Vector & q, Matrix & Q) {
        Vector e_track, e_obst;
        e_track.resize(x.rows()/3);
        e_obst.resize(x.rows()/3);
        int x_dim = x.rows();
        e_track << x.block (x_dim/3,0, x_dim/3, 1);
        e_obst << x.tail(x_dim/3);
        Q.setZero();
        q.setZero();
        Q.block (x_dim/3,x_dim/3,x_dim/3,x_dim/3) = 2.0*Q_track_f_;
        q.block(x_dim/3,0,x_dim/3,1) << 2.0*Q_track_f_*e_track;
        double l_obs = (f_ * exp(- 0.5 * e_obst.transpose() * W_ * e_obst));
        q.tail(x_dim/3) << - W_ * l_obs * e_obst;
        Q.bottomRightCorner(x_dim/3,x_dim/3) << W_ * W_ * l_obs*(e_obst * e_obst.transpose()) - l_obs * W_ * Matrix::Identity(x_dim/3,x_dim/3);
        q_0 = finalCost(x,u);
        //std::cout << "Q_f: " << Q << "\n Q_track_f: " << Q_track_f_ << "\n e_track: " << e_track  <<"\n q: " << q << "\n q0: " << q_0 << std::endl;
        return q_0;
    }

    void linear_dynamics(const Vector & x, const Vector & u, const int & k , Matrix & A, Matrix & B, std::vector<Matrix> & Sigma) {
        A = A_;
        B = B_;
        Sigma = Sigma_[k];
        //std::cout << "A: " << A << "\n B: " << B << std::endl;
    }
};

}
#endif /* POINTMASSLINEAROBSTACLE_H_ */
