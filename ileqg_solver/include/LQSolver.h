#ifndef LQSOLVER_H_
#define LQSOLVER_H_

#include "types.h"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <math.h>

namespace ileqg {

struct LQProblem {
    // Initial state
    Vector x0;

    // Dynamics
    std::vector<Matrix> A, B;
    std::vector< std::vector<Matrix> > Sigma;
    // Cost
    std::vector<Matrix> Q, R;
    std::vector<Vector> q, r;
    std::vector<double> q0;

    LQProblem(int T, int x_dim, int u_dim, int s_dim) {
        // resize variables
        A.resize(T);
        B.resize(T);
        Sigma.resize(T);
        Q.resize(T);
        R.resize(T);
        q.resize(T);
        r.resize(T);
        q0.resize(T);
        x0.resize(x_dim);
        x0.setZero();
        for (int i=0; i<T; i++) {
            A[i].resize(x_dim, x_dim);
            A[i].setZero();
            B[i].resize(x_dim, u_dim);
            B[i].setZero();
            Sigma[i].resize(s_dim);
             // x_dim,x_dim
            for (int j = 0 ; j < s_dim; j++) {
                Sigma[i][j].resize(x_dim,x_dim);
                Sigma[i][j].setZero();
            }
            Q[i].resize(x_dim, x_dim);
            Q[i].setZero();
            R[i].resize(x_dim, u_dim);
            R[i].setZero();
            q[i].resize(x_dim);
            q[i].setZero();
            r[i].resize(u_dim);
            r[i].setZero();
            q0[i] = 0.0;
        }
    }
};

struct LQSolution {
    // Optimal linear feedback solution
    std::vector<Vector> ff;
    std::vector<Vector> x;
    std::vector<Matrix> K;

    LQSolution(int T, int x_dim, int u_dim) {
        ff.resize(T);
        K.resize(T);
        x.resize(T);

        for (int i=0; i<T; i++) {
            ff[i].resize(u_dim);
            ff[i].setZero();
            K[i].resize(u_dim, x_dim);
            K[i].setZero();
            x[i].resize(x_dim);
            x[i].setZero();
        }
    }
};


class LQSolver
{

public :
    static double solve(const LQProblem & prob, LQSolution & sol, const std::vector<double> & theta) {
        return solve(prob.A.size(), prob.A, prob.B, prob.Sigma, prob.Q, prob.R, prob.q, prob.r, prob.q0, sol.ff, sol.K, theta, true);
    }

    static double getCost(const LQProblem & prob, LQSolution & sol, const std::vector<double> & theta) {
        std::cout << "getCost" << std::endl;
        return solve(prob.A.size(), prob.A, prob.B, prob.Sigma, prob.Q, prob.R, prob.q, prob.r, prob.q0, sol.ff, sol.K, theta, false);
    }

    // SOLVER: returns expected cost
    static double solve(int time_horizon, const std::vector<Matrix>& A, const std::vector<Matrix>& B, const std::vector< std::vector<Matrix> >& Sigma,
                        const std::vector<Matrix>& Q, const std::vector<Matrix>& R, const std::vector<Vector>& q, const std::vector<Vector>& r,
                        const std::vector<double>& q0, std::vector<Vector>& ff, std::vector<Matrix>& K, const std::vector<double> theta, bool optimize) {
        size_t x_dim = A[0].rows();
        size_t u_dim = B[0].cols();
        size_t s_dim = Sigma[0].size();
        /*
        Matrix W(x_dim, x_dim);
        Vector w(x_dim);
        double w0;

        W << Q[time_horizon-1];
        w << q[time_horizon-1];
        w0 = q0[time_horizon-1];
        */

        // Cost to go function = x^T W x + x^T w + w_0

        std::vector<Matrix> W, W_s;
        std::vector<Vector> w, w_s;
        std::vector <double> w0,w0_s;

        W.resize (s_dim);
        W_s.resize (s_dim);
        w.resize (s_dim);
        w_s.resize (s_dim);
        w0.resize (s_dim);
        w0_s.resize (s_dim);
        //initialize all
        for (size_t i = 0; i < s_dim; i++) {
            W[i] = Matrix(x_dim,x_dim);
            W_s[i] = Matrix(x_dim,x_dim);
            w[i] = Vector(x_dim);
            w_s[i] = Vector(x_dim);
        }

        for (size_t i = 0; i < s_dim; i++) {
            W[i] << Q[time_horizon-1];
            w[i] << q[time_horizon-1];
            w0[i] = q0[time_horizon-1];
        }

        // Shortcuts
        Vector g(u_dim);
        Matrix G(u_dim, x_dim);
        Matrix H(u_dim, u_dim);
        Matrix H_inverse(u_dim, u_dim);


       Matrix tmp = Matrix(x_dim,x_dim);
       Matrix invTerm = Matrix(x_dim,x_dim);

        // Bellman recursion
        for (int i=time_horizon-2 ; i >= 0 ; i --) {
            for (size_t j = 0; j < s_dim; j++) {
                W_s[j] =  W[j];
                w_s [j] = w [j];
                w0_s [j] = w0 [j];

                //invTerm << (Sigma[i][j].inverse() -  theta[j] * W[j]);
                //invTerm.inverse ();

                tmp << Eigen::MatrixXd::Identity(x_dim,x_dim) - theta[j] * W[j] * Sigma[i][j]; // sigma;
                double det = tmp.determinant();

                //W_s[j] = W_s[j] + theta[j] * W[j] *invTerm * W[j];
                W_s[j] = tmp.inverse() * W_s[j];
                //w_s[j] = w_s[j] + theta[j] * W[j] *invTerm * w[j];
                w_s[j] = tmp.inverse() * w_s[j];

                if (theta[j] != 0.0 ) {
                     //w0_s[j] = w0_s[j] + theta[j] * w[j].transpose() *invTerm * w[j] - 0.5 * log (det) /theta[j];
                    w0_s[j] = w0_s[j] + theta[j] * w[j].transpose() *Sigma[i][j] * w_s[j] - 0.5 * log (det) /theta[j];
                }
            }

            Matrix W_s_av (x_dim,x_dim);
            W_s_av << Eigen::MatrixXd::Zero (x_dim,x_dim);
            Vector w_s_av (x_dim);
            w_s_av << Eigen::VectorXd::Zero(x_dim);

            for (size_t j = 0; j < s_dim; j++) {
                W_s_av += W_s[j];
                w_s_av += w_s[j];
            }

            W_s_av = (1.0/s_dim) * W_s_av;
            w_s_av = (1.0/s_dim) * w_s_av;

            // Update shortcuts
            g = r[i] + B[i].transpose() * w_s_av;
            G = B[i].transpose()*W_s_av*A[i];
            H = R[i] + B[i].transpose()*W_s_av*B[i];

            if (optimize) {
                // Compute optimal control
                Eigen::EigenSolver<Matrix> es(H, true);
                Vector eigenvals = es.eigenvalues().real();

                //std::cout << "Eigenvalues: " << eigenvals << std::endl;

                // Check for negative (or too small) eigenvalues and correct them
                for (size_t e=1 ; e < u_dim ; e++) {
                    if (eigenvals(e) < 1e-15) {
                        eigenvals(e) = 1e-15;
                    }
                }

                H_inverse = es.eigenvectors().real() * eigenvals.asDiagonal().inverse()
                        * es.eigenvectors().real().transpose();

               ff[i] = - H_inverse * g;
               K[i] = - H_inverse * G;
            }

           for (size_t j = 0; j < s_dim; j++) {
               g = r[i] + B[i].transpose() * w_s[j];
               G = B[i].transpose()*W_s[j]*A[i];
               H = R[i] + B[i].transpose()*W_s[j]*B[i];

               W[j] = Q[i] + A[i].transpose()*W_s[j]*A[i] + K[i].transpose()*H*K[i] + K[i].transpose()*G + G.transpose()*K[i];
               w[j] = q[i] + A[i].transpose()*w_s[j] + K[i].transpose()*H*ff[i] + K[i].transpose()*g + G.transpose()*ff[i];
               w0[j] = q0[i] + w0_s[j] + 0.5*ff[i].transpose()*H*ff[i] + ff[i].transpose()*g;
           }
        }

        double cost= 0;

        for (size_t j = 0; j < s_dim; j++) {
          cost = cost + w0[j];
        }
        cost = (1.0/s_dim)*cost;

        return cost;
    }
};

}


#endif /* LQSOLVER_H_ */
