#include "ros/ros.h"
#include "ileqg_solver/FF_FB_plan.h"
#include "geometry_msgs/PoseStamped.h"
#include "eigen_conversions/eigen_msg.h"
#include <std_msgs/Float64MultiArray.h>
#include <sensor_msgs/JointState.h>
#include <IterativeLQSolver.h>
#include <PointMassLinearObstacle.h>

#include <tf/tf.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h>

#include <eigen3/Eigen/Dense>
#include <nav_msgs/Path.h>

/**
 *  Linear Quadratic MPC example
 *
 */

// This should be in a class
tf::Pose  pose;
geometry_msgs::Twist x_dot;
bool jp_msg_recieved = false;

void pose_callback(const geometry_msgs::PoseConstPtr &msg){
  jp_msg_recieved = true;
  pose.setOrigin(tf::Vector3(msg->position.x,msg->position.y,msg->position.z));
  pose.setRotation(tf::Quaternion(msg->orientation.x,
                                  msg->orientation.y,
                                  msg->orientation.z,
                                  msg->orientation.w));
}

void velocity_callback(const geometry_msgs::TwistConstPtr &msg){
    x_dot.linear.x = msg->linear.x;
    x_dot.linear.y = msg->linear.y;
    x_dot.linear.z = msg->linear.z;
    x_dot.angular.x = msg->angular.x;
    x_dot.angular.y = msg->angular.y;
    x_dot.angular.z = msg->angular.z;
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "lwr_linear_obst_node");
  ros::NodeHandle nh("lwr_linear_obst_node");

  /* Node parameters */
  double control_rate = 100;
  ros::Rate loop_rate(control_rate);


  /* Publishers and subscribers */
  ros::Publisher ff_fb_pub  = nh.advertise<ileqg_solver::FF_FB_plan>(
                               "/ff_fb_plan", 1);

  ros::Publisher viz_pub_rob  = nh.advertise<nav_msgs::Path>(
                              "/ileqg_rob", 10);

  ros::Publisher viz_pub_obst  = nh.advertise<nav_msgs::Path>(
                              "/ileqg_obst", 10);

  ros::Publisher viz_pub_goal  = nh.advertise<nav_msgs::Path>(
                              "/ileqg_goal", 10);

  tf::TransformListener tf_listener;

  /* Iterative solver parameters */
  unsigned int time_horizon = 250;
  unsigned int  x_dim = 6, u_dim = 3; // 3d position + velocity, 3d force
  unsigned int s_dim = 2;                                 // at the end effector
  double sampling_time = 1.0/50.0;

  // IterativeLQSolverParams( minLineSearchParameter,
  //                          initialLineSearchParameter,
  //                          minRelativeCostImprovement,
  //                          maxNumberIterations )

  std::vector<double> theta (s_dim);
  //theta = {-1e2,1e-5};
  //theta = {-1e-20,1e-19};
  theta = {-1e-1,1e-5};
  //theta = {0,0};
  /// Inizialize
  std::vector< std::vector< ileqg::Matrix > > Sigma;

  Sigma.resize(time_horizon);

  for (size_t i=0; i<Sigma.size(); i++) {
      Sigma[i].resize(s_dim);
       // x_dim,x_dim
       for (size_t j = 0 ; j < s_dim; j++) {
           // 0 variance. the eps=1e-6 is there to make it invertible
           Sigma[i][j] = Eigen::MatrixXd::Identity(x_dim*3,x_dim*3)*1e-6;
       }
       // Goal uncertainty
       Sigma[i][0].block(x_dim,x_dim, x_dim,x_dim) << Eigen::MatrixXd::Identity(x_dim,x_dim);
       // Obstacle uncertainty
       Sigma[i][1].bottomRightCorner(x_dim,x_dim) << Eigen::MatrixXd::Identity(x_dim,x_dim);
  }


  ileqg::IterativeLQSolverParams sparams(1e-30, 1, 1e-5, 1);
  ros::Time solution_time;

  /* Allocate/Initialize message variables */
  ileqg_solver::FF_FB_plan ff_fb_plan;
  ff_fb_plan.times.resize(time_horizon);
  ff_fb_plan.ff.resize(time_horizon);
  ff_fb_plan.fb.resize(time_horizon);
  ff_fb_plan.xd.resize(time_horizon);
  ff_fb_plan.xd_dot.resize(time_horizon);
  ff_fb_plan.x_external.resize(time_horizon);
  ff_fb_plan.x_external_dot.resize(time_horizon);
  ileqg::Matrix fb_matrix(6,18);
  fb_matrix.setZero();

  nav_msgs::Path viz_path;
  viz_path.poses.resize(time_horizon);
  nav_msgs::Path viz_path_obst;
  viz_path_obst.poses.resize(time_horizon);
  nav_msgs::Path viz_path_goal;
  viz_path_goal.poses.resize(time_horizon);


  /**             PROBLEM DEFINITION:
   *
   *  - Robot dynamics (3D cartesian admittance): M x_r_ddot + D x_r_dot = u
   *
   *          Robot state (2nd order) -> x_rob =  [x_r^T x_r_dot^T]^T;
   *              --> Linear dynamics    x_rob_dot = A_rob x_rob + B_rob u;
   *
   *  - Goal dynamics: x_d_dot = A_d(x_d_center - x_rob)
   *
   *  - Cost:  e^T Q_track_f e + sum_0^{time_horizon} e^T Q_track e + u^T R u
   *
   *               where e = (x_d- x_rob)
   *
   *    -- For convenience, we solve the problem with an augmented state
   *       x = [x_rob^T e^T]^T --
   *
   * - Initial state: x_initial = [x_rob^T  0^T]^T
   *
   **/

  /*** DYNAMICS ***/
  // Robot dynamic parameters M x_r_ddot + D x_r_dot = u
  ileqg::Matrix mass(x_dim/2, x_dim/2), damping(x_dim/2,x_dim/2),
                  massInverse(x_dim/2, x_dim/2);
  mass = ileqg::Matrix::Identity(x_dim/2, x_dim/2)* 10;
  damping = ileqg::Matrix::Identity(x_dim/2, x_dim/2)* 30;
  massInverse = mass.inverse();

  // Robot dynamics x_rob_dot = A_rob x_rob + B_rob u;
  ileqg::Matrix A_rob(x_dim,x_dim), B_rob(x_dim,u_dim), A_d(x_dim,x_dim);
  A_rob.setZero();
  A_rob.topRightCorner(x_dim/2, x_dim/2).setIdentity();
  A_rob.bottomRightCorner(x_dim/2, x_dim/2) << - massInverse * damping;
  B_rob << ileqg::Matrix::Zero(x_dim/2, x_dim/2) , massInverse;

  // Goal dynamics  x_d_dot = A_d(x_d_center - x_rob)
  A_d.setZero();
  A_d.topLeftCorner(x_dim/2,x_dim/2) << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  ileqg::Vector x_d_center(x_dim), x_obstacle_center(x_dim);

  // Goal attractor
  x_d_center << 1.14, 0.15, 1.24, 0, 0, 0;
  // Obstacle attractor
  x_obstacle_center << 1.16, -0.28, 1.14 ,0 ,0, 0;


  /*** COST ***/
  // Cost parameters
  ileqg::Matrix Q_track(x_dim, x_dim),
                Q_track_f(x_dim, x_dim),
                R(u_dim, u_dim);
  ileqg::Matrix W(x_dim,x_dim);
  W.setZero();
  W.topLeftCorner(x_dim/2,x_dim/2) << 1*ileqg::Matrix::Identity(x_dim/2, x_dim/2);
  Q_track_f = ileqg::Matrix::Identity(x_dim, x_dim) * 1e4;
  Q_track = (1.0 / sampling_time) * ileqg::Matrix::Identity(x_dim, x_dim) * 1e4;
  Q_track.bottomRightCorner(x_dim/2,x_dim/2).setZero();
  R = ileqg::Matrix::Identity(u_dim, u_dim) * 1e2;

  /*** INITIAL STATE ***/
  // Note that the error part of the augmented state
  // x = [x_rob^T e^T]^T at the beginning is always 0
  ileqg::Vector x_initial(x_dim*3);

  /***PROBLEM AND SOLVER***/
  // Define problem
  std::shared_ptr<ileqg::PointMassLinearObstacle> pointMassLinearProblem(
    new ileqg::PointMassLinearObstacle(time_horizon, A_rob, B_rob, A_d, A_d,
                        x_d_center, x_obstacle_center, Q_track, Q_track_f, R, Sigma, W, 6e4));


  // Initial nominal trajectory (0 trajectory)
  ileqg::LQSolution lqsol(time_horizon, x_dim*3, u_dim);
  // Define iterative solver
  ROS_INFO ("Before");
  ileqg::IterativeLQSolver iterative_solver(time_horizon, x_dim*3, u_dim,
                                            sampling_time,
                                            pointMassLinearProblem,
                                            sparams, x_initial, lqsol, theta);

  tf::StampedTransform transform;
  tf::StampedTransform transform_obst;
    try{
    ros::Time now = ros::Time::now();
    tf_listener.waitForTransform("/base_link", "/lwr_right_7_link",now,ros::Duration(5.0));
  }
  catch (tf::TransformException ex){
    ROS_ERROR("%s",ex.what());
    ros::Duration(1.0).sleep();
  }


  try{
  ros::Time now = ros::Time::now();
  tf_listener.lookupTransform("/base_link", "/lwr_right_7_link",
                           ros::Time(0), transform);

  tf_listener.lookupTransform("/base_link", "/obstacle", ros::Time(0), transform_obst);

}
catch (tf::TransformException ex){
  ROS_ERROR("%s",ex.what());
  ros::Duration(1.0).sleep();
}

/* Set initial state */
x_initial.setZero();
// Robot
x_initial.head(x_dim) << transform.getOrigin().getX(),
                         transform.getOrigin().getY(),
                         transform.getOrigin().getZ();
// Tracking error
x_initial.block(x_dim,0,x_dim,1) << 0,0,0;
// Difference to obstacle
x_initial.tail(x_dim) << transform_obst.getOrigin().getX() - transform.getOrigin().getX(),
                         transform_obst.getOrigin().getY() - transform.getOrigin().getY(),
                         transform_obst.getOrigin().getZ() - transform.getOrigin().getZ();

iterative_solver.setInitialState(x_initial);

// Get time
solution_time = ros::Time::now();
// Solve problem
iterative_solver.solve(lqsol);
// Populate ff_fb_plan message
for (size_t i=0 ; i < time_horizon ; i++) {
  // Time stamp
  ff_fb_plan.times[i].data = solution_time +
                             ros::Duration((double)(i*sampling_time));
  ff_fb_plan.ff[i].force.x = lqsol.ff[i](0);
  ff_fb_plan.ff[i].force.y = lqsol.ff[i](1);
  ff_fb_plan.ff[i].force.z = lqsol.ff[i](2);
  ff_fb_plan.ff[i].torque.x = 0;
  ff_fb_plan.ff[i].torque.y = 0;
  ff_fb_plan.ff[i].torque.z = 0;
  ff_fb_plan.xd[i].position.x = lqsol.x[i](0);
  ff_fb_plan.xd[i].position.y = lqsol.x[i](1);
  ff_fb_plan.xd[i].position.z = lqsol.x[i](2);
  ff_fb_plan.xd[i].orientation.x = 0;
  ff_fb_plan.xd[i].orientation.y = 0;
  ff_fb_plan.xd[i].orientation.z = 0;
  ff_fb_plan.xd[i].orientation.w = 1;
  ff_fb_plan.xd_dot[i].linear.x = lqsol.x[i](3);
  ff_fb_plan.xd_dot[i].linear.y = lqsol.x[i](4);
  ff_fb_plan.xd_dot[i].linear.z = lqsol.x[i](5);
  ff_fb_plan.xd_dot[i].angular.x = 0;
  ff_fb_plan.xd_dot[i].angular.y = 0;
  ff_fb_plan.xd_dot[i].angular.z = 0;
  ff_fb_plan.x_external[i].position.x = lqsol.x[i](12);
  ff_fb_plan.x_external[i].position.y = lqsol.x[i](13);
  ff_fb_plan.x_external[i].position.z = lqsol.x[i](14);
  ff_fb_plan.x_external[i].orientation.x = 0;
  ff_fb_plan.x_external[i].orientation.y = 0;
  ff_fb_plan.x_external[i].orientation.z = 0;
  ff_fb_plan.x_external[i].orientation.w = 1;
  ff_fb_plan.x_external_dot[i].linear.x = lqsol.x[i](15);
  ff_fb_plan.x_external_dot[i].linear.y = lqsol.x[i](16);
  ff_fb_plan.x_external_dot[i].linear.z = lqsol.x[i](17);
  ff_fb_plan.x_external_dot[i].angular.x = 0;
  ff_fb_plan.x_external_dot[i].angular.y = 0;
  ff_fb_plan.x_external_dot[i].angular.z = 0;

  fb_matrix.block<3,3>(0,0) << lqsol.K[i].block<3,3>(0,0);
  fb_matrix.block<3,3>(0,6) << lqsol.K[i].block<3,3>(0,6);
  fb_matrix.block<3,3>(0,12) << lqsol.K[i].block<3,3>(0,12);
  tf::matrixEigenToMsg(fb_matrix, ff_fb_plan.fb[i]);

  //std::cout  << "Feedback: "  << std::endl << lqsol.K[i] << std::endl;

  // For visualization purposes
  viz_path.poses[i].pose.position.x = lqsol.x[i](0);
  viz_path.poses[i].pose.position.y = lqsol.x[i](1);
  viz_path.poses[i].pose.position.z = lqsol.x[i](2);

  viz_path_goal.poses[i].pose.position.x = lqsol.x[i](6) + lqsol.x[i](0);
  viz_path_goal.poses[i].pose.position.y = lqsol.x[i](7) + lqsol.x[i](1) ;
  viz_path_goal.poses[i].pose.position.z = lqsol.x[i](8) + lqsol.x[i](2);

  viz_path_obst.poses[i].pose.position.x = lqsol.x[i](12) + lqsol.x[i](0);
  viz_path_obst.poses[i].pose.position.y = lqsol.x[i](13) + lqsol.x[i](1);
  viz_path_obst.poses[i].pose.position.z = lqsol.x[i](14) + lqsol.x[i](2);
  }


    ff_fb_pub.publish(ff_fb_plan);

    viz_path.header.frame_id = "/base_link";
    viz_pub_rob.publish(viz_path);

    viz_path_obst.header.frame_id = "/base_link";
    viz_pub_obst.publish(viz_path_obst);

    viz_path_goal.header.frame_id = "/base_link";
    viz_pub_goal.publish(viz_path_goal);

    ros::spinOnce();
    loop_rate.sleep();

  return 0;
}
