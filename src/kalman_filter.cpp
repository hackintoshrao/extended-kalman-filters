#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_ * x_ ;
  Matrix Ft = F_.trasponse();
  P_ = F_ * P_* Ft + Q_ ;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd H_pred(3);

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // Obtain the rho, phi and rho-dot (radars measurement space)
  // from the predicted state values.
  float rho_pred = sqrt(px * px + py * py);
  //check for division by zero
  if(rho_pred < .00001) {
   px += .001;
   py += .001;
   rho_pred = sqrt(px * px + py * py);
  }

  float phi_pred = atan2(py,px);
  // Reducing the phi_value to range -pi to +pi.
  while(phi_pred > M_PI) phi_pred -= M_PI;
  while(phi_pred < -M_PI) phi_pred += M_PI;
  
  float rho_dot_pred = (px*vx+py*vy)/ rho_pred;


  VectorXd H_func(3);

  H_pred  << rho_pred, phi_pred, rho_dot_pred;
 // obtain the error
  VectorXd  y = z - H_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

 //new estimate
 x_ = x_ + (K * y);
 long x_size = x_.size();
 MatrixXd I = MatrixXd::Identity(x_size, x_size);
 P_ = (I - K * H_) * P_;
}
