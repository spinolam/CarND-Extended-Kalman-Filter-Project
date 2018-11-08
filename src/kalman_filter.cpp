#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter()
{
}

KalmanFilter::~KalmanFilter()
{
}

void KalmanFilter::Predict()
{
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  const MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  const VectorXd z_pred = H_ * x_;
  const VectorXd y = z - z_pred;
  const MatrixXd Ht = H_.transpose();
  const MatrixXd S = H_ * P_ * Ht + R_;
  const MatrixXd K = P_ * Ht * S.inverse();

  // new estimate
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}

double normalizeTheta(const double theta_in)
{
  double theta_norm = theta_in;
  
  while (theta_norm > M_PI)
  {
    theta_norm -= (2 * M_PI);
  }
  
  while (theta_norm < -M_PI)
  {
    theta_norm += (2 * M_PI);
  }

  return theta_norm; 
}

// z - actual measurement in polar coordinates
// h - predicted value x_ in polar coordinates (Hx)
void KalmanFilter::UpdateEKF(const VectorXd &z, const VectorXd &h)
{
  const MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  const MatrixXd PHt = P_ * H_.transpose();
  const MatrixXd S = H_ * PHt + R_;
  const MatrixXd K = PHt * S.inverse();
  VectorXd y = z - h; // dif between current radar measurement z and predicted x_ (in polar coordinates)
  y(1) = normalizeTheta(y(1));
  
  // new estimate
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}
