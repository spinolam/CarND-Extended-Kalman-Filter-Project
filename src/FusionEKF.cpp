#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

FusionEKF::FusionEKF()
{
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd::Zero(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 0, 0.0225;
  H_laser_ << 1.0, 0, 0, 0, 0, 1.0, 0, 0;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

  ekf_.x_ = VectorXd::Zero(4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000;
  ekf_.F_ = MatrixXd::Identity(4, 4);
  ekf_.F_(0, 2) = 1;
  ekf_.Q_ = MatrixXd::Zero(4, 4);
}

FusionEKF::~FusionEKF()
{
}

void FusionEKF::UpdateQ(const double dt, const double noise_ax, const double noise_ay)
{
  const double dt2 = dt * dt;
  const double dt3 = dt2 * dt;
  const double dt4 = dt3 * dt;

  ekf_.Q_(0, 0) = dt4 / 4.0 * noise_ax;
  ekf_.Q_(0, 2) = dt3 / 2.0 * noise_ax;
  ekf_.Q_(1, 1) = dt4 / 4.0 * noise_ay;
  ekf_.Q_(1, 3) = dt3 / 2.0 * noise_ay;
  ekf_.Q_(2, 0) = dt3 / 2.0 * noise_ax;
  ekf_.Q_(2, 2) = dt2 * noise_ax;
  ekf_.Q_(3, 1) = dt3 / 2.0 * noise_ay;
  ekf_.Q_(3, 3) = dt2 * noise_ay;
}

void FusionEKF::UpdateF(const double dt)
{
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  if (not is_initialized_)
  {
    ekf_.x_ = tools.MeasurementToCartesian(measurement_pack);
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }
  
  const float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  // dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  UpdateF(dt);
  UpdateQ(dt, 9, 9);
  ekf_.Predict();

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    VectorXd z = tools.GetRadarMeasurement(measurement_pack); // current measurement in polar
    VectorXd h = tools.CartesianToPolar(ekf_.x_); // predicted in polar
    ekf_.H_ = tools.CalculateJacobian(tools.PolarToCartesian(z)); // Jacobian of current measurement 
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(z, h);
  }
  else
  {
    VectorXd z = tools.GetLidarMeasurement(measurement_pack);
    ekf_.H_ = H_laser_; 
    ekf_.R_ = R_laser_;
    ekf_.Update(z);
  }

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
