#include <iostream>
#include <math.h>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools()
{
}

Tools::~Tools()
{
}

VectorXd Tools::GetRadarMeasurement(const MeasurementPackage &measurement)
{
  VectorXd m = VectorXd::Zero(3);
  if (measurement.sensor_type_ == MeasurementPackage::RADAR)
  {
    m << measurement.raw_measurements_[0], measurement.raw_measurements_[1], measurement.raw_measurements_[2];
  }
  else
  {
    std::cout << "Invalid measurement request!" << std::endl;
  }
  return m;
}

VectorXd Tools::GetLidarMeasurement(const MeasurementPackage &measurement)
{
  VectorXd m = VectorXd::Zero(2);
  if (measurement.sensor_type_ == MeasurementPackage::LASER)
  {
    m << measurement.raw_measurements_[0], measurement.raw_measurements_[1];
  }
  else
  {
    std::cout << "Invalid measurement request!" << std::endl;
  }
  return m;
}

VectorXd Tools::MeasurementToCartesian(const MeasurementPackage &measurement)
{
  VectorXd cartesian = VectorXd::Zero(4);

  if (measurement.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Convert radar from polar to cartesian coordinates and initialize state.
    VectorXd polar = VectorXd::Zero(3);
    polar << measurement.raw_measurements_[0], measurement.raw_measurements_[1], measurement.raw_measurements_[2];
    cartesian = PolarToCartesian(polar);
  }
  else if (measurement.sensor_type_ == MeasurementPackage::LASER)
  {
    cartesian << measurement.raw_measurements_[0], measurement.raw_measurements_[1], 0, 0;
  }
  else
  {
    // Invalid measurement type
    std::cout << "Invalid measurement package!" << std::endl;
  }

  return cartesian;
}

VectorXd Tools::CartesianToPolar(const VectorXd &vec)
{
  const double px = vec(0);
  const double py = vec(1);
  const double vx = vec(2);
  const double vy = vec(3);

  const double range = sqrt(px * px + py * py);
  const double theta = atan2(py, px); // returns value between -pi and pi
  const double range_rate = (px * vx + py * vy) / std::max(0.0001, range);
    
  VectorXd polar = VectorXd::Zero(3);
  polar << range, theta, range_rate;
  return polar;
}

VectorXd Tools::PolarToCartesian(const VectorXd &vec)
{
  VectorXd cartesian = VectorXd::Zero(4);

  const double range = vec(0);
  const double theta = vec(1);
  const double range_rate = vec(2);

  const double px = range * cos(theta);
  const double py = range * sin(theta);
  const double vx = range_rate * cos(theta);
  const double vy = range_rate * sin(theta);
  cartesian << px, py, vx, vy;
  return cartesian;
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth)
{
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() == 0 or ground_truth.size() == 0 or (estimations.size() != ground_truth.size()))
  {
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i)
  {
    VectorXd residual = ground_truth[i] - estimations[i];
    VectorXd sqr_residual = residual.array() * residual.array();
    rmse += sqr_residual;
  }

  // // calculate the mean
  rmse = rmse / estimations.size();

  // // calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
  MatrixXd Hj = MatrixXd::Zero(3, 4);

  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // compute the Jacobian matrix
  const double px2 = px * px;
  const double py2 = py * py;
  const double px2_plus_py2 = px2 + py2;
  const double sqr_px2_plus_py2 = sqrt(px2_plus_py2);

  // check division by zero
  if (fabs(px2_plus_py2) <= 0.0001)
  {
    cout << "divsion by zero!" << endl;
    return Hj;
  }

  Hj(0, 0) = px / sqr_px2_plus_py2;
  Hj(0, 1) = py / sqr_px2_plus_py2;
  Hj(1, 0) = -(py / px2_plus_py2);
  Hj(1, 1) = px / px2_plus_py2;
  Hj(2, 0) = py * (vx * py - vy * px) / (px2_plus_py2 * sqr_px2_plus_py2);
  Hj(2, 1) = px * (vy * px - vx * py) / (px2_plus_py2 * sqr_px2_plus_py2);
  Hj(2, 2) = px / sqr_px2_plus_py2;
  Hj(2, 3) = py / sqr_px2_plus_py2;

  return Hj;
}
