#ifndef TOOLS_H_
#define TOOLS_H_
#include "measurement_package.h"
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * Helper functions to extract sensor measurement reading
  */
  VectorXd GetRadarMeasurement(const MeasurementPackage &measurement);
  VectorXd GetLidarMeasurement(const MeasurementPackage &measurement);
  VectorXd MeasurementToCartesian(const MeasurementPackage &measurement);

  /**
  * A helper method to convert cartesian to polar coordinates.
  */
  VectorXd CartesianToPolar(const VectorXd &vec);

  /**
  * A helper method to convert polar to cartesian coordinates.
  */
  VectorXd PolarToCartesian(const VectorXd &vec);
  
  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

};

#endif /* TOOLS_H_ */
