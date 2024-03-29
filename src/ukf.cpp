#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  is_initialized_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  previous_timestamp_ = 0;

  timeStep_ = 1;


  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.30;   
  
  // Process noise standard deviation yaw acceleration in rad/s^2

  std_yawdd_ = M_PI/8.0; 
  
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;    
  weights_ = VectorXd(2*n_aug_+1);
  
  Xsig_pred_ = MatrixXd(n_x_, n_aug_*2+1);

  NIS_radar_count_ = 0;
  NIS_radar_total_ = 0;
  NIS_laser_count_ = 0;
  NIS_laser_total_ = 0;
  NIS_radar_threshold_ = 7.815;
  NIS_laser_threshold_ = 5.991;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "UKF: " << endl;

        P_ << std_radr_*std_radr_, 0, 0, 0, 0,
              0, std_radr_*std_radr_, 0, 0, 0,
              0, 0, 0.9*0.9, 0, 0, 
              0, 0, 0, std_radphi_*std_radphi_, 0,
              0, 0, 0, 0, std_radphi_*std_radphi_;

        cout << " Time Step = " << timeStep_++ << endl;
        if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
          /**
          Convert radar from polar to cartesian coordinates and initialize state.
          */
          
          cout << "***** RADAR (initializing) ******" << endl;
          double rho = measurement_package.raw_measurements_[0];
          double phi = measurement_package.raw_measurements_[1];
          
          phi = atan2(sin(phi), cos(phi));

          x_ << rho * cos(phi), rho * sin(phi), 5.5, 0, 0; 


          cout << "x_ initialized = " << endl << x_ << endl;
          cout << "P_ initialized = " << endl << P_ << endl;

        }

        else if (measurement_package.sensor_type_ == MeasurementPackage::LASER) {
          /**
          Initialize state.
          */
          
          cout << "***** LASER (initializing) *****" << endl;
          x_ << measurement_package.raw_measurements_[0], 
                      measurement_package.raw_measurements_[1], 0, 0, 0;
          cout << "x_ initialized = " << endl << x_ << endl;
          cout << "P_ initialized = " << endl << P_ << endl;
        
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
  }

  if (previous_timestamp_==0)
      previous_timestamp_ = measurement_package.timestamp_;
  double dt = (double)(measurement_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  
  previous_timestamp_ = measurement_package.timestamp_;

  cout << endl << "**********************************************************************************" <<endl;
  cout << " Time Step = " << timeStep_++ << endl;
  //cout << endl << "dt : " << dt << endl;

  Prediction(dt);

  if (measurement_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true) {
      // Radar updates
      cout << endl << "***** RADAR update *****" << endl;

      UpdateRadar(measurement_package);

  } 
  else if (measurement_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true) {
    // Laser updates
    cout << endl << "***** LASER update *****" << endl;

    UpdateLidar(measurement_package);
  }

  // print the output
  cout << endl << "x_ updated = " << endl << x_ << endl;
  cout << endl << "P_ updated = " << endl << P_ << endl;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  VectorXd x_aug = VectorXd(n_aug_);

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  //cout << "Xsig_aug = " << endl << Xsig_aug << endl;

  //Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p +=  0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p += 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p += nu_a*delta_t;

    yaw_p += 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p += nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //cout << endl << "Xsig_pred = " << endl << Xsig_pred_ << std::endl;

   // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // cout << endl << "Weights : " << endl << weights_ << endl;

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  cout << endl << "x_ predicted = " << endl << x_ << endl;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //cout << "x_diff(3) before normalization: " << endl << x_diff(3) << endl;
    // angle normalization
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3) ) );
  //cout << "x_diff(3) after normalization: " << endl << x_diff(3) << endl;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }


  //print result
  
  std::cout << endl << "P_ predicted = " << endl << P_ << std::endl;



}



/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    
    // measurement model
    Zsig(0,i) = p_x;                       
    Zsig(1,i) = p_y;                                 
  }

  //cout << "Zsig = " << endl << Zsig << endl;

  VectorXd z = VectorXd(n_z);
  z << measurement_package.raw_measurements_[0], measurement_package.raw_measurements_[1];
  
  MatrixXd R = MatrixXd(n_z, n_z);
  R <<  std_laspx_*std_laspx_,0,
        0, std_laspy_*std_laspy_;
  Update(Zsig, z, R);
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} measurement_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_package) {
  
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  VectorXd z = VectorXd(n_z);
  z << measurement_package.raw_measurements_[0], measurement_package.raw_measurements_[1], measurement_package.raw_measurements_[2];

   //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  
  Update(Zsig, z, R);
}

/**
 *  Updates the state and the state covariance matrix using a passed-in Z Sigma measurements, the actual measurement, 
 *   and the noise covariance, from either radar or lidar.
 *  @param {const MatrixXd&} Zsig, {const VectorXd&} z, {const MatrixXd&} R
 */
void UKF::Update(const MatrixXd& Zsig, const VectorXd& z, const MatrixXd& R) {
  
  int n_z = z.size();
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //cout << "z_pred : " << endl << z_pred <<endl;

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  // cross correlation matrix Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    if (n_z==3) {
      //angle normalization  
      // while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      // while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1) ) );      
    }
    S += weights_(i) * z_diff * z_diff.transpose();

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_; 
    //angle normalization
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3) ) );
  

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
 
  S += R;

  //print result
  //cout << "z_pred : " << endl << z_pred << endl;
  //cout << "S : " << endl << S << endl;
  //cout << "Tc : " << endl << Tc << endl;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //residual
  VectorXd z_diff = z - z_pred;

  if (n_z==3)
  {
    //angle normalization
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1) ) );
  }

  double NIS = z_diff.transpose() * S.inverse() * z_diff;
  cout << endl << "NIS : " << NIS<<endl;
  if (n_z == 2) // Laser
  {
    NIS_laser_total_++;
    if (NIS > NIS_laser_threshold_)
      NIS_laser_count_++;
    cout << "NIS laser percentage over threshold: " << (float)NIS_laser_count_ / NIS_laser_total_ *100.0 << "%" << endl;
  }
  if (n_z == 3) // Radar
  {
    NIS_radar_total_++;
    if (NIS > NIS_radar_threshold_)
      NIS_radar_count_++;
    cout << "NIS radar percentage over threshold: " << (float)NIS_radar_count_ / NIS_radar_total_ *100.0 << "%" << endl;
  }  

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();

}

