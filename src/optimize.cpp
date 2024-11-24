#include "cavaan.h"
#include <float.h>
#include <cmath>

arma::vec optim2(arma::vec theta, Model* model, double (*objective)(arma::vec, Model*),
    arma::vec (*gradient)(arma::vec, Model*)) {

  const double b1a = 1;
  double b1, b1b, bh, val_k, val_l = DBL_MAX; // l = k - 1
  arma::vec grad_b0, grad_b1, grad_bh, d, theta_b0 = theta, theta_b1;

  for (int k = 0; k < model->optimizerInfo.maxIterations; k++) {
    val_l   = val_k;
    val_k   = objective(theta_b0, model);
    grad_b0 = gradient(theta_b0, model);
    d       = (-1) * normalize(grad_b0); // d / sqrt(d'd)
 
    if (k && model->optimizerInfo.threshold > std::abs(val_k - val_l)) break;
    Rcpp::Rcout << "Iteration: " << k + 1 << ", LogLik: " << val_l << 
      ", Change: " << val_k - val_l << "\n";
    
    b1b = -2 * val_k / arma::dot(grad_b0, d);

    b1 = b1b < b1a ? b1b : b1a;

    theta_b1 = theta + b1 * d;
    grad_b1 = gradient(theta_b1, model);
   

    theta_b0 = grad_b0 * b1 / (grad_b1 - grad_b0);
  
  }
  
  return theta;
}


arma::vec optim(arma::vec theta, Model* model, double (*objective)(arma::vec, Model*),
                 arma::vec (*gradient)(arma::vec, Model*)) {

  const double max_step_size = 1.0; // Max allowable step size
  double val_k, val_l = DBL_MAX;   // Initialize objective values
  arma::vec grad_b0, d, theta_b0 = theta;

  for (int k = 0; k < model->optimizerInfo.maxIterations; k++) {
    // Evaluate objective and gradient at current point
    val_l = val_k;
    val_k = objective(theta_b0, model);
    grad_b0 = gradient(theta_b0, model);

    // Compute descent direction (negative normalized gradient)
    d = (-1) * grad_b0 / normalize(grad_b0);

    // Print iteration details
    Rcpp::Rcout << "Iteration: " << k + 1 << ", Objective: " << val_k << 
      ", Change: " << std::abs(val_k - val_l) << "\n";

    // Check convergence
    if (k > 0 && (model->optimizerInfo.threshold > std::abs(val_k - val_l))) {
      Rcpp::Rcout << "Finished: Iteration: " << k + 1 << ", Objective: " << val_k << 
        ", Change: " << std::abs(val_k - val_l) << "\n" << "Threshold: " <<
        model->optimizerInfo.threshold << "\n";
      break;
    }


    // Perform backtracking line search for step size
    double b1 = max_step_size;
    double alpha = 0.5;  // Reduction factor
    double beta = 0.1;   // Armijo condition constant

    while (objective(theta_b0 + b1 * d, model) > val_k + beta * b1 * arma::dot(grad_b0, d)) {
      b1 *= alpha; // Reduce step size
      if (b1 < 1e-10) { // Avoid too small step sizes
        Rcpp::Rcout << "Step size too small, stopping optimization.\n";
        return theta_b0;
      }
    }

    // Update parameters
    theta_b0 = theta_b0 + b1 * d;
  }

  return theta_b0;
}


arma::vec optimBFGS(arma::vec theta, Model* model, double (*objective)(arma::vec, Model*),
                    arma::vec (*gradient)(arma::vec, Model*)) {
    const int maxIterations = model->optimizerInfo.maxIterations;
    const double threshold = model->optimizerInfo.threshold;

    int n = theta.n_elem;
    arma::mat H_inv = arma::eye<arma::mat>(n, n); // Initialize inverse Hessian approximation as identity matrix
    arma::vec grad_k, grad_k1, s_k, y_k, theta_k1;
    double val_k, val_k1;

    // Evaluate initial objective and gradient
    val_k = objective(theta, model);
    grad_k = gradient(theta, model);

    for (int k = 0; k < maxIterations; k++) {
        // Compute the search direction
        arma::vec d = -H_inv * grad_k;

        // Perform a line search to find the step size
        double alpha = 1.0; // Initial step size
        double alpha_decay = 0.5; // Step size reduction factor
        double c1 = 1e-4; // Armijo condition parameter

        while (true) {
            theta_k1 = theta + alpha * d;
            val_k1 = objective(theta_k1, model);
            if (val_k1 <= val_k + c1 * alpha * arma::dot(grad_k, d)) {
                break; // Armijo condition satisfied
            }
            alpha *= alpha_decay;
            if (alpha < 1e-10) {
                Rcpp::Rcout << "Step size too small, stopping optimization.\n";
                return theta; // Return current parameters if line search fails
            }
        }

        // Compute the new gradient
        grad_k1 = gradient(theta_k1, model);

        // Check for convergence
        if (arma::norm(grad_k1) < threshold) {
            Rcpp::Rcout << "Converged after " << k + 1 << " iterations.\n";
            return theta_k1;
        }

        // Update the inverse Hessian approximation using the BFGS formula
        s_k = theta_k1 - theta; // Parameter change
        y_k = grad_k1 - grad_k; // Gradient change

        double rho = 1.0 / arma::dot(y_k, s_k); // Ensure denominator is not zero
        if (std::abs(rho) > 1e-10) { // Avoid numerical issues
            arma::mat outer_s = s_k * s_k.t();
            arma::mat outer_y = y_k * y_k.t();
            arma::mat H_y = H_inv * y_k;

            H_inv += rho * outer_s - rho * H_y * H_y.t() * rho;
        }

        // Update parameters and gradient for the next iteration
        theta = theta_k1;
        grad_k = grad_k1;
        val_k = val_k1;

        // Diagnostic output
        Rcpp::Rcout << "Iteration: " << k + 1 << ", Objective: " << val_k << ", Gradient Norm: " << arma::norm(grad_k) << "\n";
    }

    Rcpp::Rcout << "Reached maximum iterations without convergence.\n";
    return theta;
}
