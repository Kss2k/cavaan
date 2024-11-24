#ifndef model_h
#define model_h


#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <stack>
#include <map>
#include <cmath>
#include <stdexcept>
#include <set>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


#define BETA_STAR 0
#define GAMMA_STAR 1
#define PHI 2


class Expression {
public:
  Expression(const std::string& expr);

  double evaluate(const std::map<std::string, double>& variables);
  double evaluateR(const Rcpp::List& vars);
  std::vector<std::string> getVariables();

private:
  enum TokenType { NUMBER, VARIABLE, FUNCTION, OPERATOR, PARENTHESIS };

  struct Token {
    TokenType type;
    std::string value;
  };

  std::vector<Token> tokens; // Original tokens
  std::vector<Token> rpn;  // Reverse Polish Notation tokens

  void tokenize(const std::string& expr);
  void parseToRPN();
  double evaluateRPN(const std::map<std::string, double>& variables);

  int getPrecedence(const std::string& op);
  bool isRightAssociative(const std::string& op);


  std::set<std::string> functions = { "sin", "cos", "log", "abs", "sqr" };
};


// Matrices for a single group in a SEM
typedef struct {
  arma::mat BStar;
  arma::mat GammaStar;
  arma::mat Phi;
  arma::mat BStarInv;
  arma::mat G;
  arma::mat S;
  arma::mat Sigma;
  int p;
} MatricesGroup;


// ParTable
typedef struct {
  std::vector<std::string> lhs;
  std::vector<std::string> op;
  std::vector<std::string> rhs;
  arma::vec est;
  std::vector<std::string> label;
  std::vector<int>  matrix;
  std::vector<int>  group;
  std::vector<int>  row;
  std::vector<int>  col;
  std::vector<bool> free;
  std::vector<bool> fill;
  std::vector<bool> continueFromLast;
  std::vector<bool> isEquation;
  std::vector<Expression*> expressions;
  int nfree;
} ParTable;


// GradientMatricesParam
typedef struct {
  std::vector<MatricesGroup*> matricesGroups; // vector with elem for each group
} GradientMatricesParam;


typedef struct {
  int iterations;
  int maxIterations;
  double threshold;
  bool convergence;
  arma::vec par;
} OptimizerInfo;


// Model
typedef struct {
  std::vector<MatricesGroup*> matricesGroups;
  std::vector<GradientMatricesParam*> gradientMatricesParams; // vector with elem for each param
  ParTable *parTable;
  int ngroups;
  int p;
  OptimizerInfo optimizerInfo;
} Model;


// model.cpp
MatricesGroup *createMatricesGroup(Rcpp::List matrices);
MatricesGroup *copyMatricesGroup(MatricesGroup *matricesGroup, bool fillZero);
Model *createModel(Rcpp::List model);
ParTable *createParTable(Rcpp::List model);
void fillModel(Model *model, const arma::vec &theta, bool replace, bool calcSigma);
void fillMatricesGroups(std::vector<MatricesGroup*> matricesGroups, ParTable *parTable, 
    const arma::vec &theta, bool calcSigma, bool fillConst);
int countFree(std::vector<bool> free);


// gradient.cpp
void getBaseGradients(Model *model);
// arma::vec getGradientModel(const arma::vec &theta, Model *model);
arma::vec getGradientModelSimple(const arma::vec &theta, Model *model);
arma::vec getGradientModelGeneral(const arma::vec &theta, Model *model);
arma::vec normalize(arma::vec x);
double getLogLikModel(const arma::vec &theta, Model *model);


// parse_eqations.cpp
Expression *createExpression(std::string expr);
double evaluateExpression(Expression *expr_ptr, Rcpp::List vars);


#endif
