#include "cavaan.h"
#include <cctype>


Expression::Expression(const std::string& expr) {
  tokenize(expr);
  parseToRPN();
}


bool isalpha_(char x) {
  return isalpha(x) || x == '_';
}


void Expression::tokenize(const std::string& expr) {
  size_t i = 0;
  while (i < expr.length()) {
    if (isspace(expr[i])) {
      ++i;
    } else if (isdigit(expr[i]) || expr[i] == '.') {
      // Parse number
      size_t start = i;
      while (i < expr.length() && (isdigit(expr[i]) || expr[i] == '.')) ++i;
      tokens.push_back({ NUMBER, expr.substr(start, i - start) });
    } else if (isalpha_(expr[i])) {
      // Parse identifier (variable or function)
      size_t start = i;
      while (i < expr.length() && (isalpha_(expr[i]) || isdigit(expr[i]))) ++i;
      std::string ident = expr.substr(start, i - start);
      if (functions.find(ident) != functions.end()) {
        tokens.push_back({ FUNCTION, ident });
      } else {
        tokens.push_back({ VARIABLE, ident });
      }
    } else if (expr[i] == '(' || expr[i] == ')') {
      tokens.push_back({ PARENTHESIS, std::string(1, expr[i]) });
      ++i;
    } else {
      // Operators
      char op = expr[i];
      if (op == '+' || op == '-' || op == '*' || op == '/' || op == '^') {
        tokens.push_back({ OPERATOR, std::string(1, op) });
        ++i;
      } else {
        throw std::runtime_error(std::string("Invalid character: ") + expr[i]);
      }
    }
  }
}


int Expression::getPrecedence(const std::string& op) {
  if (op == "+" || op == "-")
    return 1;
  else if (op == "*" || op == "/")
    return 2;
  else if (op == "^")
    return 3;
  else
    return 0;
}


bool Expression::isRightAssociative(const std::string& op) {
  return op == "^";
}


void Expression::parseToRPN() {
  std::stack<Token> opStack;
  for (const Token& token : tokens) {
    if (token.type == NUMBER || token.type == VARIABLE) {
      rpn.push_back(token);
    } else if (token.type == FUNCTION) {
      opStack.push(token);
    } else if (token.type == OPERATOR) {
      while (!opStack.empty() && (opStack.top().type == OPERATOR || opStack.top().type == FUNCTION)) {
        Token o2 = opStack.top();
        if (o2.type == FUNCTION || (getPrecedence(token.value) < getPrecedence(o2.value)) ||
          (getPrecedence(token.value) == getPrecedence(o2.value) && !isRightAssociative(token.value))) {
          opStack.pop();
          rpn.push_back(o2);
        } else {
          break;
        }
      }
      opStack.push(token);
    } else if (token.type == PARENTHESIS) {
      if (token.value == "(") {
        opStack.push(token);
      } else if (token.value == ")") {
        while (!opStack.empty() && opStack.top().value != "(") {
          rpn.push_back(opStack.top());
          opStack.pop();
        }
        if (opStack.empty()) {
          throw std::runtime_error("Mismatched parentheses");
        }
        opStack.pop(); // Pop the '('
        if (!opStack.empty() && opStack.top().type == FUNCTION) {
          rpn.push_back(opStack.top());
          opStack.pop();
        }
      }
    }
  }
  while (!opStack.empty()) {
    if (opStack.top().value == "(" || opStack.top().value == ")") {
      throw std::runtime_error("Mismatched parentheses");
    }
    rpn.push_back(opStack.top());
    opStack.pop();
  }
}


double Expression::evaluateRPN(const std::map<std::string, double>& variables) {
  std::stack<double> evalStack;
  for (const Token& token : rpn) {
    if (token.type == NUMBER) {
      evalStack.push(std::stod(token.value));
    } else if (token.type == VARIABLE) {
      auto it = variables.find(token.value);
      if (it != variables.end()) {
        evalStack.push(it->second);
      } else {
        throw std::runtime_error("Undefined variable: " + token.value);
      }
    } else if (token.type == OPERATOR) {
      if (evalStack.size() < 2) {
        throw std::runtime_error("Insufficient values in expression");
      }
      double rhs = evalStack.top(); evalStack.pop();
      double lhs = evalStack.top(); evalStack.pop();
      if (token.value == "+") {
        evalStack.push(lhs + rhs);
      } else if (token.value == "-") {
        evalStack.push(lhs - rhs);
      } else if (token.value == "*") {
        evalStack.push(lhs * rhs);
      } else if (token.value == "/") {
        evalStack.push(lhs / rhs);
      } else if (token.value == "^") {
        evalStack.push(std::pow(lhs, rhs));
      }
    } else if (token.type == FUNCTION) {
      if (evalStack.empty()) {
        throw std::runtime_error("Insufficient values in expression");
      }
      double arg = evalStack.top(); evalStack.pop();
      if (token.value == "sin") {
        evalStack.push(std::sin(arg));
      } else if (token.value == "cos") {
        evalStack.push(std::cos(arg));
      } else if (token.value == "log") {
        evalStack.push(std::log(arg));
      } else if (token.value == "abs") {
        evalStack.push(std::abs(arg));
      } else if (token.value == "sqr") {
        evalStack.push(arg * arg);
      } else {
        throw std::runtime_error("Unknown function: " + token.value);
      }
    }
  }
  if (evalStack.size() != 1) {
    throw std::runtime_error("Invalid expression");
  }
  return evalStack.top();
}


double Expression::evaluate(const std::map<std::string, double>& variables) {
  return evaluateRPN(variables);
}


double Expression::evaluateR(const Rcpp::List& vars) {
  std::map<std::string, double> variables;
  Rcpp::CharacterVector names = vars.names();
  for (int i = 0; i < vars.size(); ++i) {
    std::string name = Rcpp::as<std::string>(names[i]);
    double value = Rcpp::as<double>(vars[i]);
    variables[name] = value;
  }
  return evaluate(variables);
}


Expression *createExpression(std::string expr) {
  Expression* expr_ptr = new Expression(expr);
  return expr_ptr;
}


double evaluateExpression(Expression *expr_ptr, Rcpp::List vars) {
  return expr_ptr->evaluateR(vars);
}


std::vector<std::string> Expression::getVariables() {
  std::vector<std::string> variables;

  for (int i = 0; i < (int)(tokens.size()); i++) {
    // enum TokenType VARIABLE_TYPE = VARIABLE;
    if (tokens[i].type != VARIABLE) continue;
    variables.push_back(tokens[i].value);
  }

  return variables;
}


// [[Rcpp::export]]
Rcpp::StringVector getVariablesEquation(std::string expr) {
  Expression *expr_ptr = createExpression(expr);
  std::vector<std::string> variables = expr_ptr->getVariables();
  Rcpp::StringVector out(variables.begin(), variables.end());
  return out;
}
