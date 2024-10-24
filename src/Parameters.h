#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// Define a "Variable" that can store and return different data types
struct Variable {

  // Enumerate the admissible Variable value types
  enum ValueType {
    NONE_VARIABLE,
    DOUBLE_VARIABLE,
    INT_VARIABLE,
    BOOL_VARIABLE,
    STRING_VARIABLE,
    VECTOR_VARIABLE
  };

  // explicit conversion from specified type to Variable
  Variable(void) { type = ValueType::NONE_VARIABLE; }
  Variable(double              value) { double_value = value; type = ValueType::DOUBLE_VARIABLE; }
  Variable(int                 value) { int_value    = value; type = ValueType::INT_VARIABLE;    }
  Variable(bool                value) { bool_value   = value; type = ValueType::BOOL_VARIABLE;   }
  Variable(std::string         value) { string_value = value; type = ValueType::STRING_VARIABLE; }
  Variable(std::vector<double> value) { vector_value = value; type = ValueType::VECTOR_VARIABLE; }

  void check_type_match(ValueType requested) { if (type != requested) std::cerr << "ERROR in `Variable`; the requested value does not match the Variable's type" << std::endl; }

  // implicit conversion from Variable to specified return type
  operator double()              { check_type_match(ValueType::DOUBLE_VARIABLE); return double_value; }
  operator int()                 { check_type_match(ValueType::INT_VARIABLE);    return int_value;    }
  operator bool()                { check_type_match(ValueType::BOOL_VARIABLE);   return bool_value;   }
  operator std::string()         { check_type_match(ValueType::STRING_VARIABLE); return string_value; }
  operator std::vector<double>() { check_type_match(ValueType::VECTOR_VARIABLE); return vector_value; }

  ValueType           type;
  double              double_value;
  int                 int_value;
  bool                bool_value;
  std::string         string_value;
  std::vector<double> vector_value;
  
}; // Variable

// Define a "Parameters" object that stores Variable-valued key-value pairs
typedef std::unordered_map<std::string,Variable> Parameters;

#endif /* PARAMETERS_H */
