#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm> // For std::find_if, std::erase
#include <cctype>    // For std::isspace
#include "curveReader.h"

// Function to trim leading and trailing whitespace
std::string trim(std::string s) {
    // Trim leading whitespace
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));

    // Trim trailing whitespace
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
    return s;
}

bool parseBool(const std::string& str) {
    return (str == "1" || str == "true" || str == "True" || str == "TRUE");
}

template <typename T>
void printVec(std::string name, std::vector<T>& vec) {
  std::cerr << name << " "; 
  for( auto v : vec ) std::cerr << v << " ";
  std::cerr << "\n";
}

void CurveReader::printCurveInfo(CurveInfo& c) {
  printVec("x", c.x);
  printVec("y", c.y);
  printVec("isOnCurve", c.isOnCurve);
  printVec("isMdlVtx", c.isMdlVtx);
}

CurveReader::CurveInfo CurveReader::readCurveInfo(const std::string& filename) {
  CurveReader::CurveInfo c;
  std::ifstream file(filename);
  std::string line;

  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return c;
  }

  // Read header
  std::getline(file, line);

  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string token;
    double ignored;
    std::getline(ss, token, ','); c.x.push_back(std::stod(token));
    std::getline(ss, token, ','); c.y.push_back(std::stod(token));
    std::getline(ss, token, ','); ignored = std::stod(token); //z
    std::getline(ss, token, ','); c.isOnCurve.push_back(parseBool(trim(token)));
    std::getline(ss, token, ','); ignored = std::stod(token); //angle
    std::getline(ss, token, ','); c.isMdlVtx.push_back(parseBool(trim(token)));
  }

  return c;
}
