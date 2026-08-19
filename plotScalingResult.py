import csv
import matplotlib.pyplot as plt

#Parse the scaling result csv file
x = []
y1 = []
y2 = []
modeCount = 0;
with open('scalingTestResult.csv', mode = 'r') as infile:
  csvFile = csv.DictReader(infile, fieldnames = ['mode', 'numSplines', 'paraCoords', '1stDerivDuration', '2ndDerivDuration'])
  for row in csvFile:
    x.append(int(row['numSplines']))
    y1.append(float(row['paraCoords'])/float(row['1stDerivDuration']))
    y2.append(float(row['paraCoords'])/float(row['2ndDerivDuration']))
    if (row['mode'] == 'uniform'):
      modeCount += 1

#Plotting 1st deriv result uniform
plt.plot(x[:modeCount], y1[:modeCount])
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("1st Derivative Result with Uniform Spline Length")
plt.show()

#Plotting 2nd deriv result uniform
plt.plot(x[:modeCount], y2[:modeCount])
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("2nd Derivative Result with Uniform Spline Length")
plt.show()

#Plotting 1st deriv result gaussian
plt.plot(x[:modeCount], y1[modeCount, len(y1)])
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("1st Derivative Result with Variable Spline Length Sampled from a Gaussian Distribution")
plt.show()

#Plotting 2nd deriv
plt.plot(x[:modeCount], y2[modeCount, len(y2)])
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("2nd Derivative Result with Variable Spline Length Sampled from a Gaussian Distribution")
plt.show()


