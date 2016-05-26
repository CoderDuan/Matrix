# Matrix
Implemented a template class named Mat in C++. <br>
Support calculating the determinent value, transpose matrix, inverse matrix(if exists) of the matrix and some other algorithm.
To use this Mat class, just include "mat.h"
```cpp
#include "mat.h"
```
For initialize, you can use an array pointer to initialize. 
If you use default constructor, all elements in the matrix will be initialized to 0.
```cpp
Mat<4,3> mat; // A mat with 4 rows and 3 column. By default, the type of each element is float.
              // Elements of mat is initialized to 0
```
or
```cpp
Mat<4,3> mat = Mat<4,3>(); // elements of mat is initialized to 0
```
or
```cpp
float data[4] = {1,2,3,4};
Mat<2,2> mat(data);
```
use printData() method to print the matrix
```cpp
float data[4] = {1,2,3,4};
Mat<2,2> mat(data);
mat.printData();
```
output will be
>1    2<br>
>3    4<br>

other methods:
```cpp
Mat<2,4,float> mat(data);
Mat<4,2> trans = mat.transpose();
```
```cpp
Mat<4,4,float> mat(data);
Mat<4,4> invr = mat.inverse();
```
```cpp
Mat<4,4,float> mat(data);
float deter = mat.determinant();
```








