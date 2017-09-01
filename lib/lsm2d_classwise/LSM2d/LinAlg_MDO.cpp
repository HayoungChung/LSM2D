#ifndef __LinAlg_MDO_CPP__
#define __LinAlg_MDO_CPP__
// a function that transposes a matrix

#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>

using namespace std;

typedef vector<vector<double> > Matrix;
typedef vector<double> Vector;

// Declaration
Matrix transpose(Matrix);
double dot(Vector&, Vector&);
Vector dot(Matrix&, Vector&);
Matrix dot(Matrix&, Matrix&);
void print(Matrix&);
void print(Vector&);

template <typename _Scalar>
void print(vector<vector<vector<_Scalar> > >&);

Matrix transpose(Matrix A){
    int n = A.size(), m = A[0].size();

    Matrix output;
    output.resize(m, vector<double>(n,0.0));

    for (int ii = 0; ii < n; ii++){
        for (int jj = 0; jj < m; jj++){
            output[jj][ii] = A[ii][jj];
        }
    }
    return output;
}

double dot(Vector& input1, Vector& input2){
    int n1 = input1.size(), n2 = input2.size();
    assert(n1==n2);
    double output = 0;
    for (int ii = 0; ii < n1; ii++){
        output += input1[ii]*input2[ii];
    }    
    return output;
}

Vector dot(Matrix& input1, Vector& input2){
    int n = input1.size(), m = input1[0].size();
    int n2 = input2.size();
    assert(m==n2);
    vector <double> output(n,0);

    for (int ii = 0; ii < n; ii++){
        for (int jj =0; jj < m; jj++){
            output[ii] += input1[ii][jj]*input2[jj];
        }
    }    
    return output;
}

Matrix dot(Matrix& input1, Matrix& input2){
    int n1 = input1.size(), m1 = input1[0].size();
    int n2 = input2.size(), m2 = input2[0].size();
    assert(m1==n2) ;
    
    Matrix output;
    output.resize(n1, vector<double>(m2, 0.0));

    for (int ii = 0; ii < n1; ii++){
        for (int jj = 0; jj < m2; jj++){
            for (int kk = 0; kk < m1; kk++){
                output[ii][jj] += input1[ii][kk]*input2[kk][jj];
            }
        }
    }    
    return output;
}

void print(Matrix& A){
    int n = A.size(), m = A[0].size();
    cout << "size: (" << n << ", "<< m << ")\n";
    for (int ii = 0; ii < n; ii++){
        for (int jj = 0; jj < m; jj++){
            cout << A[ii][jj] << " ";
        }
        cout << "\n";
    }
}

void print(Vector& V){
    int n = V.size();
    cout << "size: (" << n << ")\n[";
    for (int ii = 0; ii < n; ii++){
            cout << V[ii] << " ";
    }
    cout << "]\n";
}

template <typename _Scalar>
void print(vector<vector<vector<_Scalar> > >& T){
    int n = T.size(), m = T[0].size(), l = T[0][0].size();
    for (int ii = 0; ii < n; ii ++){
        for (int jj = 0; jj < m; jj ++){
            for (int kk = 0; kk < l; kk ++){
                cout << "(" << ii << ","<< jj << "," << kk << "): " << T[ii][jj][kk] << endl;
            }
        }
    } 
}

Matrix operator*(Matrix& input, double a){
    int n = input.size(), m = input[0].size();
    Matrix output;
    output.resize(n, Vector(m,0.0));
    
    for (int i = 0; i < n; i ++){
        for (int j = 0; j < m; j ++){
            output[i][j] = input[i][j] * a;        
        }
    }
    return output;
}

Vector operator*(Vector& input, double a){
    int n = input.size();
    Vector output;
    output.resize(n, 0);
    for (int i = 0; i < n; i ++){
        output[i] = input[i] * a;
    }
    return output;
}

Matrix operator+(Matrix& input, double a){
    int n = input.size(), m = input[0].size();
    Matrix output;
    output.resize(n, Vector(m,0.0));

    for (int i = 0; i < n; i ++){
        for (int j = 0; j < m; j ++){
            output[i][j] = input[i][j] + a;
        }
    }
    return output;
}

Matrix operator+(Matrix& input, Matrix& input2){
    int n = input.size(), m = input[0].size();
    Matrix output;
    output.resize(n, Vector(m,0.0));
    for (int i = 0; i < n; i ++){
        for (int j = 0; j < m; j ++){
            output[i][j] = input[i][j] + input2[i][j];
        }
    }
    return output;
}

Vector operator+(Vector& input, double a){
    int n = input.size();
    Vector output;
    output.resize(n, 0);

    for (int i = 0; i < n; i ++){
        output[i] = input[i] + a;
    }
    return output;
}

Vector operator+(Vector& input, Vector& input2){
    int n = input.size();
    Vector output;
    output.resize(n, 0.0);
    for (int i = 0; i < n; i ++){
            output[i] = input[i] + input2[i];
    }
    return output;
}


#if __DEBUGFLAG__
int main(){
    Matrix A;
    A = {{1,2,3},{4,5,6},{7,8,9}};
    Matrix B = transpose(A);

    vector<double> a1 = {1, 2, 3};
    vector<double> b1 = {1, 2, 3};
    double c1 = dot(a1,b1);
    cout << c1 << endl; // verified

    vector<vector<double> > a2 = {{1, 2, 3}, {2,3,4}, {4,5,6}};
    vector<double> b2 = {1, 2, 3};
    vector<double> c2 = dot(a2,b2);
    // for (int ii = 0; ii < 3; ii++) cout << c2[ii] << endl; // verified
    print(c2);

    vector<vector<double> > a3 = {{1, 2, 3}, {2,3,4}, {4,5,6}};
    vector<vector<double> > b3 = {{1, 2, 3}, {2,3,4}, {4,5,6}};
    vector<vector<double> > c3 = dot(a3,b3);
    // for (int ii = 0; ii < 3; ii++) {
    //     cout << "\n";
    //     for (int jj = 0; jj < 3; jj++) {
    //         cout << c3[ii][jj] << " "; 
    //     }
    //     cout << "\n";
    // }
    print(c3);
}
#endif

#endif