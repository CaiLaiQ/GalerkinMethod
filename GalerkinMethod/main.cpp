//
//  main.cpp
//  DifEquationTask4
//
//  Created by Равил Файрушин on 20.03.16.
//  Copyright © 2016 Равил Файрушин. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

void clearFile(char *fileName) {
    ofstream myfile;
    myfile.open (fileName, ios::out | ios::trunc);
    myfile.close();
}
void writeTextToFile(char *fileName, char *text) {
    ofstream myfile;
    myfile.open (fileName, ios::out | ios::app);
    myfile<<text<<' '<<endl;
    myfile.close();
}
void writeNumberToFile(char *fileName, double N) {
    ofstream myfile;
    myfile.open (fileName, ios::out | ios::app);
    myfile<<N<<' '<<endl;
    myfile.close();
}
void writeToFile(char *fileName, double *x, int length) {
    ofstream myfile;
    myfile.open (fileName, ios::out | ios::app);
    for (int i = 0; i < length; i++) {
        myfile<<x[i]<<' ';
    }
    myfile<<endl;
    myfile.close();
}

double sum(const function<double(double)> &f, double *nodes, int dim) {
    double sum = 0;
    for (int i = 0; i < dim; i++) {
        sum += f(nodes[i]);
    }
    return sum;
}

/*метод Гаусса-Зейделя*/
// Условие окончания
bool converge(double *xk, double *xkp, int n, double eps)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
    }
    if(sqrt(norm) >= eps)
        return false;
    return true;
}
void calculateSolution(double *a, double *x, double *b, int n, double eps) {
    double *p = new double[n];
    do
    {
        for (int i = 0; i < n; i++)
            p[i] = x[i];
        
        for (int i = 0; i < n; i++)
        {
            double var = 0;
            for (int j = 0; j < i; j++)
                var += (a[i*n+j] * x[j]);
            for (int j = i + 1; j < n; j++)
                var += (a[i*n+j] * p[j]);
            x[i] = (b[i] - var) / a[i*n+i];
        }
    }
    while (!converge(x, p, n, eps));
}
/**/

int main() {
    int i,n;
    
    /*Создаем узлы x[i] для интегрирования*/
    double alpha,beta;
    double h = 0.01;
    alpha = -1.0;
    beta = 1.0;
    int M = (beta-alpha)/h+1;
    double *x = new double[M];
    for (int i = 0; i < M; i++) {
        x[i] = (beta-alpha)/2.0*cos(M_PI*(2.0*i - 1)/(2.0 * M)) + (beta + alpha)/2.0;
    }
    /*конец*/
    
    /*Здесь будут объявлены функции*/
    function<double (double)> fx = [](double _x)->double {
        return -M_PI*_x + 1.0/sqrt(1.0 -_x*_x);
    };
    function<double (double,double)> Kxt = [] (double _x, double _t)-> double {
        return (_t+_x)*(_t+_x)*_t;
    };
    function<double (int, double)> Tn0x = [] (int n, double _x) -> double {
        return cos(n * acos(_x));
    };
    function<double (int, double)> Tnx = [&] (int n, double _x) -> double {
        double x1 = 2.0*_x/(beta-alpha) - (beta +alpha)/(beta-alpha);
        return Tn0x(n,x1);
    };
    
    function<double(double)> func1 = [Tnx,&i,&n] (double xi) -> double {
        return Tnx(i,xi)*Tnx(n,xi);
    };
    function<double(double)> func2 = [&] (double xi) -> double {
        double sum = 0;
        for (int l = 0; l < M; l++) {
            sum += Kxt(xi,x[l]) * Tnx(n,x[l]);
        }
        return sum;
    };
    function<double(double)> func3 = [&] (double xm) -> double {
        return Tnx(i,xm)*func2(xm)*sqrt((beta-xm)*(xm-alpha));
    };
    function<double(double)> func4 = [&] (double xm) -> double {
        return Tnx(i,xm);
    };
    function<double(double)> func5 = [&] (double xm) -> double {
        return Tnx(i,xm)*xm*sqrt((beta-xm)*(xm-alpha));
    };
    
    /*Конец объявлений функций*/
    
    /*Здесь мы создаем матрицу с[n][n], вектор-столбец b[n] и заполняем их*/
    int N = 10;
    double *c = new double[N*N];
    double *b = new double[N];
    double integral1, integral2;
    for (i = 0; i < N; i++) {
        for (n = 0; n < N; n++) {
            integral1 = M_PI/M * sum(func1, x, M);
            integral2 = -(M_PI/M)*(M_PI/M) * sum(func3, x, M);
            c[i*N+n] = integral1 + integral2;
            //            cout<<c[i*N + n]<<' ';
        }
        integral1 = M_PI/M * sum(func4, x, M);
        integral2 = -M_PI*M_PI/M * sum(func5, x, M);
        b[i] = integral1 + integral2;
        //        cout<<b[i]<<' ';
    }
    
    /* конец*/
    char *fileName = strdup("/Users/Ravil/Documents/projects/Differential equtions/DifEquationTask3/DifEquationTask3/nodes.txt");
    clearFile(fileName);
    writeNumberToFile(fileName, N);
    writeToFile(fileName, c, N * N);
    writeToFile(fileName, b, N);
    
    return 0;
}



















