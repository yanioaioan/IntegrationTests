#ifndef ONEDSOLVER_H
#define ONEDSOLVER_H


class OneDSolver
{
public:
    OneDSolver();
    ~OneDSolver();

    double Solver1D(double h,double y0,int method,double (*fcn)(double));

};

#endif // ONEDSOLVER_H
