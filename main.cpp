//#include <iostream>

//#include <boost/numeric/odeint.hpp>

//using namespace boost::numeric::odeint;

////[ rhs_function
///* The type of container used to hold the state vector */
//typedef std::vector< double > state_type;

//const double gam = 0.15;

///* The rhs of x' = f(x) */
//void harmonic_oscillator( const state_type &x , state_type &dxdt , const double /* t */ )
//{
//    dxdt[0] = x[1];
//    dxdt[1] = -x[0] - gam*x[1];
//}
////]



////[ rhs_class
///* The rhs of x' = f(x) defined as a class */
//class harm_osc {

//    double m_gam;

//public:
//    harm_osc( double gam ) : m_gam(gam) { }

//    void operator() ( const state_type &x , state_type &dxdt , const double /* t */ )
//    {
//        dxdt[0] = x[1];
//        dxdt[1] = -x[0] - m_gam*x[1];
//    }
//};
////]


//int main()
// {
//    //[ state_initialization
//    state_type x(4);
//    x[0] = 1.0; // start at x=1.0, p=0.0
//    x[1] = 0.0;
//    //]


////    //[ integration
////    size_t steps = integrate( harmonic_oscillator ,
////            x , 0.0 , 10.0 , 0.1 );
////    //]


////    //[ define_const_stepper
////    runge_kutta4< state_type > stepper;
////    integrate_const( stepper , harmonic_oscillator , x , 0.0 , 10.0 , 0.01 );
////    //]





//    //[ integration_class
//    harm_osc ho(0.15);
//    size_t steps = integrate( ho ,
//            x , 0.0, 10.0, 1.0 );
//    //]

//    return 0;
//}

//////////////////////////////////////////////////////////////////////////////////////////////

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>

//double rk4(double(*f)(double, double), double dx, double x, double y)
//{
//    double	k1 = dx * f(x, y),
//        k2 = dx * f(x + dx / 2, y + k1 / 2),
//        k3 = dx * f(x + dx / 2, y + k2 / 2),
//        k4 = dx * f(x + dx, y + k3);
//    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
//}

//double rate(double x, double y)
//{
//    return x * sqrt(y);
//}

//int main(void)
//{
//    double *y, x, y2;
//    double x0 = 0, x1 = 10, dx = .1;
//    int i, n = 1 + (x1 - x0)/dx;
//    y = (double *)malloc(sizeof(double) * n);

//    for (y[0] = 1, i = 1; i < n; i++)
//        y[i] = rk4(rate, dx, x0 + dx * (i - 1), y[i-1]);

//    printf("x\ty\trel. err.\n------------\n");
//    for (i = 0; i < n; i += 10) {
//        x = x0 + dx * i;
//        y2 = pow(x * x / 4 + 1, 2);
//        printf("%g\t%g\t%g\n", x, y[i], y[i]/y2 - 1);
//    }

//    return 0;
//}

///////////////////////////////////////////////////////////////////////////////////

//// Simple RK4 integration framework
//// Copyright (c) 2004, Glenn Fiedler
//// http://www.gaffer.org/articles

//#include <stdio.h>
//#include <math.h>
//#include <iostream>
//#include <fstream>
//#include <cstdlib>


//struct State
//{
//    float x;
//    float v;
//};

//struct Derivative
//{
//    float dx;
//    float dv;
//};

//float acceleration(const State &state, float t)
//{
//    const float k = 10;
//    const float b = 1;
//    return - k*state.x - b*state.v;
//}

//Derivative evaluate(const State &initial, float t)
//{
//    Derivative output;
//    output.dx = initial.v;
//    output.dv = acceleration(initial, t);
//    return output;
//}

//Derivative evaluate(const State &initial, float t, float dt, const Derivative &d)
//{
//    State state;
//    state.x = initial.x + d.dx*dt;
//    state.v = initial.v + d.dv*dt;
//    Derivative output;
//    output.dx = state.v;
//    output.dv = acceleration(state, t+dt);
//    return output;
//}

//void integrate(State &state, float t, float dt)
//{
//    Derivative a = evaluate(state, t);
//    Derivative b = evaluate(state, t, dt*0.5f, a);
//    Derivative c = evaluate(state, t, dt*0.5f, b);
//    Derivative d = evaluate(state, t, dt, c);

//    const float dxdt = 1.0f/6.0f * (a.dx + 2.0f*(b.dx + c.dx) + d.dx);
//    const float dvdt = 1.0f/6.0f * (a.dv + 2.0f*(b.dv + c.dv) + d.dv);

//    state.x = state.x + dxdt*dt;
//    state.v = state.v + dvdt*dt;
//}

//int main()
//{
//    int countersteps=100;

//    //drawing part
//    //
//    std::ofstream myfile;
//    myfile.open ("testIntegration_Simple RK4 integration framework.ppm");
//    myfile << "P3\n";
//    myfile << "1000 1000\n";
//    myfile << "255\n";//max picel value number
////    for(int i=0;i<1000;i++)
////    {
////        for(int j=0;j<1000;j++)
////        {
////            myfile <<"255"<< " 255 "<<"255 ";
////        }
////        myfile <<"\n";
////    }



////    myfile.close();

//    //

//    State state;
//    state.x = 100;
//    state.v = 0;

//    float t = 0;
//    float dt = 0.1;


//////    while (fabs(state.x)>0.001f || fabs(state.v)>0.001f)
//    while(countersteps-->0)
//    {

//        printf("%.2f, %.2f\n", state.x, state.v);

//        for(int i=0;i<100;i++)
//        {
//            myfile <<  abs(state.x)*i<<" "<< abs(state.v)*i/20 << " "<< 0 << " ";

//            for(int j=0;j<100;j++)
//            {
//                myfile << "\n";
//                myfile <<  abs(state.x)<<" "<< abs(state.v) << " "<< 0 << " ";
//            }
//        }



//        integrate(state, t, dt);
//        t += dt;
//    }

//    myfile << "\n";
//    myfile.close();


////    getc(stdin);

//    return 0;
//}



//State Definition and Integration of physics' based simulation START
///////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

typedef boost::array< double , 3 > state_type;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}

void write_lorenz( const state_type &x , const double t )
{
    cout << /*t <<*/ '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
}
///////////////////////////////////////////////////////////////////////////////////
//State Definition and Integration of physics' based simulation END




//A 1st order 1D DE solver. Integration START
///////////////////////////////////////////////////////////////////////////////////

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "OneDSolver.h"

#include <fstream>
//#include <cstdlib>
//#include <stdio.h>
//#include <math.h>
//#include <iostream>


#define EULER    0
#define MODIFIED_EULER    1
#define HEUNS    2
#define MIDPOINT 3
#define RUNGEKUTTA_4 4
#define RUNGEKUTTA_6 5


double EvalFcn(double x)
{
   return(-0.05 * x);
}


//A 1st order 1D DE solver. Integration END
///////////////////////////////////////////////////////////////////////////////////
/// \brief main
/// \param argc
/// \param argv
/// \return
///
///
int main(int argc, char **argv)
{


//    State Definition and Integration of physics' based simulation
//    state_type x = { 10.0 , 1.0 , 1.0 }; // initial conditions
//    integrate( lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz );





//     A 1st order 1D DE solver. Integration
//     The darkest the color of the image->the more accurate the highest the order of integration -> the more accurate the results

        std::ofstream myfile;
        myfile.open (" A_1st_order_ 1D_DE_solver_Integration.ppm");
        myfile << "P3\n";
        myfile << "100 60\n";
        myfile << "255\n";//max picel value number


       double t;
       double dt=0.1;    /* Step size            */
       double T=10;     /* Simulation duration  */
       double y = 1;     /* Initial value        */


       OneDSolver s;
       for (t=0.1;t<T;t+=dt)
       {
          printf("%g %g\n",t,y);
          y = s.Solver1D(dt,y,EULER,(double (*)(double))EvalFcn);


          myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";

            for (int i=0;i<10;i++)
            {
                myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";
                myfile << "\n";
            }

          counter++;

       }


       for (t=0.1;t<T;t+=dt) {
          printf("%g %g\n",t,y);
          y = s.Solver1D(dt,y,MODIFIED_EULER,(double (*)(double))EvalFcn);


          myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";

          for (int i=0;i<10;i++)
          {
              myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";
              myfile << "\n";
          }

          counter++;
       }


       for (t=0.1;t<T;t+=dt) {
          printf("%g %g\n",t,y);
          y = s.Solver1D(dt,y,HEUNS,(double (*)(double))EvalFcn);


          myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";

          for (int i=0;i<10;i++)
          {
              myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";
              myfile << "\n";
          }
       }


       for (t=0.1;t<T;t+=dt) {
          printf("%g %g\n",t,y);
          y = s.Solver1D(dt,y,MIDPOINT,(double (*)(double))EvalFcn);


          myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";

          for (int i=0;i<10;i++)
          {
              myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";
              myfile << "\n";
          }
       }

       for (t=0.1;t<T;t+=dt) {
          printf("%g %g\n",t,y);
          y = s.Solver1D(dt,y,RUNGEKUTTA_4,(double (*)(double))EvalFcn);


          myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";

          for (int i=0;i<10;i++)
          {
              myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";
              myfile << "\n";
          }
       }

       for (t=0.1;t<T;t+=dt) {
          printf("%g %g\n",t,y);
          y = s.Solver1D(dt,y,RUNGEKUTTA_6,(double (*)(double))EvalFcn);


          myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";

          for (int i=0;i<10;i++)
          {
              myfile <<  int(abs(y)*200)<<" "<< int(abs(y)*200) << " "<< 0 << " ";
              myfile << "\n";
          }
       }

       myfile << "\n";
       myfile.close();




}
