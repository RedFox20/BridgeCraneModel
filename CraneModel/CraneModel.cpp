#include "CraneModel.h"
#include <cmath>

namespace crane3d
{
    static constexpr double g = 9.81; // gravity constant, 9.81m/s^2

    //////////////////////////////////////////////////////////////////////

	// @return Velocity
	inline double vel(double x1, double x2, double deltaTime)
	{
		return (x2 - x1) / deltaTime; // v = dx / dt
	}

	//////////////////////////////////////////////////////////////////////

    Model::Model()
    {
    }

	ModelState Model::Update(double deltaTime, double Frail, double Fcart, double Fline)
    {
        // time derivative - velocity
        // second derivative - acceleration

		// Complete Non-Linear model
		double N1 = CartNetAccel(Fcart);
		double N2 = RailNetAccel(Frail);
		double N3 = LineNetAccel(Fline);

		double μ1 = PayloadCartRatio;
		double μ2 = PayloadRailCartRatio;

		// x1..x10 as per 3DCrane mathematical model description
		double x1 = Yw;
		double x2 = Yw_vel; // velocity derivative of Yw
		double x3 = Xw;
		double x4 = Xw_vel;
		double x5 = Alfa;
		double x6 = Alfa_vel;
    	double x7 = Beta;
		double x8 = Beta_vel;
		double x9 = R;
		double x10 = R_vel;

		// NOTE: variable names left exactly the same as the original math. model descr. 
		double s5 = sin(x5);
		double s7 = sin(x7);
		double c5 = cos(x5);
		double c7 = cos(x7);

		double V5 = c5 * s5*x8*x8*x9 - 2 * x10*x6 + g*c5*c7;
		double V6 = 2 * x8*(c5*x6*x9 + s5 * x10) + g*s7;
		double V7 = s5 * s5*x8*x8*x9 + g * s5*c7 + x6*x6*x9;

		// derived x1..x10 as per 3DCrane mathematical model description
		double d1 = x2;
		double d2 = N1 + μ1 * c5 * N3;
		double d3 = x4;
    	double d4 = N2 + μ2 * s5 * s7 * N3;
		double d5 = x6;
		double d6 = (s5 * N1 - c5 * s7 * N2 + (μ1 - μ2 * s7*s7) * c5*s5*N3 * V5) / x9;
		double d7 = x8;
		double d8 = -(c7*N2 + μ2 * s5*c7*s7*N3 + V6) / (s5*x9);
		double d9 = x10;
		double d10 = -c5 * N1 - s5 * s7*N2 - (1 + μ1 * c5*c5 + μ2 * s5*s5*s7*s7)*N3 + V7;

		// map the derived state to the current internal state:
		Yw   = d1;
		Xw   = d3;
		Alfa = d5;
		Beta = d7;
		R    = d9;
		Yw_vel   = d2;
		Xw_vel   = d4;
		Alfa_vel = d6;
		Beta_vel = d8;
		R_vel    = d10;

		// and finally, return the observable current state
		return GetState();
    }

	ModelState Model::GetState() const
	{
		ModelState s;
		s.Alfa = Alfa;
		s.Beta = Beta;
		s.RailOffset = Xw;
		s.CartOffset = Yw;
		s.LiftLine = R;
		return s;
	}

    //////////////////////////////////////////////////////////////////////
}
