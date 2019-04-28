#include "CraneModel.h"
#include <cmath>

namespace crane3d
{
    // Simulation Constants
    // X: Right
    // Y: Out of the screen
    // Z: Up
    static constexpr double mc = 1;
    static constexpr double mw = 1.155;
    static constexpr double ms = 2.20;
    static constexpr double g = 9.81; // gravity constant, 9.81m/s^2
    static constexpr double Tx = 100; // friction
    static constexpr double Ty = 82;
    static constexpr double Tz = 75;

    static constexpr double My = mw + ms;
    static constexpr double M = mw * mw + mc * ms + mw * ms + mc * mw;

    static constexpr double mi1 = mc / mw;
    static constexpr double mi2 = mc / (mw + ms);

    static constexpr double T1 = Ty / mw;
    static constexpr double T2 = Tx / (mw + mc);   //poprawiono z ms na mc - GRI_10:00
    static constexpr double T3 = Tz / mc;

    static constexpr double Tsx = 5 / (mw + ms);
    static constexpr double Tsy = 7.5 / mw;
    static constexpr double Tsz = 10 / mc;

    ModelOld::ModelOld(const SimulationState& initialState) : S{initialState}
    {
        S.ALFA = (3.1415926 / 2) - initialState.ALFA; // Rotate 180 deg
    }

    SimulationState ModelOld::Get() const
    {
        SimulationState state = S;
        state.ALFA = (3.1415926 / 2) - S.ALFA; // Rotate 180 deg
        return state;
    }

    static double sign(double x)
    {
        if (x == 0) { return 0; }
        if (x > 0) { return 1; }
        if (x < 0) { return -1; }
        return 0;
    }
    
    SimulationState ModelOld::GetDerivative(double Fx, double Fy, double Fr) const
    {
        const double RR = S.Rstale;

        // @todo Reconstruct this further by using 3DCrane user manual physics model description
        const double u1 = Fy / mw;
        const double u2 = Fx / (mw + mc);
        const double u3 = Fr / mc;

        const double ca = cos(S.ALFA);
        const double sa = sin(S.ALFA);
        const double cb = cos(S.VBETA);
        const double sb = sin(S.VBETA);

        const double t = S.SAMPLE_TIME_ARG;
        const double vx = S.VXWOZEK;
        const double vy = S.VYWOZEK;
        const double va = S.VALFA;
        const double vr = S.VRLINKI;
        const double r = S.RLINKI;
        const double sit = sign(t);
        const double svx = sign(vx);
        const double svy = sign(vy);

        const double sasb = sa*sb;
        const double casb = ca*sb;
        const double cami1 = ca*mi1;
        const double cami2 = ca*mi2;

        const double t3_Tsz = (-t * T3 - Tsz * sit);
        const double t2_Tsx = (vy * T2 + Tsx * svy);
        const double mi2sasb = mi2*sasb;

        const double un2_t3 = u3 - t3_Tsz;

        const double cami1_un2_t3 = cami1*un2_t3;
        const double mi2sasb_un2_t3 = mi2sasb*un2_t3;

        SimulationState d;
        d.XWOZEK = vx;
        d.YWOZEK = vy;
        d.ALFA = va;
        d.BETA = r;

        d.VXWOZEK = u1 - vx * T1 - Tsy * svx + cami1_un2_t3;
        d.VYWOZEK = u2 - vy * T2 - Tsx * svy + mi2sasb_un2_t3;
        d.VALFA = (-(vx * T1 + Tsy * svx)*sa
                + u1 * sa
                + ca*vr * sa*r * r
                - ca*u2 * sb
                + ca*cb*g
                - casb*mi2sasb_un2_t3
                + casb*t2_Tsx
                + sa*cami1_un2_t3
                - 2 * t * va) / vr;
        d.VBETA = -(cb*u2
                + g * sb
                + 2 * vr * ca*va * r
                + cb*mi2sasb_un2_t3
                + 2 * t * sa*r 
                - cb*t2_Tsx
                ) / (vr * sa);
        
        if (RR == 1)
        {
            d.VRLINKI = 0;
            d.SAMPLE_TIME_ARG = 0;
        }
        else
        {
            d.RLINKI = t;
            d.VRLINKI = ca*(vx * T1 + Tsy * svx)
                    - ca*u1 + vr * sa*sa*r * r
                    - u2 * sasb + sa*cb*g
                    - mi2sasb*sa*sb*u3
                    + mi2sasb*sa*sb*t3_Tsz
                    + t2_Tsx*sasb
                    - mi1 * u3
                    + mi1 * u3 * sa*sa
                    + mi1 * t3_Tsz
                    - mi1 * t3_Tsz*sa*sa
                    + vr * va * va
                    - u3
                    - t * T3
                    - Tsz * sit;
        }

        return d;
    }

    //////////////////////////////////////////////////////////////////////

    Model::Model()
    {
		
    }

    void Model::Update(double deltaTime, const ModelState& s, double Frail, double Fcart, double Fline)
    {
        double sinA = sin(Alfa), cosA = cos(Alfa);
        double sinB = sin(Beta), cosB = cos(Beta);

        S = Fline - Tr;
		// second derivates
		double lineNetAccel = LineNetAccel(Fline);
		double aX = -lineNetAccel * sinA * sinB;
		double aY = -lineNetAccel * cosA;
		double aZ = lineNetAccel * sinA * cosB - g;

		double aXw = RailNetAccel(Frail) + RailAccel(Frail)*lineNetAccel * sinA * sinB;
		double aYw = CartNetAccel(Fcart) + CartAccel(Fcart)*lineNetAccel * cosA;

        // time derivative - velocity
        // second derivative - acceleration

		double N1 = CartNetAccel(Fcart);
		double N2 = RailNetAccel(Frail);
		double N3 = LineNetAccel(Fline);

		double μ1 = PayloadCartRatio;
		double μ2 = PayloadRailCartRatio;

		double s5 = sin(s.x5);
		double s7 = sin(s.x7);

		double c5 = cos(s.x5);
		double c7 = cos(s.x7);

		double V5 = c5 * s5*s.x8*s.x8*s.x9 - 2 * s.x10*s.x6 + g*c5*c7;
		double V6 = 2 * s.x8*(c5*s.x6*s.x9 + s5 * s.x10) + g*s7;
		double V7 = s5 * s5*s.x8*s.x8*s.x9 + g * s5*c7 + s.x6*s.x6*s.x9;

		ModelState d;
		d.x1 = s.x2;
		d.x2 = N1 + μ1 * c5 * N3;
		d.x3 = s.x4;
    	d.x4 = N2 + μ2 * s5 * s7 * N3;
		d.x5 = s.x6;
		d.x6 = (s5 * N1 - c5 * s7 * N2 + (μ1 - μ2 * s7*s7) * c5*s5*N3 * V5) / s.x9;
		d.x7 = s.x8;
		d.x8 = -(c7*N2 + μ2 * s5*c7*s7*N3 + V6) / (s5*s.x9);
		d.x9 = s.x10;
		d.x10 = -c5 * N1 - s5 * s7*N2 - (1 + μ1 * c5*c5 + μ2 * s5*s5*s7*s7)*N3 + V7;


    }

    //////////////////////////////////////////////////////////////////////
}
