#include "CraneModel/CraneModel.h"
#include <cstdio>
#include <cstdlib>

int main()
{
    using namespace crane3d;

    SimulationState initial;
    Model craneModel {initial};

    double Fx = 0.0; // force driving the rail
    double Fy = 1.0; // force along the rail
    double Fr = 0.0; // force driving the cable
    SimulationState derived = craneModel.GetDerivative(Fx, Fy, Fr);

    system("pause");
}
