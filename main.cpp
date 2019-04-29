#include "CraneModel/Model.h"
#include <cstdio>
#include <cstdlib>

int main()
{
    using namespace crane3d;

    double Frail = 0.0; // force driving the rail
    double Fcart = 0.0; // force along the rail
    double Fline = 0.0; // force driving the cable

	Model model;
	ModelState state = model.Update(0.016, Frail, Fcart, Fline);
	state.Print();

    system("pause");
}
