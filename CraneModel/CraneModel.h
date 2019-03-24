#pragma once

namespace crane3d
{
    // @todo Reconstruct this using 3DCrane user manual physics model description
    struct SimulationState
    {
        double XWOZEK = 0.0;
        double VXWOZEK = 0.0;
        double YWOZEK = 0.0;
        double VYWOZEK = 0.0;
        double ALFA = 0.0;
        double VALFA = 0.0;
        double BETA = 0.0;
        double VBETA = 0.0;
        double RLINKI = 0.0;
        double VRLINKI = 0.0;
        double SAMPLE_TIME_ARG = 0.0;
        double Rstale = 0.0;
    };

    class Model
    {
        SimulationState S;

    public:

        explicit Model(const SimulationState& initialState);

        SimulationState Get() const;

        /**
         * @param Fx force driving the rail with cart
         * @param Fy force driving the cart along the rail
         * @param Fr force controlling the length of the lift-line
         */
        SimulationState GetDerivative(double Fx, double Fy, double Fr) const;
    };



}
