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

    class ModelOld
    {
        SimulationState S;

    public:

        explicit ModelOld(const SimulationState& initialState);

        SimulationState Get() const;

        /**
         * @param Fx force driving the rail with cart
         * @param Fy force driving the cart along the rail
         * @param Fr force controlling the length of the lift-line
         */
        SimulationState GetDerivative(double Fx, double Fy, double Fr) const;
    };

    struct Vec3d { double X, Y, Z; };

    /**
     * NEW model state passed to update
     */
    struct ModelParameters
    {
        double Alfa = 0.0; // α pendulum measured alfa angle
        double Beta = 0.0; // β pendulum measured beta angle

        double LiftLine = 0.0; // R new lift-line length
        double OffsetX = 0.0;  // Xw distance of the rail with the cart from the center of the construction frame
        double OffsetY = 0.0;  // Yw distance of the cart from the center of the rail
    };

	struct ModelState
	{
		double x1 = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0, x5 = 0.0;
		double x6 = 0.0, x7 = 0.0, x8 = 0.0, x9 = 0.0, x10 = 0.0;

		ModelState() = default;
	};

    class Model
    {
        static constexpr double Mpayload = 1.000; // Mc mass of the payload
        static constexpr double Mcart    = 1.155; // Mw mass of the cart
        static constexpr double Mrail    = 2.200; // Ms mass of the moving rail
        static constexpr double Mrailcart    = Mcart + Mrail;    // Mcart + Mrail    (Mw + Ms)
        static constexpr double Mcartpayload = Mcart + Mpayload; // Mcart + Mpayload (Mw + Mc)

        // Coordinates of the payload
        // X: outermost movement of the rail, considered as forward
        // Y: left-right movement of the cart
        // Z: up-down movement of the payload
        Vec3d PayloadPos { 0.0, 0.0, 0.0 };

        double Xw = 0.0; // distance of the rail with the cart from the center of the construction frame
        double Yw = 0.0; // distance of the cart from the center of the rail;

        double R = 0.0; // length of the lift-line
        double S = 0.0; // reaction force in the lift-line acting on the cart
        double Alfa = 0.0; // angle between y axis (cart moving left-right) and the lift-line
        double Beta = 0.0; // angle between negative direction on the z axis and the projection
                           // of the lift-line onto the xz plane

        // Friction forces
        static constexpr double Tx = 100.0; // rail friction
        static constexpr double Ty = 82.0; // cart friction
        static constexpr double Tr = 75.0; // liftline friction 

    public:

        Model();

        /**
         * @param deltaTime Time since last update
         * @param s current state of the model
         * @param Frail force driving the rail with cart (Fx)
         * @param Fcart force driving the cart along the rail (Fy)
         * @param Fline force controlling the length of the lift-line (Fr)
         */
        void Update(double deltaTime, const ModelState& s, double Frail, double Fcart, double Fline);

    private:

		// μ1 = Mc / Mw
		static constexpr double PayloadCartRatio = Mpayload / Mcart;
		// μ2 = Mc / (Mw + Ms)
		static constexpr double PayloadRailCartRatio = Mpayload / Mrailcart;

		// u1 = Fy / Mw
		// u1 = Fcart / Mcart
		static double CartAccel(double Fcart) { return Fcart / Mcart; }
		// u2 = Fx / (Mw + Mc)
		// u2 = Frail / Mcartpayload
		static double RailAccel(double Frail) { return Frail / Mcartpayload; }
		// u3 = Fr / Mc
		// u3 = Fline / Mpayload 
		static double LineAccel(double Fline) { return Fline / Mpayload; }

		// T1 = Ty / Mw
		static constexpr double CartFrictionAccel = Ty / Mcart;
		// T2 = Tx / (Mw + Mc)
		static constexpr double RailFrictionAccel = Tx / Mcartpayload;
		// T3 = Tr / Mc
		static constexpr double LineFrictionAccel = Tr / Mpayload;

		// N1 = u1 - T1
		static double CartNetAccel(double Fcart) { return CartAccel(Fcart) - CartFrictionAccel; }
		// N2 = u2 - T2
		static double RailNetAccel(double Frail) { return RailAccel(Frail) - RailFrictionAccel; }
		// N3 = u3 - T3
		static double LineNetAccel(double Fline) { return LineAccel(Fline) - LineFrictionAccel; }

		// s = S / Mc
		static double LiftReactionAccel(double Sline) { return Sline / Mpayload; }

    };

}
