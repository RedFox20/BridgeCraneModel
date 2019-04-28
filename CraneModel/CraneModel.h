#pragma once

namespace crane3d
{
	/**
	 * Output state of the model
	 */
	struct ModelState
	{
		double Alfa = 0.0; // α pendulum measured alfa angle
		double Beta = 0.0; // β pendulum measured beta angle

		double LiftLine = 0.0; // R lift-line length
		double RailOffset = 0.0;  // Xw distance of the rail with the cart from the center of the construction frame
		double CartOffset = 0.0;  // Yw distance of the cart from the center of the rail

		// Payload 3D coordinates
		double PayloadX = 0.0;
		double PayloadY = 0.0;
		double PayloadZ = 0.0;
	};

	// Coordinate system of the Crane model
	// X: outermost movement of the rail, considered as forward
	// Y: left-right movement of the cart
	// Z: up-down movement of the payload
    class Model
    {
        static constexpr double Mpayload = 1.000; // Mc mass of the payload
        static constexpr double Mcart    = 1.155; // Mw mass of the cart
        static constexpr double Mrail    = 2.200; // Ms mass of the moving rail
        static constexpr double Mrailcart    = Mcart + Mrail;    // Mcart + Mrail    (Mw + Ms)
        static constexpr double Mcartpayload = Mcart + Mpayload; // Mcart + Mpayload (Mw + Mc)

		// Friction forces
		static constexpr double Tx = 100.0; // rail friction
		static constexpr double Ty = 82.0; // cart friction
		static constexpr double Tr = 75.0; // liftline friction 

        double Xw = 0.0; // distance of the rail with the cart from the center of the construction frame
        double Yw = 0.0; // distance of the cart from the center of the rail;
        double R = 0.0; // length of the lift-line
        double Alfa = 0.0; // angle between y axis (cart moving left-right) and the lift-line
        double Beta = 0.0; // angle between negative direction on the z axis and the projection
                           // of the lift-line onto the xz plane

		// velocity time derivatives
		double Xw_vel = 0.0;
		double Yw_vel = 0.0;
		double R_vel = 0.0;
		double Alfa_vel = 0.0;
		double Beta_vel = 0.0;

		//double S = 0.0; // reaction force in the lift-line acting on the cart

    public:

        Model();

        /**
         * @param deltaTime Time since last update
         * @param Frail force driving the rail with cart (Fx)
         * @param Fcart force driving the cart along the rail (Fy)
         * @param Fline force controlling the length of the lift-line (Fr)
         * @return New state of the crane model
         */
		ModelState Update(double deltaTime, double Frail, double Fcart, double Fline);

		/**
		 * @return Current state of the crane:
		 *  distance of the rail, cart, length of lift-line and swing angles of the payload
		 */
		ModelState GetState() const;

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
