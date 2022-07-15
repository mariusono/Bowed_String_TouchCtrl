/*
  ==============================================================================

    ViolinString.cpp
    Created: 31 Mar 2022 8:19:26pm
    Author:  Marius Onofrei

  ==============================================================================
*/

#include <JuceHeader.h>
#include "ViolinString.h"

//==============================================================================
ViolinString::ViolinString (double freq, double fs) : fs (fs), freq (freq)
{
    // In your constructor, you should add any child components, and
    // initialise any special settings that your component needs.
    
    int ana = 3;
    
//    setOpaque(true);
    
    k = 1.0/fs; // time step
    
    L = 0.5;
    rho = 7850.0;
    r = 5e-4;
    T = 1000.0;
    E = 2e11;
    sig0 = 1;
    //sig1 = 0.01;

	sig1 = 0.015; // largest possible value in slider !! 
    
    A = double_Pi*r*r; // [m^2]
    I = double_Pi*r*r*r*r/4;
    K = sqrt(E*I/(rho*A));
    

    // Bowing params
    
    FB = 0.2; // [N]
    vB = 0.2; // bow velocity (m/s)
    a = 80.0; // friction law free parameter (1/m^2) (a) % decreasing sig increases the stick time
    tol = 1e-7; // tolerance for Newton-Raphson method
    A_NR = sqrt(2*a)*exp(0.5);
    
    inp_bow_pos_x = 0.5;
    
    stickFact = 1;
    
	//fN_var = 0.0;
	//vB_var = 0.0;
    
    fN_var = FB;
    vB_var = vB;
    
    // Derived params for num sim
    c = freq*2*L; // wave speed term
    h = sqrt((c * c * k * k + 4 * sig1 * k + sqrt(pow((c * c * k * k + 4 * sig1 * k),2) + 16 * K * K * k * k)) / 2.0);
    N = floor(L / h); // N is number of intervals !!! So N+1 points
    h = L / N;
    
    // Bow interpolation and spreading function init:
    I_B.resize(N+1, 0);
    J_B.resize(N+1, 0);

    // INITIALIZE STRING STATE VECTORS
    uVecs.reserve(3);

    for (int i = 0; i < 3; ++i)
        uVecs.push_back(std::vector<double>(N+1, 0));

    u.resize(3, nullptr);

    for (int i = 0; i < u.size(); ++i)
        u[i] = &uVecs[i][0];
    
}

ViolinString::~ViolinString()
{
}


void ViolinString::paint (juce::Graphics& g)
{
    // clear the background
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
    
    // choose your favourite colour
    g.setColour(Colours::cyan);
    
    // draw the state
    g.strokePath(visualiseState (g, 50000), PathStrokeType(2.0f));

}

void ViolinString::paint_crossSec (juce::Graphics& g)
{
    // clear the background
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
    
    // choose your favourite colour
    g.setColour(Colours::cyan);
    
    // draw the state
    if (paintFctIdx == 0)
    {
        g.strokePath(visualiseState (g, 50000), PathStrokeType(2.0f));
    }
    else if (paintFctIdx == 0)
    {
        g.strokePath(visualiseState_crossSec(g, 50000), PathStrokeType(2.0f));
    }
}


Path ViolinString::visualiseState_crossSec (Graphics& g, double visualScaling)
{
    // String-boundaries are in the vertical middle of the component
    double stringBoundaries = getWidth() / 2.0;


//    // String-boundaries are in the vertical middle of the component
//    double stringBoundaries = getHeight() / 2.0;

    // initialise path
    Path stringPath;

    // start path
    stringPath.startNewSubPath (0, -u[1][0] * visualScaling + stringBoundaries);

    double spacing = getHeight() / static_cast<double>(N);
    double y = spacing;

    for (int l = 1; l <= N; l++) // if you don't save the boundaries use l < N
    {
        // Needs to be -u, because a positive u would visually go down
        float newX = -u[1][l] * visualScaling + stringBoundaries;

        // if we get NAN values, make sure that we don't get an exception
        if (isnan(newX))
            newX = 0;

        stringPath.lineTo (y, newX);
        y += spacing;
    }
    // if you don't save the boundaries, and add a stringPath.lineTo (x, getWidth()) here to end the statedrawing

    return stringPath;
    
}


Path ViolinString::visualiseState (Graphics& g, double visualScaling)
{
    
    double stringBoundaries_y  = getWidth()/2;

    // initialise path
    Path stringPath;

    // start path
    //stringPath.startNewSubPath (stringBoundaries_y,u[1][0] * visualScaling);
    stringPath.startNewSubPath (stringBoundaries_y,0);

    double spacing = getHeight() / static_cast<double>(N);
//    double y = spacing;
    double y = 0;

    for (int l = 0; l <= N; l++) // go through all the points ! (N+1 points, N intervals)
    {
        // Needs to be -u, because a positive u would visually go down
        float newX = u[1][l] * visualScaling + getWidth()/2;

        // if we get NAN values, make sure that we don't get an exception
        if (isnan(newX))
            newX = 0;

        stringPath.lineTo (newX,y);
        y += spacing;
    }

    return stringPath;
    
    
    
//
//
//    // String-boundaries are in the vertical middle of the component
//    double stringBoundaries = getWidth() / 2.0;
//
//
////    // String-boundaries are in the vertical middle of the component
////    double stringBoundaries = getHeight() / 2.0;
//
//    // initialise path
//    Path stringPath;
//
//    // start path
//    stringPath.startNewSubPath (0, -u[1][0] * visualScaling + stringBoundaries);
//
//    double spacing = getHeight() / static_cast<double>(N);
//    double y = spacing;
//
//    for (int l = 1; l <= N; l++) // if you don't save the boundaries use l < N
//    {
//        // Needs to be -u, because a positive u would visually go down
//        float newX = -u[1][l] * visualScaling + stringBoundaries;
//
//        // if we get NAN values, make sure that we don't get an exception
//        if (isnan(newX))
//            newX = 0;
//
//        stringPath.lineTo (y, newX);
//        y += spacing;
//    }
//    // if you don't save the boundaries, and add a stringPath.lineTo (x, getWidth()) here to end the statedrawing
//
//    return stringPath;
}

void ViolinString::resized()
{

}


double ViolinString::process()
{

    double I_B_J_B = 0;
    double I_B_u = 0;
    double I_B_uPrev = 0;
    double I_B_dxx_u = 0;
    double I_B_dxx_uPrev = 0;
    double I_B_dxxxx_u = 0;
    int idx_p1;
    int idx_p2;
    int idx_m1;
    int idx_m2;

	/*
	if (counter == 0)
	{
		for (int l = 0; l <= N; ++l)
		{
			juce::Logger::getCurrentLogger()->outputDebugString("I_B: " + String(I_B[l]));
			juce::Logger::getCurrentLogger()->outputDebugString("J_B: " + String(J_B[l]));
		}
		init = false;
	}
	counter++;
	*/

	/*
	if (counter == 0)
	{
		for (int l = 0; l <= N; ++l)
		{
			juce::Logger::getCurrentLogger()->outputDebugString("uNext: " + String(u[0][l]));
			juce::Logger::getCurrentLogger()->outputDebugString("u: " + String(u[1][l]));
			juce::Logger::getCurrentLogger()->outputDebugString("uPrev: " + String(u[2][l]));
		}
		init = false;
	}
	*/
    for (int idx = 2; idx < N-1; ++idx)
    {
        idx_p1 = idx + 1;
        idx_p2 = idx + 2;
        idx_m1 = idx - 1;
        idx_m2 = idx - 2;

        I_B_J_B = I_B_J_B + I_B[idx] * J_B[idx];
        I_B_u = I_B_u + I_B[idx] * u[1][idx];
        I_B_uPrev = I_B_uPrev + I_B[idx] * u[2][idx];

        I_B_dxx_u = I_B_dxx_u + I_B[idx] * (u[1][idx_p1] - 2. * u[1][idx] + u[1][idx_m1]) * (1 / (h * h));
        I_B_dxx_uPrev = I_B_dxx_uPrev + I_B[idx] * (u[2][idx_p1] - 2. * u[2][idx] + u[2][idx_m1]) * (1 / (h * h));

        I_B_dxxxx_u = I_B_dxxxx_u
            + I_B[idx] * (u[1][idx_p2] - 4. * u[1][idx_p1] + 6. * u[1][idx] - 4. * u[1][idx_m1] + u[1][idx_m2]) * (1 / (h * h * h * h));
    }

	/*
	if (counter == 0)
	{
		juce::Logger::getCurrentLogger()->outputDebugString("I_B_J_B: " + String(I_B_J_B));
		juce::Logger::getCurrentLogger()->outputDebugString("I_B_u: " + String(I_B_u));
		juce::Logger::getCurrentLogger()->outputDebugString("I_B_uPrev: " + String(I_B_uPrev));
		juce::Logger::getCurrentLogger()->outputDebugString("I_B_dxx_u: " + String(I_B_dxx_u));
		juce::Logger::getCurrentLogger()->outputDebugString("I_B_dxx_uPrev: " + String(I_B_dxx_uPrev));
		juce::Logger::getCurrentLogger()->outputDebugString("I_B_dxxxx_u: " + String(I_B_dxxxx_u));
		init = false;
	}
	*/
    
    double q = (-2 / k) * (1 / k) * (I_B_u - I_B_uPrev)
                + (2 / k) * vB_var
                + 2 * sig0 * vB_var
                - (c * c) * I_B_dxx_u
                + K * K * I_B_dxxxx_u
                - (2 * sig1) * (1 / k) * I_B_dxx_u
                + (2 * sig1) * (1 / k) * I_B_dxx_uPrev;

    if (fN_var > 0)
    {
        double eps = 1;
        //w_rnd_last = -1 + (1 - (-1)).*rand(1);

        int iter_check = 0;
        vrel_last = -vB; // should this be vB_var ?
        // Newton-Raphson iterative scheme
        while ((eps > tol))
        {
            ++iter_check;

            g1 = (2/k + 2*sig0) * vrel_last + I_B_J_B * fN_var * sqrt(2 * a) * exp(0.5) * vrel_last * exp(- a *vrel_last * vrel_last) + q;

            g1_d_vrel = I_B_J_B * (fN_var * sqrt(2 * a) * exp(0.5) * exp(-a * vrel_last * vrel_last) * (1 + 2 * a * vrel_last * vrel_last)) / (rho * A) + 2 * sig0 + 2 / k;
            
            vrel = vrel_last - g1 / g1_d_vrel;
            eps = abs(vrel - vrel_last);

            vrel_last = vrel;
            if (iter_check == 99)
                //Logger::getCurrentLogger()->outputDebugString("iter_check: (" + String(iter_check) + ")");
                break;

        }
        vrel = vrel_last;

        Fbow = fN_var * exp(0.5) * sqrt(2 * a) * vrel *exp(-a * vrel * vrel);

        //Fbow = fN_var * A_fr_model * vrel * exp(-sig * vrel * vrel);
        phi_vrel = exp(0.5) * sqrt(2 * a) * vrel * exp(-a * vrel * vrel);
    }
    else
    {
        Fbow = 0;
        phi_vrel = exp(0.5) * sqrt(2 * a) * vrel * exp(-a * vrel * vrel);
        vrel_last = -vB;
    }

	/*
	if (counter == 0)
	{
		juce::Logger::getCurrentLogger()->outputDebugString("Fbow: " + String(Fbow));
	}
    */

    // Update equations
    for (int idx = 2; idx < N-1 ; ++idx)
    {
        idx_p1 = idx + 1;
        idx_p2 = idx + 2;
        idx_m1 = idx - 1;
        idx_m2 = idx - 2;

         // maybe you can reused the definitions from above I_S_etc

        u[0][idx] = (k * k / (1 + sig0 * k)) * ((c * c / (h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1])
                + (2 * sig1 / (k * h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1] - u[2][idx_p1] + 2 * u[2][idx] - u[2][idx_m1])
                - (K * K) * (1 / (h * h * h * h)) * (u[1][idx_p2] - 4 * u[1][idx_p1] + 6 * u[1][idx] - 4 * u[1][idx_m1] + u[1][idx_m2])
                - J_B[idx] * Fbow / (rho * A)
                + (2 / (k * k)) * u[1][idx] - (1 - sig0 * k) * u[2][idx] / (k * k));

		//u[0][idx] = (k * k / (1 + sig0 * k)) *  ((c * c / (h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1])
		//	+ (2 * sig1 / (k * h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1] - u[2][idx_p1] + 2 * u[2][idx] - u[2][idx_m1]));

		//u[0][idx] = (2 * sig1 / (k * h * h)) ;
		/*
		if (counter == 0)
		{
			juce::Logger::getCurrentLogger()->outputDebugString("u[0][idx]: " + String(u[0][idx]));
			juce::Logger::getCurrentLogger()->outputDebugString("sig1: " + String(sig1));
			juce::Logger::getCurrentLogger()->outputDebugString("k: " + String(k));
			juce::Logger::getCurrentLogger()->outputDebugString("h: " + String(h));
		}
		*/
    }
    
    u[0][1] = (k * k / (1 + sig0 * k)) * ((c * c / (h * h)) * (u[1][2] - 2 * u[1][1] + u[1][0])
        + (2 * sig1 / (k * h * h)) * (u[1][2] - 2 * u[1][1] + u[1][0] - u[2][2] + 2 * u[2][1] - u[2][0])
        - (K * K) * (1 / (h * h * h * h)) * (u[1][3] - 4 * u[1][2] + 6 * u[1][1] - 4 * u[1][0] - u[1][1])
        - J_B[1] * Fbow / (rho * A)
        + (2 / (k * k)) * u[1][1] - (1 - sig0 * k) * u[2][1] / (k * k));

    //int idx = N - 2; // this was a bug i think
    int idx = N-1;
    idx_p1 = idx + 1;
    idx_p2 = idx + 2;
    idx_m1 = idx - 1;
    idx_m2 = idx - 2;

    u[0][idx] = (k * k / (1 + sig0 * k)) * ((c * c / (h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1])
        + (2 * sig1 / (k * h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1] - u[2][idx_p1] + 2 * u[2][idx] - u[2][idx_m1])
        - (K * K) * (1 / (h * h * h * h)) * (-u[1][idx] - 4 * u[1][idx_p1] + 6 * u[1][idx] - 4 * u[1][idx_m1] + u[1][idx_m2])
        - J_B[idx] * Fbow / (rho * A)
        + (2 / (k * k)) * u[1][idx] - (1 - sig0 * k) * u[2][idx] / (k * k));

    uTmp = u[2];
    u[2] = u[1];
    u[1] = u[0];
    u[0] = uTmp;

	// u[0] -> uNext, u[1] -> uCurrent, u[2] -> uPrev ! 

	F_fr = Fbow; // check this ! 
    
    int idx_out = round(N * 0.7); // THIS CAN BE VARIED !! MAYBE USE AN INTERPOLATION FCT...
    double out = u[0][idx_out];

	/*
	if (counter == 0)
	{
		for (int l = 0; l <= N; ++l)
		{
			juce::Logger::getCurrentLogger()->outputDebugString("uNext: " + String(u[0][l]));
			juce::Logger::getCurrentLogger()->outputDebugString("u: " + String(u[1][l]));
			juce::Logger::getCurrentLogger()->outputDebugString("uPrev: " + String(u[2][l]));
		}
		init = false;
	}
	*/
	counter++;
    return out;
}
    
