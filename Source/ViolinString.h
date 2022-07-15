/*
  ==============================================================================

    ViolinString.h
    Created: 31 Mar 2022 8:19:26pm
    Author:  Marius Onofrei

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
*/
class ViolinString  : public juce::Component
{
public:
    ViolinString(double freq, double fs);
    ~ViolinString() override;


    void paint (juce::Graphics&) override;
    void paint_crossSec (juce::Graphics&);
    void resized() override;
    
    // function to draw the state of the string
    Path visualiseState (Graphics& g, double visualScaling);
    Path visualiseState_crossSec (Graphics& g, double visualScaling);

    
    double process();
//    void newtonRaphson();
    //void updateStates();

    void setBow (bool val) { isBowing = val; };
    void setVb (double val) { vB_var = val; }
    void setFb(double val) { fN_var = val; }
    void setSig0(double val) { sig0 = val; }
    void setSig1(double val) { sig1 = val; }
    void setFrParam(double val) { a = val; }
    void setStickFact(double val) { stickFact = val; }
    void setI_B(std::vector<double> val) { I_B = val; }
    void setJ_B(std::vector<double> val) { J_B = val; }
    void setPaintFctIdx(int val) { paintFctIdx = val; }

    double getSig0() { return sig0; }
    double getSig1() { return sig1; }
    double getVb() { return vB_var; }
    double getFb() { return fN_var; }
    double getFrParam() { return a; }
    double getStickFact() { return stickFact;  }
    double getN() { return N;  }
    double get_h() { return h;  }
	double getF_fr() { return F_fr; }
    std::vector<double> getI_B() { return I_B;  }
    std::vector<double> getJ_B() { return J_B;  }
    int get_paintFctIdx() { return paintFctIdx;  }

    //std::vector<double*> get_u() {return u;}
    std::vector<double*> get_u() {return u;}
	std::vector<std::vector<double>> get_uVecs() { return uVecs; }
    
    bool isActive() { return active; };
	bool isBowOn() { return isBowing; };
    void activate() { active = true; };
    void deactivate() { active = false; };
private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ViolinString)
    
    // paint stuff
    int paintFctIdx;
    
    // init params
    double fs, freq, mass, sig0, sig1, vB, A_NR, tol, a, k, h, c;
    double FB;
    double L,rho,r,T,E;
    double A,I,K;
    double inp_bow_pos_x; // in perc
    double N;
    
    double fN_var, vB_var;
    double Fbow;
    double phi_vrel;
    
    double stickFact;
    
    // BOWING POSITION INIT // this might go in MainComponent..
    /// Parameters you may want to modulate:
    double bp = 0.25; // in percentage
    int bP;
    double alpha_bow;
    
    // Bowing interpolation fct and spreading fct
    std::vector<double> I_B; // interpolant grid for bowing pos
    std::vector<double> J_B; // spreading function grid for bowing pos
    
    // pointers to STRING states
    std::vector<double*> u;
    
    // Vector states
    std::vector<std::vector<double>> uVecs;
    double* uTmp = nullptr;
    
    
    // INITIALIZE PARAMS FOR NEWTON-RAPHSON
    double vrel;
    double vrel_last;
    double g1;
    double g1_d_vrel;
    
    
    // process params
    double F_fr;
    
    // NR params
    double eps, numerator, denominator;
    
    // State params
    double u_n, u_nm1, u_np1; // velocity and displacement of mass

    
//    gamma, k, s0, s1, B, kappa, h, N, lambdaSq, muSq, kOh, gOh, a, BM, Vb, Fb, pickup, tol, q, qPrev, bp, b, eps;
    bool isBowing = false;
//    std::vector<double> u, uPrev, uNext;
    bool active = false;
    
    unsigned long countGlobal;
//    unsigned long t = 0;

	// FOR DEBUG STUFF
	int counter = 0;
	bool init = true;
};
