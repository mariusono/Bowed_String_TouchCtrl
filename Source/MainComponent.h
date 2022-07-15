#pragma once

#include <JuceHeader.h>
#include "ViolinString.h"
//#include "SenselWrapper.h"


#ifdef  _WIN64
#pragma warning (disable:4996)
#endif


#if defined(WIN32)
# include <windows.h>
# include <conio.h>
#else
# include "conio.h"
# include <unistd.h>
# define Sleep(x) usleep((x) * 1000)
#endif

#include <HD/hd.h>
#include <HDU/hduVector.h>
#include <HDU/hduError.h>
#include <stdio.h>



//==============================================================================
/*
	This component lives inside our window, and this is where you should put all
	your controls and content.
*/
class MainComponent : public juce::AudioAppComponent,
	public juce::HighResolutionTimer,
	public juce::Timer,
	//public juce::MultiTimer,
	public juce::Button::Listener,
	public juce::Slider::Listener,
	public juce::MouseListener
{
public:
	//==============================================================================
	MainComponent();
	~MainComponent() override;

	//==============================================================================
	void prepareToPlay(int samplesPerBlockExpected, double sampleRate) override;
	void getNextAudioBlock(const juce::AudioSourceChannelInfo& bufferToFill) override;
	void releaseResources() override;

	//==============================================================================
	void paint(juce::Graphics& g) override;
	void paintOverChildren(juce::Graphics& g) override;
	void resized() override;


	// function to draw the state of the string
//    Path visualiseState (Graphics& g, double visualScaling, std::vector<std::vector<double>> uVecs, int Npoints, int strNo);
//    Path visualiseState_crossSec (Graphics& g, double visualScaling, std::vector<std::vector<double>> uVecs, int Npoints, int strNo);

	Path visualiseState(Graphics& g, float visualScaling, ViolinString* string, float startPos_x);
	Path visualiseState_crossSec(Graphics& g, double visualScaling, ViolinString* string, double startPos_x, double startPos_y, std::vector<double> direction);



	// MultiTimer callback
	void timerCallback() override;

	// Hi resolution timer callback
	void hiResTimerCallback() override;

	int idx_mass = 0;
	double wheelDeltaY = 0.0;

	// Slider stuff
	void sliderValueChanged(Slider* slider) override;

	// Button callback
	void buttonClicked(juce::Button* button) override; // [2]

	int timerGraphics = 1;
	int timerLogger = 2;

private:
	//==============================================================================

	double fs;
	double bufferSize;


	float minOut;
	float maxOut;
	int numStrings;
	int octave;
	int polyphony;

	double globalCurrentSample;

	//bool includeHaptics = true;


	// For graphics
	float opa_level = 1.0; // Opacity
	std::vector<double> gPositionRaw;
	float widthScreen = 800;
	float heightScreen = 600;

	//std::vector<double> gBowEnd_1;
	//std::vector<double> gBowEnd_2;




	// For something else I think
	bool flagMouseUp = true;

	// SLIDER STUFF
	Slider globalDampingSlider;
	Label globalDampingLabel;

	Slider dampingSlider;
	Label dampingLabel;

	Slider freqDampingSlider;
	Label freqDampingLabel;

	Slider frParamSlider;
	Label frParamLabel;

	Slider volumeSlider;
	Label volumeLabel;

	// BUTTON STUFF
	TextButton hapticFrictionButton;
	Label hapticFrictionLabel;

	double gDampingVal;
	double gFrParam = 80;
	double gSig0 = 1;
	double gSig1 = 0.008;
	double gVolume;
	double gStickFact;


	// Phantom stuff
	HHD hHD; // Device initialization instance
	hduVector3Dd gPositionCallback; // Global position vector
	//hduVector3Dd gPositionPrev; // Global position vector

	hduVector3Dd gBowEnd_1;
	hduVector3Dd gBowEnd_2;


	/*
	// Sensel Stuff
	OwnedArray<Sensel> sensels;
	int amountOfSensels = 1;
	*/

	double hiResCallbackFreq;

	double gIdx_map;
	bool flagBow;

	// STRING INSTANCES
	OwnedArray<ViolinString> ViolinStrings; // a vector of unique pointers maaaybe it's this: std::vector<std::unique_ptr<ViolinString>>
	std::vector<ViolinString*> activeStrings; // why the * ? Pointer !

	// More stuff..
	double N;

	std::vector<double> preOutputVec;
	std::vector<bool> keyDownVec;
	std::vector<bool> bowOnVec;

	std::vector<double> gPositionScaled;



	std::vector<double> gPosition;
	std::vector<double> gPositionPrev;

	/// Parameters you may want to modulate:
	double FB_Max;

	//double maxVb;
	//double maxFB;

	double bp; // in percentage
	int bP;
	double alpha_bow;

	// Bowing interpolation fct and spreading fct
	//std::vector<double> I_B; // interpolant grid for bowing pos
	//std::vector<double> J_B; // spreading function grid for bowing pos

	double x_inp_var;
	double y_inp_var;




	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainComponent)
};
