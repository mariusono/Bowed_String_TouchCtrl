#include "MainComponent.h"

// OTHER FUNCTIONS

double linearMapping(double rangeIn_top, double rangeIn_bottom, double rangeOut_top, double rangeOut_bottom, double value) {
	double newValue = rangeOut_bottom + ((rangeOut_top - rangeOut_bottom) * (value - rangeIn_bottom) / (rangeIn_top - rangeIn_bottom));
	return newValue;
}


double exponentialMapping(double rangeIn_top, double rangeIn_bottom, double rangeOut_top, double rangeOut_bottom, double fac, double value)
{
	// make sure values passed to function are within the rangeIn_bottom rangeIn_top interval !!
	// maybe add an error exception here..
	// first map value to rangeIn to 0 - 1
	double valueMapped = 0.0 + ((1.0 - 0.0) * (value - rangeIn_bottom) / (rangeIn_top - rangeIn_bottom));

	// map to an exponential curve between 0 and 1 with a factor fac
	double mapToExp = (exp(valueMapped * fac) - 1) / (exp(fac) - 1);

	// map back to desired output range
	double newValue = rangeOut_bottom + ((rangeOut_top - rangeOut_bottom) * (mapToExp - 0.0) / (1.0 - 0.0));

	return newValue;
}

double clamp(double in, double min, double max)
{
	if (in > max)
		return max;
	else if (in < min)
		return min;
	else
		return in;
}

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}






double distanceBetweenBowAndString_M1(hduVector3Dd string_1_p1, hduVector3Dd string_1_p2, hduVector3Dd bow_p1, hduVector3Dd bow_p2)
{
	// // METHOD 1
// Calculate distance between line - negative values for dist means it touched ! 

	hduVector3Dd e1; // dir string
	hduVector3Dd e2; // dir bow

	e1[0] = string_1_p2[0] - string_1_p1[0];
	e1[1] = string_1_p2[1] - string_1_p1[1];
	e1[2] = string_1_p2[2] - string_1_p1[2];

	e2[0] = bow_p2[0] - bow_p1[0];
	e2[1] = bow_p2[1] - bow_p1[1];
	e2[2] = bow_p2[2] - bow_p1[2];

	if (sgn(e2[1]) != sgn(e1[1]))
	{
		e2[0] = bow_p1[0] - bow_p2[0];
		e2[1] = bow_p1[1] - bow_p2[1];
		e2[2] = bow_p1[2] - bow_p2[2];
	}

	hduVector3Dd n; // normal vec , cross product between e1 and e2
	n[0] = e1[1] * e2[2] - e1[2] * e2[1];
	n[1] = e1[2] * e2[0] - e1[0] * e2[2];
	n[2] = e1[0] * e2[1] - e1[1] * e2[0];

	double n_norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

	//float dist = 0;
	double dist = (n[0] * (string_1_p1[0] - bow_p1[0]) + n[1] * (string_1_p1[1] - bow_p1[1]) + n[2] * (string_1_p1[2] - bow_p1[2])) / n_norm;

	return dist;
}

struct intersectionPoints
{
	double xC, yC, zC, xD, yD, zD;
};

intersectionPoints distanceBetweenBowAndString_M2(hduVector3Dd string_1_p1, hduVector3Dd string_1_p2, hduVector3Dd bow_p1, hduVector3Dd bow_p2)
{
	// // METHOD 2
// Calculate distance between line - negative values for dist means it touched ! 

	intersectionPoints pointsIntersect;

	double xA = string_1_p1[0];
	double yA = string_1_p1[1];
	double zA = string_1_p1[2];
	double xB = bow_p1[0];
	double yB = bow_p1[1];
	double zB = bow_p1[2];


	double ax = string_1_p2[0] - string_1_p1[0];
	double ay = string_1_p2[1] - string_1_p1[1];
	double az = string_1_p2[2] - string_1_p1[2];

	double bx = bow_p2[0] - bow_p1[0];
	double by = bow_p2[1] - bow_p1[1];
	double bz = bow_p2[2] - bow_p1[2];

	double xC = (ay*ay*bx*bx*xA + ax * ax*by*by*xB + az * az*bx*bx*xA + ax * ax*bz*bz*xB + ay * ay*bz*bz*xA + az * az*by*by*xA - ax * ay*bx*bx*yA + ax * ay*bx*bx*yB - ax * ay*bz*bz*yA + ax * ay*bz*bz*yB + ax * ax*bx*by*yA - ax * az*bx*bx*zA - ax * ax*bx*by*yB + ax * az*bx*bx*zB - ax * az*by*by*zA + ax * az*by*by*zB + ax * ax*bx*bz*zA - ax * ax*bx*bz*zB - ax * ay*bx*by*xA - ax * ay*bx*by*xB - ax * az*bx*bz*xA - ax * az*bx*bz*xB - 2 * ay*az*by*bz*xA + ax * az*by*bz*yA - ax * az*by*bz*yB + ax * ay*by*bz*zA - ax * ay*by*bz*zB) / (ax*ax*by*by + ax * ax*bz*bz - 2 * ax*ay*bx*by - 2 * ax*az*bx*bz + ay * ay*bx*bx + ay * ay*bz*bz - 2 * ay*az*by*bz + az * az*bx*bx + az * az*by*by);
	double yC = (ax*ax*by*by*yA + ax * ax*bz*bz*yA + ay * ay*bx*bx*yB + az * az*bx*bx*yA + az * az*by*by*yA + ay * ay*bz*bz*yB - ax * ay*by*by*xA + ax * ay*by*by*xB - ax * ay*bz*bz*xA + ax * ay*bz*bz*xB + ay * ay*bx*by*xA - ay * ay*bx*by*xB - ay * az*bx*bx*zA + ay * az*bx*bx*zB - ay * az*by*by*zA + ay * az*by*by*zB + ay * ay*by*bz*zA - ay * ay*by*bz*zB + ay * az*bx*bz*xA - ay * az*bx*bz*xB - ax * ay*bx*by*yA - ax * ay*bx*by*yB - 2 * ax*az*bx*bz*yA - ay * az*by*bz*yA - ay * az*by*bz*yB + ax * ay*bx*bz*zA - ax * ay*bx*bz*zB) / (ax*ax*by*by + ax * ax*bz*bz - 2 * ax*ay*bx*by - 2 * ax*az*bx*bz + ay * ay*bx*bx + ay * ay*bz*bz - 2 * ay*az*by*bz + az * az*bx*bx + az * az*by*by);
	double zC = (ax*ax*by*by*zA + ay * ay*bx*bx*zA + ax * ax*bz*bz*zA + ay * ay*bz*bz*zA + az * az*bx*bx*zB + az * az*by*by*zB - ax * az*by*by*xA + ax * az*by*by*xB - ax * az*bz*bz*xA + ax * az*bz*bz*xB - ay * az*bx*bx*yA + ay * az*bx*bx*yB + az * az*bx*bz*xA - ay * az*bz*bz*yA - az * az*bx*bz*xB + ay * az*bz*bz*yB + az * az*by*bz*yA - az * az*by*bz*yB + ay * az*bx*by*xA - ay * az*bx*by*xB + ax * az*bx*by*yA - ax * az*bx*by*yB - 2 * ax*ay*bx*by*zA - ax * az*bx*bz*zA - ax * az*bx*bz*zB - ay * az*by*bz*zA - ay * az*by*bz*zB) / (ax*ax*by*by + ax * ax*bz*bz - 2 * ax*ay*bx*by - 2 * ax*az*bx*bz + ay * ay*bx*bx + ay * ay*bz*bz - 2 * ay*az*by*bz + az * az*bx*bx + az * az*by*by);

	double xD = (ay*ay*bx*bx*xA + ax * ax*by*by*xB + az * az*bx*bx*xA + ax * ax*bz*bz*xB + ay * ay*bz*bz*xB + az * az*by*by*xB - ax * ay*bx*bx*yA + ax * ay*bx*bx*yB + ax * ax*bx*by*yA - ax * az*bx*bx*zA - ax * ax*bx*by*yB + ax * az*bx*bx*zB + az * az*bx*by*yA - az * az*bx*by*yB + ax * ax*bx*bz*zA - ax * ax*bx*bz*zB + ay * ay*bx*bz*zA - ay * ay*bx*bz*zB - ax * ay*bx*by*xA - ax * ay*bx*by*xB - ax * az*bx*bz*xA - ax * az*bx*bz*xB - 2 * ay*az*by*bz*xB - ay * az*bx*bz*yA + ay * az*bx*bz*yB - ay * az*bx*by*zA + ay * az*bx*by*zB) / (ax*ax*by*by + ax * ax*bz*bz - 2 * ax*ay*bx*by - 2 * ax*az*bx*bz + ay * ay*bx*bx + ay * ay*bz*bz - 2 * ay*az*by*bz + az * az*bx*bx + az * az*by*by);
	double yD = (ax*ax*by*by*yA + ay * ay*bx*bx*yB + ax * ax*bz*bz*yB + az * az*bx*bx*yB + az * az*by*by*yA + ay * ay*bz*bz*yB - ax * ay*by*by*xA + ax * ay*by*by*xB + ay * ay*bx*by*xA - ay * ay*bx*by*xB + az * az*bx*by*xA - az * az*bx*by*xB - ay * az*by*by*zA + ay * az*by*by*zB + ax * ax*by*bz*zA - ax * ax*by*bz*zB + ay * ay*by*bz*zA - ay * ay*by*bz*zB - ax * az*by*bz*xA + ax * az*by*bz*xB - ax * ay*bx*by*yA - ax * ay*bx*by*yB - 2 * ax*az*bx*bz*yB - ay * az*by*bz*yA - ay * az*by*bz*yB - ax * az*bx*by*zA + ax * az*bx*by*zB) / (ax*ax*by*by + ax * ax*bz*bz - 2 * ax*ay*bx*by - 2 * ax*az*bx*bz + ay * ay*bx*bx + ay * ay*bz*bz - 2 * ay*az*by*bz + az * az*bx*bx + az * az*by*by);
	double zD = (ax*ax*by*by*zB + ax * ax*bz*bz*zA + ay * ay*bx*bx*zB + ay * ay*bz*bz*zA + az * az*bx*bx*zB + az * az*by*by*zB - ax * az*bz*bz*xA + ax * az*bz*bz*xB + ay * ay*bx*bz*xA - ay * ay*bx*bz*xB + az * az*bx*bz*xA - ay * az*bz*bz*yA - az * az*bx*bz*xB + ay * az*bz*bz*yB + ax * ax*by*bz*yA - ax * ax*by*bz*yB + az * az*by*bz*yA - az * az*by*bz*yB - ax * ay*by*bz*xA + ax * ay*by*bz*xB - ax * ay*bx*bz*yA + ax * ay*bx*bz*yB - 2 * ax*ay*bx*by*zB - ax * az*bx*bz*zA - ax * az*bx*bz*zB - ay * az*by*bz*zA - ay * az*by*bz*zB) / (ax*ax*by*by + ax * ax*bz*bz - 2 * ax*ay*bx*by - 2 * ax*az*bx*bz + ay * ay*bx*bx + ay * ay*bz*bz - 2 * ay*az*by*bz + az * az*bx*bx + az * az*by*by);

	pointsIntersect.xC = xC;
	pointsIntersect.yC = yC;
	pointsIntersect.zC = zC;
	pointsIntersect.xD = xD;
	pointsIntersect.yD = yD;
	pointsIntersect.zD = zD;

	return pointsIntersect;
}





int find_if_point_on_string(intersectionPoints pointsIntersect, hduVector3Dd string_1_p1, hduVector3Dd string_1_p2)
{
	// CHECK TO SEE IF C IS ON THE STRING ! 
	//dotproduct = dot(p12 - p11, C - p11);
	hduVector3Dd C_min_string_1_p1; // normal vec , cross product between e1 and e2

	C_min_string_1_p1[0] = pointsIntersect.xC - string_1_p1[0];
	C_min_string_1_p1[1] = pointsIntersect.yC - string_1_p1[1];
	C_min_string_1_p1[2] = pointsIntersect.zC - string_1_p1[2];

	hduVector3Dd e1;
	e1[0] = string_1_p2[0] - string_1_p1[0];
	e1[1] = string_1_p2[1] - string_1_p1[1];
	e1[2] = string_1_p2[2] - string_1_p1[2];


	double dotproduct = e1[0] * C_min_string_1_p1[0] + e1[1] * C_min_string_1_p1[1] + e1[2] * C_min_string_1_p1[2];
	double squareddist = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);

	int pointOnString;

	if ((dotproduct < 0) || (dotproduct > squareddist))
	{
		pointOnString = 0;
	}
	else
	{
		pointOnString = 1;
	}
	//Logger::getCurrentLogger()->outputDebugString("pointOnString: (" + String(pointOnString) + ")");

	return pointOnString;

}


int find_if_point_on_bow(intersectionPoints pointsIntersect, hduVector3Dd bow_p1, hduVector3Dd bow_p2)
{
	// CHECK TO SEE IF C IS ON THE STRING ! 
	//dotproduct = dot(p12 - p11, C - p11);
	hduVector3Dd D_min_bow_p1; // normal vec , cross product between e1 and e2

	D_min_bow_p1[0] = pointsIntersect.xD - bow_p1[0];
	D_min_bow_p1[1] = pointsIntersect.yD - bow_p1[1];
	D_min_bow_p1[2] = pointsIntersect.zD - bow_p1[2];

	hduVector3Dd e1;
	e1[0] = bow_p2[0] - bow_p1[0];
	e1[1] = bow_p2[1] - bow_p1[1];
	e1[2] = bow_p2[2] - bow_p1[2];


	double dotproduct = e1[0] * D_min_bow_p1[0] + e1[1] * D_min_bow_p1[1] + e1[2] * D_min_bow_p1[2];
	double squareddist = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);

	int pointOnString;

	if ((dotproduct < 0) || (dotproduct > squareddist))
	{
		pointOnString = 0;
	}
	else
	{
		pointOnString = 1;
	}
	//Logger::getCurrentLogger()->outputDebugString("pointOnString: (" + String(pointOnString) + ")");

	return pointOnString;

}


hduVector3Dd intersect_points_dir(intersectionPoints points)
{
	hduVector3Dd e_CD;
	double d_alt = sqrt((points.xD - points.xC) * (points.xD - points.xC) + (points.yD - points.yC) * (points.yD - points.yC) + (points.zD - points.zC) * (points.zD - points.zC));
	e_CD[0] = (points.xD - points.xC) / d_alt;
	e_CD[1] = (points.yD - points.yC) / d_alt;
	e_CD[2] = (points.zD - points.zC) / d_alt;

	return e_CD;
}

hduVector3Dd get_fr_force_direction(hduVector3Dd e1, hduVector3Dd e_CD)
{
	hduVector3Dd e1_e_CD_cross;

	e1_e_CD_cross[0] = e1[1] * e_CD[2] - e1[2] * e_CD[1];
	e1_e_CD_cross[1] = e1[2] * e_CD[0] - e1[0] * e_CD[2];
	e1_e_CD_cross[2] = e1[0] * e_CD[1] - e1[1] * e_CD[0];

	double e1_e_CD_cross_norm = sqrt(e1_e_CD_cross[0] * e1_e_CD_cross[0] + e1_e_CD_cross[1] * e1_e_CD_cross[1] + e1_e_CD_cross[2] * e1_e_CD_cross[2]);

	double e1_norm = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
	double e_CD_norm = sqrt(e_CD[0] * e_CD[0] + e_CD[1] * e_CD[1] + e_CD[2] * e_CD[2]);

	hduVector3Dd direction;
	direction[0] = e1_e_CD_cross[0] / (e_CD_norm*e1_norm);
	direction[1] = e1_e_CD_cross[1] / (e_CD_norm*e1_norm);
	direction[2] = e1_e_CD_cross[2] / (e_CD_norm*e1_norm);

	return direction;
}



// DEBUG FLAGS
bool soundFlag = true;
bool hapticStringStiffnessFlag = true;
bool hapticFrictionFlag = true;

// BUTTON FLAG
bool includeHaptics = true;

// global variables:
int gNumStrings = 4;

//double lengthBow = 200; // units in touch device space
double lengthBow = 450; // units in touch device space


std::vector<double> gOutputStrings; // does the JUCE output have to be floats ? 
double gOutput = 0; // global output val
hduVector3Dd gF;

//float gF_fr = 0; // global fr force val
std::vector<double> gF_fr(gNumStrings);

// Bound coord for Touch device
double zMin = -35.0;
double zMax = 35.0;
double yMin = 0.0;


// moar global stuff

double maxVb = 0.4;
//double maxFB = 2.0;
double maxFB = 1.0;
double FB = 0;

std::vector<double> FB_vec(gNumStrings);


// TOUCH DEVICE STUFF

std::vector<hduVector3Dd> string_p1(gNumStrings);
std::vector<hduVector3Dd> string_p2(gNumStrings);

hduVector3Dd bow_p1; // point to the left of gimbal center
hduVector3Dd bow_p2; // point to the right of gimbal center

std::vector<std::vector<double>> gDirection_vec(gNumStrings, std::vector<double>(3, 0));
std::vector<std::vector<double>> gDirection_vec_print(gNumStrings, std::vector<double>(3, 0));


std::vector<double> gActiveString(gNumStrings);
std::vector<hduVector3Dd> f_vec(gNumStrings); // force for each string.. for mapping FB individually

/* not needed actually
// PROJECTIONS FOR PAINT
std::vector<double> bowProj_p1_xz_scr1(3);
std::vector<double> bowProj_p2_xz_scr1(3);
std::vector<double> bowProj_p1_xy_scr2(3);
std::vector<double> bowProj_p2_xy_scr2(3);

std::vector<std::vector<double>> string_p1_xz_scr1(gNumStrings, std::vector<double>(3));
std::vector<std::vector<double>> string_p2_xz_scr1(gNumStrings, std::vector<double>(3));
std::vector<std::vector<double>> string_p1_xy_scr2(gNumStrings, std::vector<double>(3));
std::vector<std::vector<double>> string_p2_xy_scr2(gNumStrings, std::vector<double>(3));
*/

HDCallbackCode HDCALLBACK StringCallback(void *data)
{

	hdBeginFrame(hdGetCurrentDevice());


	HDdouble usableWorkspace[6];
	HDdouble maxWorkspace[6];
	hdGetDoublev(HD_USABLE_WORKSPACE_DIMENSIONS, usableWorkspace);
	hdGetDoublev(HD_MAX_WORKSPACE_DIMENSIONS, maxWorkspace);

	//Logger::getCurrentLogger()->outputDebugString("usableWorkspace[0]: (" + String(usableWorkspace[0]) + ") usableWorkspace[1]: (" + String(usableWorkspace[1]) + ") usableWorkspace[2]: (" + String(usableWorkspace[2]) + ") usableWorkspace[3]: (" + String(usableWorkspace[3]) + ") usableWorkspace[4]: (" + String(usableWorkspace[4]) + ") usableWorkspace[5]: (" + String(usableWorkspace[5]) + ")");

	//string_p1_xz_scr1;
	int ana = 3;

	HDdouble frame[16];
	hdGetDoublev(HD_CURRENT_TRANSFORM, frame);

	hduVector3Dd position;
	hdGetDoublev(HD_CURRENT_POSITION, position); // x: left - right, y: bottom - up, z: away from u - towards u

	hduVector3Dd gimbalAngles;
	hdGetDoublev(HD_CURRENT_GIMBAL_ANGLES, gimbalAngles);

	double theta_1 = -asin(frame[2]);
	//theta_1 = double_Pi - theta_1;
	double psi_1 = atan2(frame[6] / cos(theta_1), frame[10] / cos(theta_1));
	double phi_1 = atan2(frame[1] / cos(theta_1), frame[0] / cos(theta_1));
	double angle_with_xz_plane = phi_1;
	if (phi_1 - (95.0*double_Pi / 180.0) < (-75.0 * double_Pi / 180.0))
	{
		angle_with_xz_plane = (-75.0 * double_Pi / 180.0);
	}
	else if (phi_1 - (95.0*double_Pi / 180.0) > (75.0 * double_Pi / 180.0))
	{
		angle_with_xz_plane = (75.0 * double_Pi / 180.0);
	}
	else
	{
		angle_with_xz_plane = (phi_1 - (95.0*double_Pi / 180.0));
	}

	phi_1 = clamp(phi_1, -60.0, 60.0);

	//Logger::getCurrentLogger()->outputDebugString("theta_1: (" + String(theta_1 * 180 / double_Pi) + ") psi_1: (" + String(psi_1 * 180 / double_Pi) + ") phi_1: (" + String(phi_1 * 180 / double_Pi) + ")");
	//Logger::getCurrentLogger()->outputDebugString("frame[0]: (" + String(frame[0]) + ") frame[1]: (" + String(frame[1]) + ") frame[2]: (" + String(frame[2]) + ")");



	//Logger::getCurrentLogger()->outputDebugString("frame[12]: (" + String(frame[12]) + ") frame[13]: (" + String(frame[13]) + ") frame[14]: (" + String(frame[14]) + ")");
	//Logger::getCurrentLogger()->outputDebugString("position[0]: (" + String(position[0]) + ") position[1]: (" + String(position[1]) + ") position[2]: (" + String(position[2]) + ")");
	//Logger::getCurrentLogger()->outputDebugString("gimbalAngles[0]: (" + String(gimbalAngles[0] * 180 / double_Pi) + ") gimbalAngles[1]: (" + String(gimbalAngles[1] * 180 / double_Pi) + ") gimbalAngles[2]: (" + String(gimbalAngles[2] * 180 / double_Pi) + ")");


	//Logger::getCurrentLogger()->outputDebugString("theta_1: (" + String(theta_1 * 180 / double_Pi) + ") psi_1: (" + String(psi_1 * 180 / double_Pi) + ") phi_1: (" + String(phi_1 * 180 / double_Pi) + ") gimbalAngles[2]: (" + String(gimbalAngles[2] * 180 / double_Pi) + ")");


	string_p1[0][0] = -75.0;
	string_p1[0][1] = -10.0;
	string_p1[0][2] = -90.0;

	string_p2[0][0] = -75.0;
	string_p2[0][1] = -10.0;
	string_p2[0][2] = 90.0;

	string_p1[1][0] = -25.0;
	string_p1[1][1] = 25.0;
	string_p1[1][2] = -90.0;

	string_p2[1][0] = 25.0;
	string_p2[1][1] = 25.0;
	string_p2[1][2] = 90.0;

	string_p1[2][0] = 25.0;
	string_p1[2][1] = 25.0;
	string_p1[2][2] = -90.0;

	string_p2[2][0] = 25.0;
	string_p2[2][1] = 25.0;
	string_p2[2][2] = 90.0;

	string_p1[3][0] = 75.0;
	string_p1[3][1] = -10.0;
	string_p1[3][2] = -90.0;

	string_p2[3][0] = 75.0;
	string_p2[3][1] = -10.0;
	string_p2[3][2] = 90.0;

	double d = lengthBow / 2;

	/*
	bow_p1[0] = -d * cos(-(phi_1 - double_Pi / 2))*sin(double_Pi / 2) + position[0];
	bow_p1[1] = d * sin(-(phi_1 - double_Pi / 2)) + position[1];
	bow_p1[2] = -d * cos(-(phi_1 - double_Pi / 2))*cos(double_Pi / 2) + position[2];

	bow_p2[0] = d * cos(-(phi_1 - double_Pi / 2))*sin(double_Pi / 2) + position[0];
	bow_p2[1] = -d * sin(-(phi_1 - double_Pi / 2)) + position[1];
	bow_p2[2] = d * cos(-(phi_1 - double_Pi / 2))*cos(double_Pi / 2) + position[2];
	*/

	bow_p1[0] = -d * cos(-(angle_with_xz_plane))*sin(double_Pi / 2) + position[0];
	bow_p1[1] = d * sin(-(angle_with_xz_plane)) + position[1];
	bow_p1[2] = -d * cos(-(angle_with_xz_plane))*cos(double_Pi / 2) + position[2];

	bow_p2[0] = d * cos(-(angle_with_xz_plane))*sin(double_Pi / 2) + position[0];
	bow_p2[1] = -d * sin(-(angle_with_xz_plane)) + position[1];
	bow_p2[2] = d * cos(-(angle_with_xz_plane))*cos(double_Pi / 2) + position[2];


	std::vector<double> dist(gNumStrings);
	std::vector<double> d_alt(gNumStrings);
	std::vector<hduVector3Dd> e_CD(gNumStrings);
	std::vector<int> pointOnString(gNumStrings);
	std::vector<int> pointOnBow(gNumStrings);

	//const double stringStiffness = 0.12;
	//const double stringStiffness = 0.15;
	const double stringStiffness = 0.25;
	//const double stringStiffness = 0.3;

	hduVector3Dd f_tot;
	hduVecSet(f_tot, 0.0, 0.0, 0.0);


	std::vector<hduVector3Dd> direction_vec(gNumStrings);
	for (int iString = 0; iString < gNumStrings; iString++)
	{
		dist[iString] = distanceBetweenBowAndString_M1(string_p1[iString], string_p2[iString], bow_p1, bow_p2);

		//Logger::getCurrentLogger()->outputDebugString("dist: (" + String(dist[iString]) + ")");

		intersectionPoints points = distanceBetweenBowAndString_M2(string_p1[iString], string_p2[iString], bow_p1, bow_p2);

		d_alt[iString] = sqrt((points.xD - points.xC) * (points.xD - points.xC) + (points.yD - points.yC) * (points.yD - points.yC) + (points.zD - points.zC) * (points.zD - points.zC));

		//Logger::getCurrentLogger()->outputDebugString("d_alt: (" + String(d_alt) + ")");

		e_CD[iString] = intersect_points_dir(points);

		pointOnString[iString] = find_if_point_on_string(points, string_p1[iString], string_p2[iString]);
		pointOnBow[iString] = find_if_point_on_bow(points, bow_p1, bow_p2);

		f_vec[iString] = -e_CD[iString]; // give feedback dir


		hduVector3Dd e1;
		e1[0] = string_p2[iString][0] - string_p1[iString][0];
		e1[1] = string_p2[iString][1] - string_p1[iString][1];
		e1[2] = string_p2[iString][2] - string_p1[iString][2];

		direction_vec[iString] = get_fr_force_direction(e1, e_CD[iString]);
		gDirection_vec[iString][0] = direction_vec[iString][0];
		gDirection_vec[iString][1] = direction_vec[iString][1];
		gDirection_vec[iString][2] = direction_vec[iString][2];

		if ((pointOnString[iString] == 1) && (pointOnBow[iString] == 1) && (dist[iString] < 0))
		{
			//Logger::getCurrentLogger()->outputDebugString("pointOnString: (" + String(pointOnString[iString]) + ")");

			double forceMag = stringStiffness * d_alt[iString];
			hduVecScale(f_vec[iString], f_vec[iString], forceMag);

			gActiveString[iString] = 1;

			gDirection_vec_print[iString][0] = gDirection_vec[iString][0];
			gDirection_vec_print[iString][1] = gDirection_vec[iString][1];
			gDirection_vec_print[iString][2] = gDirection_vec[iString][2];

		}
		else
		{
			//Logger::getCurrentLogger()->outputDebugString("pointOnString: (" + String(0) + ")");
			double forceMag = 0;
			hduVecScale(f_vec[iString], f_vec[iString], forceMag);

			gActiveString[iString] = 0;
		}

		hduVecAdd(f_tot,
			f_tot,
			f_vec[iString]);
	}


	gF = f_tot;



	// PLANE STUFF
	// Stiffnes, i.e. k value, of the plane.  Higher stiffness results
	// in a harder surface.
	const double planeStiffness = 0.25;
	const double popthroughForceThreshold = 1000.0;

	static int directionFlag_y = 1;
	static int directionFlag_z = 1;
	static int directionFlag_z2 = -1;

	hduVector3Dd f_tot_plane;
	hduVector3Dd f_z;
	hduVector3Dd f_z2;


	if ((position[2] <= zMin && directionFlag_z > 0) || // remember to match these bounds with the length of the string ! 
		(position[2] > zMin) && (directionFlag_z < 0))
	{
		double penetrationDistance = fabs(position[2] - zMin);
		hduVector3Dd forceDirection(0, 0, directionFlag_z);

		double k = 1.0; // higher stiffness here

		hduVector3Dd z = penetrationDistance * forceDirection;
		hduVector3Dd f = k * z;

		if (f.magnitude() > popthroughForceThreshold)
		{
			f.set(0.0, 0.0, 0.0);
			directionFlag_z = -directionFlag_z;
		}

		f_z = f;
	}
	else
	{
		f_z.set(0.0, 0.0, 0.0);
	}

	if ((position[2] <= zMax && directionFlag_z2 > 0) ||
		(position[2] > zMax) && (directionFlag_z2 < 0))
	{
		//double penetrationDistance = fabs(position[2]);
		double penetrationDistance = fabs(zMax - position[2]);
		hduVector3Dd forceDirection(0, 0, directionFlag_z2);

		double k = 1.0; // higher stiffness here

		hduVector3Dd z = penetrationDistance * forceDirection;
		hduVector3Dd f = k * z;

		if (f.magnitude() > popthroughForceThreshold)
		{
			f.set(0.0, 0.0, 0.0);
			directionFlag_z2 = -directionFlag_z2;
		}

		f_z2 = f;
	}
	else
	{
		f_z2.set(0.0, 0.0, 0.0);
	}


	// maybe add another plane somewhere.. 

	hduVecAdd(f_z2, f_z2, f_z);

	hduVecAdd(f_tot_plane, f_tot_plane, f_z2);

	hduVecAdd(gF, gF, f_tot_plane);





	// VIBRATION STUFF


	hduVector3Dd force_vibr_tot;
	hduVecSet(force_vibr_tot, 0.0, 0.0, 0.0);



	for (int iString = 0; iString < gNumStrings; iString++)
	{
		hduVector3Dd force_vibr;

		// TO DO : maybe vary the factor for each string ! 
		//double fac = exponentialMapping(maxFB, 0.0, 0.2, 1.0, -3.0, FB_vec[iString]); // TO DO : NEEDS CALIBRATION WORK
		//double fac = exponentialMapping(maxFB, 0.0, 0.2, 0.6, -3.0, FB_vec[iString]); // TO DO : NEEDS CALIBRATION WORK
		//double fac = exponentialMapping(maxFB, 0.0, 0.2, 0.6, -3.0, FB_vec[iString]); // TO DO : NEEDS CALIBRATION WORK
		//double fac = exponentialMapping(maxFB, 0.0, 1.0, 2.0, -4.0, FB_vec[iString]); // TO DO : NEEDS CALIBRATION WORK
		//double fac = exponentialMapping(maxFB, 0.0, 0.8, 1.8, -4.0, FB_vec[iString]); // TO DO : NEEDS CALIBRATION WORK
		double fac = exponentialMapping(maxFB, 0.0, 1.4, 2.6, -4.0, FB_vec[iString]); // TO DO : NEEDS CALIBRATION WORK


		hduVecScale(force_vibr, direction_vec[iString], gF_fr[iString] * fac);

		hduVecAdd(force_vibr_tot, force_vibr_tot, force_vibr);

		if (gF_fr[iString] == 0)
		{
			gDirection_vec[iString][0] = 0;
			gDirection_vec[iString][1] = 0;
			gDirection_vec[iString][2] = 0;
		}

	}

	//Logger::getCurrentLogger()->outputDebugString("FB: (" + String(FB * fac) + ") ");
	//Logger::getCurrentLogger()->outputDebugString("gF_fr: (" + String(gF_fr * fac) + ") ");
	//Logger::getCurrentLogger()->outputDebugString("fx: (" + String(force_vibr[0]) + ") fy: (" + String(force_vibr[1]) + ") fz: (" + String(force_vibr[2]) + ")");
	//Logger::getCurrentLogger()->outputDebugString("x: (" + String(direction[0]) + ") y: (" + String(direction[1]) + ") z: (" + String(direction[2]) + ")");


	hduVecAdd(force_vibr_tot, force_vibr_tot, gF);


	/*
	if (hapticFrictionFlag)
	{
		hdSetDoublev(HD_CURRENT_FORCE, force_vibr_tot); // WITH FRICTION FEEDBACK
	}
	else if (hapticStringStiffnessFlag)
	{
		hdSetDoublev(HD_CURRENT_FORCE, gF); // NO FRICTION FEEDBACK
	}
	*/

	if (includeHaptics)
	{
		hdSetDoublev(HD_CURRENT_FORCE, force_vibr_tot); // WITH FRICTION FEEDBACK
		//Logger::getCurrentLogger()->outputDebugString("with friction!");
	}
	else
	{
		hdSetDoublev(HD_CURRENT_FORCE, gF); // NO FRICTION FEEDBACK
		//Logger::getCurrentLogger()->outputDebugString("NO friction!");
	}


	//hdSetDoublev(HD_CURRENT_FORCE, gF); // NO FRICTION FEEDBACK



	hdEndFrame(hdGetCurrentDevice());

	return HD_CALLBACK_CONTINUE;

}





HDCallbackCode HDCALLBACK DevicePositionCallback(void *pUserData)
{

	HDdouble *pPosition = (HDdouble *)pUserData;

	hdBeginFrame(hdGetCurrentDevice());
	hdGetDoublev(HD_CURRENT_POSITION, pPosition);
	hdEndFrame(hdGetCurrentDevice());

	return HD_CALLBACK_CONTINUE;
}




//==============================================================================
MainComponent::MainComponent() : minOut(-1.0), maxOut(1.0), numStrings(gNumStrings), octave(0), polyphony(12) // what is the : thing here ? numbered list?
{
	// Make sure you set the size of the component after
	// you add any child components.

	// specify the number of input and output channels that we want to open
	setAudioChannels(0, 2);

	// remove mouse cursor from this
	setMouseCursor(MouseCursor::NoCursor);

	// SLIDER STUFF:sa


	addAndMakeVisible(globalDampingSlider);
	globalDampingSlider.setSliderStyle(juce::Slider::Rotary);
	globalDampingSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0); // Add a textbox..
	globalDampingSlider.setTextValueSuffix(" [-]");
	globalDampingSlider.setRange(0.0, 1.0);          // maybe scale these values afterwards
	globalDampingSlider.setValue(0.5); // with exp mapping this will result in 1.0ish
	//dampingSlider.setValue(1.5);
	globalDampingSlider.addListener(this);

	globalDampingLabel.setText("Damping Amount", juce::dontSendNotification);
	globalDampingLabel.attachToComponent(&globalDampingSlider, false); // [4]
	addAndMakeVisible(globalDampingLabel);

	addAndMakeVisible(dampingSlider);
	dampingSlider.setSliderStyle(juce::Slider::Rotary);
	dampingSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0); // Add a textbox..
	dampingSlider.setTextValueSuffix(" [-]");
	dampingSlider.setRange(0.1, 10.0);          // maybe scale these values afterwards
	dampingSlider.setValue(4.5); // with exp mapping this will result in 1.0ish
	//dampingSlider.setValue(1.5);
	dampingSlider.addListener(this);

	dampingLabel.setText("Freq Indep Damping", juce::dontSendNotification);
	dampingLabel.attachToComponent(&dampingSlider, false); // [4]
	addAndMakeVisible(dampingLabel);


	addAndMakeVisible(freqDampingSlider);
	freqDampingSlider.setSliderStyle(juce::Slider::Rotary);
	freqDampingSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0); // Add a textbox..
	freqDampingSlider.setTextValueSuffix(" [-]");
	freqDampingSlider.setRange(0.0001, 0.015);          // maybe scale these values afterwards
	freqDampingSlider.setValue(0.0111); // ADD EXP MAPPING ! with exp mapping this will result in 10.0ish
	freqDampingSlider.addListener(this);

	freqDampingLabel.setText("Freq Dep Damping", juce::dontSendNotification);
	freqDampingLabel.attachToComponent(&freqDampingSlider, false); // [4]
	addAndMakeVisible(freqDampingLabel);


	addAndMakeVisible(frParamSlider);
	frParamSlider.setSliderStyle(juce::Slider::Rotary);
	frParamSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
	frParamSlider.setTextValueSuffix(" [-]");
	frParamSlider.setRange(0.01, 15000.0);
	//frParamSlider.setValue(7500.0); // with exp mapping this will result in 100.0ish
	frParamSlider.setValue(9700.0); // with exp mapping this will result in 100.0ish
	frParamSlider.addListener(this);

	frParamLabel.setText("Friction Amount", juce::dontSendNotification);
	frParamLabel.attachToComponent(&frParamSlider, false); // [4]
	addAndMakeVisible(frParamLabel);


	addAndMakeVisible(volumeSlider);
	volumeSlider.setSliderStyle(juce::Slider::Rotary);
	volumeSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
	volumeSlider.setTextValueSuffix(" [-]");
	volumeSlider.setRange(-10.0, 10.0); // in dB. conversion done in sliderValueChanged routine
	volumeSlider.setValue(-1.0);
	volumeSlider.addListener(this);

	volumeLabel.setText("Volume", juce::dontSendNotification);
	volumeLabel.attachToComponent(&volumeSlider, false);
	addAndMakeVisible(volumeLabel);


	// BUTTON

	addAndMakeVisible(hapticFrictionButton);
	hapticFrictionButton.setButtonText("Friction Button");
	hapticFrictionButton.addListener(this); // [7]
	//hapticFrictionButton.onClick = [this] { doButton(); };

	hapticFrictionLabel.setText("Friction ON", juce::dontSendNotification);
	hapticFrictionLabel.attachToComponent(&hapticFrictionButton, true); // [4]
	addAndMakeVisible(hapticFrictionLabel);

	/*
	addAndMakeVisible(hapticFrictionLabel);
	hapticFrictionLabel.setColour(juce::Label::backgroundColourId, juce::Colours::black);
	hapticFrictionLabel.setColour(juce::Label::textColourId, juce::Colours::white);
	hapticFrictionLabel.setJustificationType(juce::Justification::centred);
	*/

	// Make component opaque
	setOpaque(true);

}




MainComponent::~MainComponent()
{
	//hapticFrictionButton.removeListener(this);

	// This shuts down the audio device and clears the audio source.
	shutdownAudio();

	//// PHANTOM STUFF
	hdStopScheduler();
	hdDisableDevice(hHD);
}

void MainComponent::buttonClicked(juce::Button* button)
{
	if (button == &hapticFrictionButton)
	{
		//Logger::getCurrentLogger()->outputDebugString("Button Clicked");

		if (includeHaptics)
		{
			includeHaptics = false;
			hapticFrictionLabel.setText("Friction OFF", juce::dontSendNotification);

			//Logger::getCurrentLogger()->outputDebugString("NO friction!");
		}
		else
		{
			includeHaptics = true;
			hapticFrictionLabel.setText("Friction ON", juce::dontSendNotification);

			//Logger::getCurrentLogger()->outputDebugString("YES friction!");
		}
	}
}


void MainComponent::sliderValueChanged(Slider* slider)
{
	if (slider == &globalDampingSlider)
	{
		double dampAmount = globalDampingSlider.getValue();

		gSig0 = linearMapping(1.0, 0.0, 10.0, 0.1, dampAmount);
		gSig0 = exponentialMapping(10.0, 0.1, 10.0, 0.1, 4, gSig0);

		gSig1 = linearMapping(1.0, 0.0, 0.015, 0.0001, dampAmount);
		gSig1 = exponentialMapping(0.015, 0.0001, 0.015, 0.0001, 2, gSig1);

		Logger::getCurrentLogger()->outputDebugString("gSig0: (" + String(gSig0) + ")");
		Logger::getCurrentLogger()->outputDebugString("gSig1: (" + String(gSig1) + ")");
	}
	/*
	if (slider == &dampingSlider)
	{
		gSig0 = dampingSlider.getValue();
		gSig0 = exponentialMapping(10.0, 0.1, 10.0, 0.1, 4, gSig0);

		Logger::getCurrentLogger()->outputDebugString("gSig0: (" + String(gSig0) + ")");
	}
	else if (slider == &freqDampingSlider)
	{
		gSig1 = freqDampingSlider.getValue();
		gSig1 = exponentialMapping(0.015, 0.0001, 0.015, 0.0001, 2, gSig1);
		////        gFrParam = exponentialMapping(1500.0, 1.0, 1500.0, 1.0, 4.0, gFrParam);
		//        gFrParam = exponentialMapping(15000.0, 1.0, 15000.0, 1.0, 8.0, gFrParam);
		//        gSig1 = linearMapping(15000.0, 1.0, 2.0, 0.0, gFrParam);

		//        gSig1 = linearMapping(15000.0, 1.0, 1, 15000.0, gFrParam);
		//        gSig1 = exponentialMapping(1.0, 15000.0, 1.0, 15000.0, -8.0, gFrParam);

		Logger::getCurrentLogger()->outputDebugString("gSig1: (" + String(gSig1) + ")");
	}
	*/
	else if (slider == &frParamSlider)
	{
		gFrParam = frParamSlider.getValue();
		////        gFrParam = exponentialMapping(1500.0, 1.0, 1500.0, 1.0, 4.0, gFrParam);
		//        gFrParam = exponentialMapping(15000.0, 1.0, 15000.0, 1.0, 8.0, gFrParam);


		gStickFact = linearMapping(15000.0, 1.0, 2.0, 0.0, gFrParam);

		gFrParam = linearMapping(15000.0, 0.01, 0.01, 15000.0, gFrParam);
		gFrParam = exponentialMapping(0.01, 15000.0, 0.01, 15000.0, -8.0, gFrParam);

		/*
		gStickFact = linearMapping(15000.0, 1.0, 2.0, 0.0, gFrParam);

		gFrParam = linearMapping(15000.0, 1.0, 1, 15000.0, gFrParam);
		gFrParam = exponentialMapping(1.0, 15000.0, 1.0, 15000.0, -8.0, gFrParam);
		*/
		Logger::getCurrentLogger()->outputDebugString("gFrParam: (" + String(gFrParam) + ")");
	}
	else if (slider == &volumeSlider)
	{
		gVolume = volumeSlider.getValue();
		gVolume = powf(10.0, gVolume / 20);
		if (gVolume < 0.317)
		{
			gVolume = 0;
		}

		Logger::getCurrentLogger()->outputDebugString("gVolume: (" + String(gVolume) + ")");

	}
}


//==============================================================================
void MainComponent::prepareToPlay(int samplesPerBlockExpected, double sampleRate)
{

	hiResCallbackFreq = 150.0;
	HighResolutionTimer::startTimer(1000.0 / hiResCallbackFreq); // 150 Hz

	// Other stuff
	globalCurrentSample = 0;

	fs = sampleRate;
	bufferSize = samplesPerBlockExpected;

	std::vector<double> freq_strings_violin{ 196.00 , 293.66 , 440.00 , 659.25 };

	int test = 0; // why use this..
	for (int i = test; i < numStrings + test; ++i)
	{
		//ViolinStrings.add(new ViolinString(440, fs));
		//ViolinStrings.add(new ViolinString(330 + 75 * i, fs));
		//ViolinStrings.add(new ViolinString(220 + 110 * i, fs));
		//ViolinStrings.add(new ViolinString(220 + 55 * i, fs));
		//ViolinStrings.add(new ViolinString(75 + 75*i, fs));

		ViolinStrings.add(new ViolinString(freq_strings_violin[i], fs));


		ViolinStrings[i]->setFb(0.0);
		ViolinStrings[i]->setVb(0.0);
		addAndMakeVisible(ViolinStrings[i]);
	}

	FB_Max = 2.0;
	FB_Max = 1.0;
	opa_level = 0.9;

	//maxVb = 0.4;
	//maxFB = 2.0;

	preOutputVec.resize(numStrings, 0);
	gOutputStrings.resize(numStrings, 0);
	keyDownVec.resize(numStrings, false);
	bowOnVec.resize(numStrings, false);


	activeStrings.resize(polyphony, nullptr);

	// Timer for graphics
	startTimerHz(15);


	gPosition.resize(3, 0);
	gPositionPrev.resize(3, 0);


	gPositionRaw.resize(2, 0);
	gPositionScaled.resize(2, 0);

	// position of bow ?
	x_inp_var = 0.3; // bowing pos in percentage
	y_inp_var = 0.5;

	setSize((int)widthScreen, (int)heightScreen);


	//I_B.resize(ViolinStrings[0]->getN() + 1,0);
	//J_B.resize(ViolinStrings[0]->getN() + 1,0);


	// PHANTOM STUFF:

	// Initialize the default haptic device.
	hHD = hdInitDevice(HD_DEFAULT_DEVICE);

	// Start the servo scheduler and enable forces.
	hdEnable(HD_FORCE_OUTPUT);
	//hdSetSchedulerRate(500); // needs to be set before u start scheduler ! 
	hdSetSchedulerRate(1000); // needs to be set before u start scheduler ! 
	//hdSetSchedulerRate(2000); // needs to be set before u start scheduler ! 
	hdStartScheduler();

	// Pass by reference global value gPosition to continuously update ! 
	hdScheduleSynchronous(DevicePositionCallback, &gPositionCallback, HD_DEFAULT_SCHEDULER_PRIORITY); // this might be redundant.. 
	//hdScheduleAsynchronous(StringCallback, 0, HD_MAX_SCHEDULER_PRIORITY); // WHY IS IT WITH ASYNCHRONOUS ?? have vibrationCallback as last in chain. keep as only one with max scheduler priority !! i.e. revert to default scheduler priority for prev force stuff.. 
	hdScheduleSynchronous(StringCallback, 0, HD_MAX_SCHEDULER_PRIORITY); //  have vibrationCallback as last in chain. keep as only one with max scheduler priority !! i.e. revert to default scheduler priority for prev force stuff.. 

	//hdScheduleAsynchronous(FrictionlessPlaneCallback, 0, HD_DEFAULT_SCHEDULER_PRIORITY); // for forces max scheduler priority is needed..
	//hdScheduleAsynchronous(VibrationCallback, 0, HD_MAX_SCHEDULER_PRIORITY); // have vibrationCallback as last in chain. keep as only one with max scheduler priority !! i.e. revert to default scheduler priority for prev force stuff.. 

}

void MainComponent::getNextAudioBlock(const juce::AudioSourceChannelInfo& bufferToFill)
{
	if (soundFlag)
	{


		// updating global model parameters every buffersize.. for efficiency. may cause some clicks
		for (int iString = 0; iString < numStrings; ++iString)
		{
			/*
			ViolinStrings[iString]->setSig0(1.0);
			ViolinStrings[iString]->setSig1(0.008);
			ViolinStrings[iString]->setFrParam(80);
			//ViolinStrings[iString]->setFrParam(2000);
			ViolinStrings[iString]->setStickFact(1);
			*/

			ViolinStrings[iString]->setSig0(gSig0);
			ViolinStrings[iString]->setSig1(gSig1);
			ViolinStrings[iString]->setFrParam(gFrParam);
			//ViolinStrings[iString]->setFrParam(2000);
			ViolinStrings[iString]->setStickFact(gStickFact);

		}

		// update volume at every buffersize to avoid clicks.. (maybe)
		double volPerBuffer = gVolume;
		//Logger::getCurrentLogger()->outputDebugString("output: (" + String(gOutputStrings[0]) + ")");
		//Logger::getCurrentLogger()->outputDebugString("gF_fr: (" + String(gF_fr) + ")");

		// Your audio-processing code goes here!
		for (int channel = 0; channel < bufferToFill.buffer->getNumChannels(); ++channel)
		{
			float *const channelData = bufferToFill.buffer->getWritePointer(channel, bufferToFill.startSample);

			if (channel == 0)
			{
				// These should be multidim arrays.. one for each string
				// resetting the interpolation and spreading operators for new bowing locations
				for (int j = 0; j < numStrings; ++j)
				{
					std::vector<double> I_B = ViolinStrings[j]->getI_B();
					std::vector<double> J_B = ViolinStrings[j]->getJ_B();
					N = ViolinStrings[j]->getN();
					double h = ViolinStrings[j]->get_h();

					bp = x_inp_var; // new bowing location !
					bP = floor(bp * N - 1); // -1 for alignment with Matlab..
					alpha_bow = bp * N - 1 - bP;

					int bP_m1 = bP - 1;
					int bP_p1 = bP + 1;
					int bP_p2 = bP + 2;

					std::fill(I_B.begin(), I_B.end(), 0); // fill with zeros

					I_B[bP_m1] = (alpha_bow * (alpha_bow - 1) * (alpha_bow - 2)) / -6.0;
					I_B[bP] = ((alpha_bow - 1) * (alpha_bow + 1) * (alpha_bow - 2)) / 2.0;
					I_B[bP_p1] = (alpha_bow * (alpha_bow + 1) * (alpha_bow - 2)) / -2.0;
					I_B[bP_p2] = (alpha_bow * (alpha_bow + 1) * (alpha_bow - 1)) / 6.0;

					for (int i = 0; i < I_B.size(); ++i)
					{
						J_B[i] = I_B[i] * (1 / h); // speed up: keep divisions out of loop !
					}

					ViolinStrings[j]->setI_B(I_B);
					ViolinStrings[j]->setJ_B(J_B);
					//                Logger::getCurrentLogger()->outputDebugString("FB: (" + String(ViolinStrings[j]->getFb()) + ")");
					//                Logger::getCurrentLogger()->outputDebugString("vB: (" + String(ViolinStrings[j]->getVb()) + ")");
				}
				for (int i = 0; i < bufferToFill.buffer->getNumSamples(); i++)
				{

					float output = 0.0;
					for (int j = 0; j < numStrings; ++j)
					{
						//                    float gainVal = exponentialMapping((float) numMasses - 1.0, 0.0, 550.0, 120.0, -1.5, (float) j); // this is a reasonable mapping..
						float gainVal = 750;
						float stringSound = ViolinStrings[j]->process() * gainVal;

						gOutputStrings[j] = stringSound;
						output = output + stringSound;

						if (ViolinStrings[j]->isBowOn())
						{
							bowOnVec[j] = true;
							gOutput = stringSound; // output for vibration feedback in pen
							gF_fr[j] = ViolinStrings[j]->getF_fr(); // output for vibration feedback in pen

							//Logger::getCurrentLogger()->outputDebugString("F_fr_test: (" + String(mass_spring_dampers[j]->getF_fr()) + ")");
						}
						else
						{
							bowOnVec[j] = false;
							gF_fr[j] = 0; // output for vibration feedback in pen
						}
					}
					int sum = std::accumulate(std::begin(bowOnVec), std::end(bowOnVec), 0); // making sure there is no vibration feedback in pen if none of the masses have the bow on
					if (sum == 0)
					{
						gOutput = 0;
						for (int j = 0; j < numStrings; ++j)
						{
							gF_fr[j] = 0;
						}
					}

					//                Logger::getCurrentLogger()->outputDebugString("output: (" + String(output) + ")");

					output = output * volPerBuffer; // add global volume contribution

					if (output > maxOut) // output is limited..
					{
						Logger::getCurrentLogger()->outputDebugString("Output is too large!");
						output = maxOut;
					}
					else if (output < minOut) {
						Logger::getCurrentLogger()->outputDebugString("Output is too small!");
						output = minOut;
					}
					channelData[i] = output;
				}
			}
			else
			{
				memcpy(channelData,
					bufferToFill.buffer->getReadPointer(0),
					bufferToFill.numSamples * sizeof(float));
			}
		}
	}
}

void MainComponent::releaseResources()
{
}

//==============================================================================
void MainComponent::paint(juce::Graphics& g)
{
	// clear the background
	g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));

	float brushStroke_crossSec = 0.03f*heightScreen;
	float brushStroke = 0.009f*heightScreen;

	//float brushStroke_crossSec = 0.03*heightScreen;
	//float brushStroke = 0.009*heightScreen;

	float lineThickness = 0.025f*heightScreen; // for bow thickness

	// DRAW STRINGS - PATHS (LINES)
	int count = 0;
	//std::vector<std::vector<double>> uVecsTest;
	std::vector<double> strMax_vals(gNumStrings, 0);
	std::vector<double> signMax_vals(gNumStrings, 0);
	//for (ViolinString* str : ViolinStrings)
	for (auto str : ViolinStrings)
	{
		if (str != nullptr)
		{
			// choose your favourite colour
			g.setColour(Colours::cyan);

			//double startPos_x = linearMapping(210, -210, widthScreen, widthScreen / 2, string_p1[count][0]);
			//float startPos_x = static_cast<float>(linearMapping(210.0 * 0.7, -210.0 * 0.7, static_cast<double>(getWidth()), static_cast<double>(getWidth()) / 2.0, string_p1[count][0]));
			float startPos_x = (float)linearMapping(210.0 * 0.7, -210.0 * 0.7, (double)widthScreen, (double)widthScreen / 2.0, string_p1[count][0]);

			//float startPos_x_scr2 = static_cast<float>(linearMapping(210.0 * 0.7, -210.0 * 0.7, static_cast<double>(getWidth()) / 2.0, 0.0, string_p1[count][0]));
			//float startPos_y_scr2 = static_cast<float>(linearMapping(-110.0 * 0.7, 205.0 * 0.7, static_cast<double>(getHeight()) / 2.0, 0.0, string_p1[count][1]));

			// draw the state
			g.strokePath(visualiseState(g, 50000.0f, str, startPos_x), PathStrokeType(brushStroke));
			//g.strokePath(visualiseState_crossSec(g, 80000, str, startPos_x_scr2, startPos_y_scr2, gDirection_vec[count]), PathStrokeType(brushStroke_crossSec));



			for (int l = 0; l <= static_cast<int>(str->getN()); l++) // if you don't save the boundaries use l < N
			{
				//strMax_vals[count] = max(fabs(str->get_uVecs()[1][l]), strMax_vals[count]);
				if (fabs(str->get_uVecs()[1][l]) > strMax_vals[count])
				{
					strMax_vals[count] = fabs(str->get_uVecs()[1][l]);
					signMax_vals[count] = sgn(str->get_uVecs()[1][l]);
				}


			}
			count++;
			//            int ana = 3;
		}
	}

	// I think the 0.7 scaling fucks up some stuff.. think of something else !

	// DRAW STRINGS - ELLIPSES 

	float spacing_x = static_cast<float>(heightScreen / 4.0);
	float spacing_y = static_cast<float>(heightScreen / 4.0);

	for (int iString = 0; iString < numStrings; iString++)
	{
		float startPos_x = static_cast<float>(linearMapping(210 * 0.7, -210 * 0.7, (double)widthScreen / 2.0, 0, string_p1[iString][0]));
		float startPos_y = static_cast<float>(linearMapping(-110 * 0.7, 205 * 0.7, (double)heightScreen / 2.0, 0, string_p1[iString][1]));

		g.setColour(Colours::red);
		g.fillEllipse(startPos_x - brushStroke_crossSec / 2 + strMax_vals[iString] * signMax_vals[iString] * 90000 * gDirection_vec_print[iString][0], // gDirection_vec is zero when the bow is not on the string.. need to keep it as the gDir of the last bow position (in the case of free vibration)
			startPos_y - brushStroke_crossSec / 2 - strMax_vals[iString] * signMax_vals[iString] * 90000 * gDirection_vec_print[iString][1],
			brushStroke_crossSec,
			brushStroke_crossSec);

		g.setColour(Colours::cyan);
		g.fillEllipse(startPos_x - brushStroke_crossSec / 2, startPos_y - brushStroke_crossSec / 2, brushStroke_crossSec, brushStroke_crossSec);
	}


	// DRAW FIRST BOW (VIEWED FROM ABOVE)  WORKS ! YEY 

	float xStart = static_cast<float>(linearMapping(210.0 * 0.7, -210.0 * 0.7, (double)widthScreen, (double)widthScreen / 2.0, bow_p1[0]));
	float yStart = static_cast<float>(linearMapping(zMax, zMin, heightScreen - 0.075*heightScreen, 0 + 0.075*heightScreen, bow_p1[2]));
	float xEnd = static_cast<float>(linearMapping(210.0 * 0.7, -210.0 * 0.7, (double)widthScreen, (double)widthScreen / 2.0, bow_p2[0]));
	float yEnd = static_cast<float>(linearMapping(zMax, zMin, (double)heightScreen - 0.075*heightScreen, 0 + 0.075*heightScreen, bow_p2[2]));

	//float xEnd = gBowEnd_2[0];
	//float yEnd = gBowEnd_2[1];

	juce::Line<float> line(juce::Point<float>(xStart, yStart),
		juce::Point<float>(xEnd, yEnd));

	if (flagMouseUp)
	{
		g.setColour(Colours::grey);
		g.setOpacity(opa_level);
	}
	else
	{
		g.setColour(Colours::orange);
		g.setOpacity(opa_level);
	}

	g.drawLine(line, lineThickness); // line thickness!

	// DRAW SECOND BOW (CROSS SECTION) (double)

	xStart = static_cast<float>(linearMapping(210.0 * 0.7, -210.0 * 0.7, (double)widthScreen / 2.0, 0.0, bow_p1[0]));
	yStart = static_cast<float>(linearMapping(-110.0 * 0.7, 205.0 * 0.7, (double)heightScreen / 2.0, 0.0, bow_p1[1]));
	yStart = yStart - lineThickness;
	xEnd = static_cast<float>(linearMapping(210.0 * 0.7, -210.0 * 0.7, (double)widthScreen / 2.0, 0.0, bow_p2[0]));
	yEnd = static_cast<float>(linearMapping(-110.0 * 0.7, 205.0 * 0.7, (double)heightScreen / 2.0, 0.0, bow_p2[1]));
	yStart = yStart - lineThickness;

	juce::Line<float> line2(juce::Point<float>(xStart, yStart),
		juce::Point<float>(xEnd, yEnd));

	if (flagMouseUp)
	{
		g.setColour(Colours::grey);
		g.setOpacity(opa_level);
	}
	else
	{
		g.setColour(Colours::orange);
		g.setOpacity(opa_level);
	}

	g.drawLine(line2, lineThickness);


	//g.fillRect(xRect,yRect,widthRect,heightRect);
//    g.fillAll(Colours::orange);
}


Path MainComponent::visualiseState(Graphics& g, float visualScaling, ViolinString* string, float startPos_x)
{

	float spacing_x = getWidth() / 2 / (numStrings + 1);

	float startPos_y = 0;

	// initialise path
	Path stringPath;

	// start path
	//stringPath.startNewSubPath (stringBoundaries_y,u[1][0] * visualScaling);
	stringPath.startNewSubPath(startPos_x, startPos_y);

	float spacing_y = getHeight() / static_cast<float>(string->getN());
	float y = 0;

	for (int l = 0; l <= static_cast<int>(string->getN()); l++) // go through all the points ! (N+1 points, N intervals)
	{

		// Needs to be -u, because a positive u would visually go down
		float newX = startPos_x + static_cast<float>(string->get_uVecs()[1][l]) * visualScaling;

		// if we get NAN values, make sure that we don't get an exception
		if (isnan(newX))
			newX = 0;
		if (isnan(y))
			y = 0;

		stringPath.lineTo(newX, y);
		y += spacing_y;
	}
	return stringPath;
}



Path MainComponent::visualiseState_crossSec(Graphics& g, double visualScaling, ViolinString* string, double startPos_x, double startPos_y, std::vector<double> direction)
{

	/*
	double spacing_y = getHeight() / 4;
	double spacing_x = getHeight() / 4;

	double startPos_y = getHeight() / 4 + spacing_y*strNo;
	double startPos_x = getWidth() / 4 + spacing_x*strNo;
	*/

	// initialise path
	Path stringPath;

	// start path
	stringPath.startNewSubPath(startPos_x, startPos_y);

	for (int l = 0; l <= static_cast<int>(string->getN()); l++) // if you don't save the boundaries use l < N
	{
		// here u need to split the mag into the vibration directions (make global var, update in StringCallback)
		float newX = static_cast<float>(startPos_x + string->get_uVecs()[1][l] * visualScaling * direction[0]);
		float newY = static_cast<float>(startPos_y - string->get_uVecs()[1][l] * visualScaling * direction[1]); // not sure why i need minus sign here..
		//float newY = static_cast<float>(startPos_y + string->get_uVecs()[1][l] * visualScaling * direction[1]);

		// if we get NAN values, make sure that we don't get an exception
		if (isnan(newY))
			newY = 0;

		if (isnan(newX))
			newX = 0;

		stringPath.lineTo(newX, newY);
	}
	// if you don't save the boundaries, and add a stringPath.lineTo (x, getWidth()) here to end the statedrawing
	return stringPath;
}







void MainComponent::paintOverChildren(juce::Graphics& g) // IMPORTANT !
{
}

void MainComponent::resized()
{
	heightScreen = getHeight();
	widthScreen = getWidth();

	/*
	int x_dampSlider = widthScreen / 16;
	int y_dampSlider = heightScreen / 2 + heightScreen / 5;
	int width_dampSlider = widthScreen / 8;
	int height_dampSlider = widthScreen / 8;

	dampingSlider.setBounds(x_dampSlider, y_dampSlider, width_dampSlider, height_dampSlider);


	int x_freqDampingSlider = x_dampSlider + widthScreen / 8;
	int y_freqDampingSlider = heightScreen / 2 + heightScreen / 5;
	int width_freqDampingSlider = widthScreen / 8;
	int height_freqDampingSlider = widthScreen / 8;

	freqDampingSlider.setBounds(x_freqDampingSlider, y_freqDampingSlider, width_freqDampingSlider, height_freqDampingSlider);


	*/

	int x_dampSlider = widthScreen / 16;
	int y_dampSlider = heightScreen / 2 + heightScreen / 5;
	int width_dampSlider = widthScreen / 8;
	int height_dampSlider = widthScreen / 8;

	globalDampingSlider.setBounds(x_dampSlider, y_dampSlider, width_dampSlider, height_dampSlider);


	int x_frParamSlider = x_dampSlider + widthScreen / 8;
	int y_frParamSlider = heightScreen / 2 + heightScreen / 5;
	int width_frParamSlider = widthScreen / 8;
	int height_frParamSlider = widthScreen / 8;

	frParamSlider.setBounds(x_frParamSlider, y_frParamSlider, width_frParamSlider, height_frParamSlider);

	int x_volumeSlider = x_frParamSlider + widthScreen / 8;
	int y_volumeSlider = heightScreen / 2 + heightScreen / 5;
	int width_volumeSlider = widthScreen / 8;
	int height_volumeSlider = widthScreen / 8;

	volumeSlider.setBounds(x_volumeSlider, y_volumeSlider, width_volumeSlider, height_volumeSlider);


	int x_button = widthScreen / 8;
	int y_button = heightScreen / 2 + heightScreen / 4 + heightScreen / 8;
	int width_button = widthScreen / 4;
	int height_button = widthScreen / 20;

	hapticFrictionButton.setBounds(x_button, y_button, width_button, height_button);
}


void MainComponent::timerCallback()
{
	repaint(); // update the graphics X times a second

	// yaw = left - right, pitch = up - down, roll = rotate stylus
	//Logger::getCurrentLogger()->outputDebugString("gBowEnd_1[0]: (" + String(gBowEnd_1[0]) + ") gBowEnd_1[1]: (" + String(gBowEnd_1[1]) + ") gBowEnd_1[2]: (" + String(gBowEnd_1[2]) + ")");
	//Logger::getCurrentLogger()->outputDebugString("gBowEnd_2[0]: (" + String(gBowEnd_2[0]) + ") gBowEnd_2[1]: (" + String(gBowEnd_2[1]) + ") gBowEnd_2[2]: (" + String(gBowEnd_2[2]) + ")");

	//Logger::getCurrentLogger()->outputDebugString("x: (" + String(gPosition[0]) + ") y: (" + String(gPosition[1]) + ") z: (" + String(gPosition[2]) + ")");
	//Logger::getCurrentLogger()->outputDebugString("xRaw: (" + String(gPositionRaw[0]) + ") yRaw: (" + String(gPositionRaw[1]) + ")");

}



void MainComponent::hiResTimerCallback()
{
	// positions mapped to pixels.. 
	gPosition[0] = gPositionCallback[0];
	gPosition[1] = gPositionCallback[1];
	gPosition[2] = gPositionCallback[2];

	//Logger::getCurrentLogger()->outputDebugString("gPosition[0]: (" + String(gPosition[0]) + ") gPosition[1]: (" + String(gPosition[1]) + ") gPosition[2]: (" + String(gPosition[2]) + ")");
	//Logger::getCurrentLogger()->outputDebugString("gPositionPrev[0]: (" + String(gPositionPrev[0]) + ") gPositionPrev[1]: (" + String(gPosition[1]) + ") gPositionPrev[2]: (" + String(gPosition[2]) + ")");


	gPositionRaw[0] = linearMapping(185.0, -185.0, widthScreen, 0.0, gPosition[0]);
	//gPositionRaw[1] = linearMapping(zMax, zMin, heightScreen - 0.02*heightScreen, 0.0 + 0.02*heightScreen, gPosition[2]);
	gPositionRaw[1] = linearMapping(zMax, zMin, heightScreen - 0.04*heightScreen, 0.0 + 0.04*heightScreen, gPosition[2]);

	// positions remapped between 0 and 1
	gPositionScaled[0] = linearMapping(widthScreen, 0.0, 1.0, 0.0, gPositionRaw[0]);
	gPositionScaled[1] = linearMapping(heightScreen - 0.075*heightScreen, 0.0 + 0.075*heightScreen, 1.0, 0.0, gPositionRaw[1]);

	// 3d positions of bow end points based on gimbal angles and position of gimbal center.
	// THERE IS SOMETHING OFF WITH THIS FUCKIN MAPPING ! FIX IT


	gBowEnd_1[0] = linearMapping(185.0, -185.0, widthScreen, 0.0, bow_p1[0]);
	gBowEnd_2[0] = linearMapping(185.0, -185.0, widthScreen, 0.0, bow_p2[0]);

	gBowEnd_1[1] = linearMapping(zMax, zMin, heightScreen - 0.02*heightScreen, 0.0 + 0.02*heightScreen, bow_p1[2]);
	gBowEnd_2[1] = linearMapping(zMax, zMin, heightScreen - 0.02*heightScreen, 0.0 + 0.02*heightScreen, bow_p2[2]);


	x_inp_var = gPositionScaled[1];
	x_inp_var = clamp(x_inp_var, 0.075, 0.925);


	//if (gPosition[1] < -0.1)
	if (std::accumulate(gActiveString.begin(), gActiveString.end(), 0) != 0)
	{
		flagBow = true;
		flagMouseUp = false;

		//Logger::getCurrentLogger()->outputDebugString("Bowing!");

	}
	//else if (gPosition[1] >= -0.1)
	else
	{
		flagBow = false;
		flagMouseUp = true;
	}

	if (flagBow)
	{
		for (int j = 0; j < numStrings; ++j)
		{
			if (gActiveString[j] == 1)
			{
				ViolinStrings[j]->setBow(true);
				ViolinStrings[j]->activate();
			}
			else
			{
				ViolinStrings[j]->setBow(false);
			}
		}
	}

	//double VbRaw = (gPosition[0] - gPositionPrev[0]) * hiResCallbackFreq;
	double VbRaw = (gPosition[0] - gPositionPrev[0]) * hiResCallbackFreq; // [pixels/s]
	double VbRawAlt = sqrt((gPosition[0] - gPositionPrev[0])*(gPosition[0] - gPositionPrev[0]) + (gPosition[1] - gPositionPrev[1])*(gPosition[1] - gPositionPrev[1])*(gPosition[2] - gPositionPrev[2])) * hiResCallbackFreq;

	if (isnan(VbRawAlt)) // I get nans here sometimes.. don't know why
		VbRawAlt = 0;

	//Logger::getCurrentLogger()->outputDebugString("VbRaw: (" + String(VbRaw) + ")");
	//Logger::getCurrentLogger()->outputDebugString("VbRawAlt: (" + String(VbRawAlt) + ")");


	//double maxVb = 0.4;
	//double Vb = linearMapping(-12.0, 12.0, -maxVb, maxVb, VbRaw);
	double Vb = linearMapping(-180.0, 180.0, -maxVb, maxVb, VbRawAlt);
	//if (Vb > -0.03 && Vb < 0.03) // limit Vb to avoid small noise during small displacements
	if (Vb > -0.05 && Vb < 0.05) // limit Vb to avoid small noise during small displacements
	{
		Vb = 0.0;
	}
	if (Vb > 0.4)
	{
		Vb = 0.4;
	}

	//Logger::getCurrentLogger()->outputDebugString("Vb: (" + String(Vb) + ")");
					//Logger::getCurrentLogger()->outputDebugString("FB: (" + String(mass_spring_dampers[idx]->getFb()) + ")");

	for (int j = 0; j < numStrings; ++j) // update bowing vel on all masses at once (better!)
	{
		ViolinStrings[j]->setVb(Vb);
	}


	for (int j = 0; j < numStrings; ++j) // update bowing vel on all masses at once (better!)
	{
		FB_vec[j] = linearMapping(10.0, 0.0, maxFB, 0.0, f_vec[j].magnitude());
		FB_vec[j] = clamp(FB_vec[j], 0.0, maxFB);
	}




	//FB = linearMapping(-0.1, -200.0, 0.0, maxFB, gPosition[1]);
	FB = linearMapping(10.0, 0.0, maxFB, 0.0, gF.magnitude());
	FB = clamp(FB, 0.0, maxFB);

	//Logger::getCurrentLogger()->outputDebugString("FB: (" + String(ViolinStrings[0]->getFb()) + ")");

	opa_level = FB / maxFB;
	if (opa_level < 0.1)
	{
		opa_level = 0.1;
	}
	else if (opa_level > 0.95)
	{
		opa_level = 0.95;
	}

	// make sure nothing happens when you don't bow..
	if (Vb == 0) {
		for (int j = 0; j < numStrings; ++j) // update bowing vel on all masses at once (better!)
		{
			FB_vec[j] = 0;
		}
	}

	// only update force now..
	for (int j = 0; j < numStrings; ++j) // update force on all masses at once (better!)
	{
		//ViolinStrings[j]->setFb(FB);
		if (gActiveString[j] == 1)
		{
			ViolinStrings[j]->setFb(FB_vec[j]);
		}
		else
		{
			ViolinStrings[j]->setFb(0.0);
		}
	}


	if (!flagBow)
	{
		//gPositionPrev[0] = gPosition[0];
		//gPositionPrev[1] = gPosition[1];
		//gPositionPrev[2] = gPosition[2];

		FB = 0;
		Vb = 0;
		for (int j = 0; j < numStrings; ++j)
		{
			ViolinStrings[j]->setBow(false);
			ViolinStrings[j]->setVb(Vb);
			ViolinStrings[j]->setFb(FB);
		}

	}

	gPositionPrev[0] = gPosition[0];
	gPositionPrev[1] = gPosition[1];
	gPositionPrev[2] = gPosition[2];
	//memcpy(&gPositionPrev, &gPosition, sizeof(gPosition)); // is this fast.. ?

}







