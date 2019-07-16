// ThorlabsImagerOCT.cpp : Defines the exported functions for the DLL application - OCT Part

#include "stdafx.h"
#include "ThorlabsImager.h"
#include <SpectralRadar.h>
#include <string>
#include <iostream>
#include <windows.h>
#include "lasercontrol.h"
using namespace std;

//Utility
std::wstring s2ws(const std::string& s)
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	std::wstring r(buf);
	delete[] buf;
	return r;
}

//DLL Internal state varibles:
static OCTDeviceHandle Dev_;
static ProbeHandle Probe_;

//OCT Probe specific info
static string octProbeName;
static double octScanSpeed; //A scans per sec

//Initialize OCT Scanner
void yOCTScannerInit(
	const char probeFilePath[]  // Probe ini path. Can be usually found at: C:\\Program Files\\Thorlabs\\SpectralRadar\\Config\\ 
	)
{
	Dev_ = initDevice();
	Probe_ = initProbe(Dev_,probeFilePath);

	//TBD do this a probe specific
	octProbeName = "Ganymede";
	octScanSpeed = 28000;
}

//Close OCT Scanner, Cleanup
void yOCTScannerClose()
{
	closeProbe(Probe_);
	closeDevice(Dev_);
}

// Scan a 3D Volume
void yOCTScan3DVolume(
	const double xStart, //Scan start position [mm]
	const double yStart, //Scan start position [mm]
	const double rangeX, //fast direction length [mm]
	const double rangeY, //slow direction length [mm]
	const double rotationAngle, //Scan angle [deg]
	const int    sizeX,  //Number of pixels on the fast direction
	const int    sizeY,  //Number of pixels on the slow direction
	const int    nBScanAvg, //Number of B scan averages (set to 1 if non)
	const char   outputDirectory[] //Output folder, make sure it exists and empty before running this function
)
{
	string outputDirectoryStr = outputDirectory;

	// Create Data Handles
	RawDataHandle Raw = createRawData();
	DataHandle BScan = createData();
	DataHandle Volume = createData();
	char message[1024];

	// Setup Scan and Start Scanning
	ScanPatternHandle Pattern = createBScanStackPattern(Probe_, rangeX, sizeX, rangeY, sizeY);
	rotateScanPattern(Pattern, rotationAngle);
	startMeasurement(Dev_, Pattern, Acquisition_AsyncFinite);

	// Looping over bscans
	for (int i = 0; i < sizeY; ++i)
	{
		for (int j = 0; j < nBScanAvg; ++j)
		{
			//Get the raw data
			getRawData(Dev_, Raw);

			//Generate file name to export data to
			char filename[1024];
			sprintf_s(filename, (outputDirectoryStr + "\\Data_Y%04d_YTotal%d_B%04d_BTotal%d_" + octProbeName + ".srr").c_str(), i + 1, sizeY, j + 1, nBScanAvg);

			//Export data
			exportRawData(Raw, RawDataExport_SRR, filename);
			if (getError(message, 512))
			{
				cerr << "\n\nAn error occurred: " << message << endl;
			}
		}
	}

	//Copy Chirp File
	string sourceChirp = "C:\\Program Files\\Thorlabs\\SpectralRadar\\Config\\Chirp.dat";
	CopyFile((s2ws(sourceChirp)).c_str(),(s2ws(outputDirectoryStr)).c_str(),FALSE);

	//Cleanup
	stopMeasurement(Dev_);
	setLaserDiode(Dev_, FALSE);
	clearRawData(Raw);
	clearScanPattern(Pattern);
}


// Photobleach a Line
void yOCTPhotobleachLine(
	const double xStart,	//Start position [mm]
	const double yStart,	//Start position [mm]
	const double xEnd,		//End position [mm]
	const double yEnd,		//End position [mm]
	const double duration,	//How many seconds to photobleach
	const double repetition //How many times galvoes should go over the line to photobleach. slower is better. recomendation: 1
	)
{
	double dwellPerPass = duration / repetition; //How many seconds to dwell per one pass of the laser beam on the B scan
	int howManyAScansPerBScan = (int)(octScanSpeed*dwellPerPass); //Total number of A Scans asqquired
	
	//Setup
	ScanPatternHandle Pattern = createBScanPatternManual(Probe_, xStart, yStart, xEnd, yEnd, howManyAScansPerBScan, FALSE);

	//Turn laser diode on
	controllaser(TRUE);

	//Photobleach
	startMeasurement(Dev_, Pattern, AcquisitionType::Acquisition_AsyncContinuous);
	Sleep((long)(1000 * duration)); //Sleep in msec
	stopMeasurement(Dev_);

	//Turn laser diode off
	controllaser(FALSE);

	//Cleanup
	clearScanPattern(Pattern);
}