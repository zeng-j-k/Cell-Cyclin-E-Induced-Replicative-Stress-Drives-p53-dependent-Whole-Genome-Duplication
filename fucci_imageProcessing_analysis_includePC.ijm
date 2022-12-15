// Matt Renshaw, 2020, CALM-STP, The Francis Crick Institute. Jingkun Zeng, The Francis Crick Institute
// Analysis of FUCCI cells
// H2B-mTurquoise, Geminin-mVenus, Cdt-mCherry, phase contrast
// H2B used to track cells using Trackmate
// Acquisition order of the channels should be H2B-mTurquoise, Geminin-mVenus, Cdt-mCherry, phase contrast

// -------------------------------------------------------------------------------------------
// --- user defined variables ---

tInterval = 1200; // time interval (seconds)
channelNames = newArray("h2b","geminin","cdt"); // acquisition order of channels; can be changed, but names must be preserved
roi = 20; // diameter of circular area (pixels) for measurements

// -------------------------------------------------------------------------------------------
// --- IMPORTANT! works on 4-channel timelapse image that is already open ---

run("Clear Results");
roiManager("reset");
title = getTitle();


getDimensions(width, height, channels, slices, frames);
if (channels != 4 || frames == 1) {
	exit("4 channel timelapse dataset required");
}


// -------------------------------------------------------------------------------------------
// --- initial image processing ---

// --- split channels and rename

selectWindow(title);
run("Split Channels");
for (i = 0; i < 3; i++) {
	selectWindow("C"+toString(i+1)+"-"+title);
	rename(channelNames[i]);
}

// --- image processing ---
for (i = 0; i < 3; i++) {	
	setBatchMode(true);
	
	// background subtraction using Minimum and Gaussian filters
	selectWindow(channelNames[i]);
	run("Subtract Background...", "rolling=50 stack");
	
//	run("Duplicate...", "title=min duplicate");

//	selectWindow("min");
//	run("Minimum...", "radius=30 stack");
//	run("Gaussian Blur...", "sigma=50 stack");
//	imageCalculator("Subtract stack", channelNames[i], "min");
//	closeWindow("min");

	// median filter
	selectWindow(channelNames[i]);
	run("Median...", "radius=5 stack");
	
	setBatchMode(false);	
}


// -------------------------------------------------------------------------------------------
// --- trackmate ---

selectWindow("h2b");
resetMinAndMax();
run("Enhance Contrast", "saturated=0.35");

run("TrackMate");

// --- TrackMate settings ---
text = "enter TrackMate settings and press Analysis to generate tracks \n \n";
text = text + "Select detector - Downsample LoG detector \n";
text = text + "Estimated blob diameter: 15.0 microns \n";
text = text + "Threshold: 2 \n";
text = text + "Downsampling factor: 4 \n";
text = text + "Initial thresholding: no thresholding \n";
text = text + "Select a view: Hyperstack Displayer \n";
text = text + "Set filters on spots: None \n";
text = text + "Select a tracker: Nearest neighbor search \n";
text = text + "Maximal linking distance: 15.0 micron \n";
text = text + "Set filters on tracks: Number of spots in track: Above 71.x \n";
text = text + "Display options: Press analysis \n";

waitForUser(text);

// --------------------------------------------------------------------------------------------
// Using TrackMate tracks measure h2b, geminin and cdt for each object

if (isOpen("Track statistics")==0) {
	waitForUser("no tracks found, press Analysis to generate tracks \n" +text);
}

// --- get all tracks ---
selectWindow("Track statistics");
trackID = Table.getColumn("TRACK_ID");
nTracks = trackID.length;
nSpots = Table.getColumn("NUMBER_SPOTS");
trackStart = Table.getColumn("TRACK_START");
trackStop = Table.getColumn("TRACK_STOP");

selectWindow("Spots in tracks statistics");
xPositions = Table.getColumn("POSITION_X");
yPositions = Table.getColumn("POSITION_Y");
tPositions = Table.getColumn("POSITION_T");

selectWindow("h2b");
//run("Hide Overlay");
print(nTracks);
Array.print(trackID);

// --- measure individual tracks ---
closeWindow("Results");
start = 0;
for (i = 0; i < nTracks; i++) {
	
	currentTrack = trackID[i];
	end = start + nSpots[i];
	print(currentTrack, start, end);

	// get X, Y and T coordinates for current track
	xCoords = Array.slice(xPositions,start,end);
	yCoords = Array.slice(yPositions,start,end);
	tCoords = Array.slice(tPositions,start,end);
	 for (t = 0; t < tCoords.length; t++) {
	 	tCoords[t] = tCoords[t]/tInterval;
	 }

	 // create ROIs for current track
	 selectWindow("h2b");
	 run("Select None");
	 roiManager("reset");
	 
	 getPixelSize(unit, pixelWidth, pixelHeight);
	 for (j = 0; j < xCoords.length; j++) {
	 	setSlice(tCoords[j]+1);
	 	run("Select None");
	 	makeOval((xCoords[j]/pixelWidth)-roi/2, (yCoords[j]/pixelHeight)-roi/2, roi, roi);
	 	roiManager("Add");
	 }

	 // measure and normalise mean intensity for all 3 channels, current track
	 for (c = 0; c < 3; c++) {
	 	setBatchMode(true);
	 	
	 	run("Select None");
	 	run("Clear Results");
	 	run("Set Measurements...", "mean stack display redirect=None decimal=3");
	 	selectWindow(channelNames[c]);
	 	count = roiManager("count");
	 	
	 	// measure
	 	for (n = 0; n < count; n++) {
	 		roiManager("select", n);
	 		run("Measure");
	 		//setResult("Track_ID", nResults-1, currentTrack);
	 	}
	 	
	 	// normalise mean intensity between min and max
	 	Table.rename("Results", "channel-Results");
	 	selectWindow("channel-Results");
	 	resultsMean = Table.getColumn("Mean");
	 	resultsTime = Table.getColumn("Slice");
	 	Array.getStatistics(resultsMean, min, max, mean, stdDev);
	 	if (isOpen("temp")) {
	 		Table.rename("temp", "Results");
	 	}
	 	
	 	for (j = 0; j < (resultsMean.length); j++) {
	 		setResult("Label", j, title);
	 		setResult("TRACK_ID", j, currentTrack);
	 		setResult("Timepoint", j, resultsTime[j]);
	 		setResult(channelNames[c]+"_Mean_RAW", j, resultsMean[j]);
	 		setResult(channelNames[c]+"_Mean_NORM", j, (resultsMean[j]-min)/(max-min));
	 	}
	 	closeWindow("channel-Results");
	 	Table.rename("Results", "temp");
	 	setBatchMode(false);
	 	
	 }

	 // --- check if cells express fucci construct at sufficient level
	 selectWindow("temp");
	 timepoints = Table.getColumn("Timepoint");
	 h2bRAW = Table.getColumn("h2b_Mean_RAW");
	 h2bNORM = Table.getColumn("h2b_Mean_NORM");
	 gemininRAW = Table.getColumn("geminin_Mean_RAW");
	 gemininNORM = Table.getColumn("geminin_Mean_NORM");
	 cdtRAW = Table.getColumn("cdt_Mean_RAW");
	 cdtNORM = Table.getColumn("cdt_Mean_NORM");
	 Array.getStatistics(cdtRAW, minCdt, maxCdt, meanCdt, stdDevCdt);
	 Array.getStatistics(gemininRAW, minGem, maxGem, meanGem, stdDevGem);
	 	 
	 if (isOpen("output")) {
	 		Table.rename("output", "Results");
	 	}
	 	row = nResults;
	
	 for (x = 0; x < timepoints.length; x++) {
	 		setResult("Label", row+x, title);
	 		setResult("TRACK_ID", row+x, currentTrack);
	 		setResult("Timepoint", row+x, timepoints[x]);
	 		setResult("geminin_Mean_RAW", row+x, gemininRAW[x]);
	 		setResult("geminin_Mean_NORM", row+x, gemininNORM[x]);
	 		setResult("cdt_Mean_RAW", row+x, cdtRAW[x]);
	 		setResult("cdt_Mean_NORM", row+x, cdtNORM[x]);
	 		setResult("h2b_Mean_RAW", row+x, h2bRAW[x]);
	 		setResult("h2b_Mean_NORM", row+x, h2bNORM[x]);	 		
	 	}

	 Table.rename("Results", "output");
	 closeWindow("temp");
	 start = end;
}

Table.rename("output", "Results");
for (i = 0; i < nResults; i++) {
	cdt = getResult("cdt_Mean_NORM", i);
	gem = getResult("geminin_Mean_NORM", i);
	//setResult("cdt/gem_RATIO", i, cdt/gem);
	setResult("fucci_TOTAL", i, (cdt+gem));
	setResult("cdt_PROP", i, (cdt/(cdt+gem)));
	setResult("geminin_PROP", i, (gem/(cdt+gem)));
}

selectWindow("Results");
saveAs("results");

// --------------------------------------------------------------------------------------------
// --- functions ---
// close any window without returning any error
function closeWindow(name) {
	if (isOpen(name)) {
	     selectWindow(name);
	     run("Close");
	}
}