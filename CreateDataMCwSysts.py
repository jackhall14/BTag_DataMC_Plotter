import ROOT
from ROOT import gDirectory, gStyle, gPad, TColor, gROOT
import os
import argparse
import re
import json
from math import sqrt, log10, floor, pow
import logging
logging.basicConfig(level=logging.INFO)

# 
# Choose between the methods for creating the input file for the plotter (--CreateFile)
# or the plotter method once that file has been created (--PlotFile)
# 

def main():
	args = get_args()

	# TDirectory channels you want to run over
	TDirNames = ["emu_OS_J2"]

	logging.info("Loading jsons ...")

	Plots = byteify(json.load(file(args.Plots)))
	logging.info("Loaded plots to create")

	Samples = byteify(json.load(file(args.Samples)))
	logging.info("Loaded samples included in plots")

	NomSamples = Samples["NomSamples"]
	SystSamples = Samples["SystSamples"]
	DataSamples = Samples["DataSamples"]

	EventVar = False
	for Type, TypeDict in Plots.items():
		logging.info("Running code over " + Type + " ... ")
		if Type == "Events":
                        EventVar = True
                else:
                        EventVar = False

		Objs = TypeDict["Objs"]
		Vars = TypeDict["Vars"]

		logging.debug("Looping over:\t" + str(Objs) + "\t and:\t" + str(Vars))
		if args.CreateFile:
			CreateSystFile(args, TDirNames, Objs, Vars, NomSamples, SystSamples, DataSamples, EventVar)
		elif args.PlotFile:
			FilePlotter(args, TDirNames, Objs, Vars, NomSamples, SystSamples, DataSamples, EventVar)
	logging.info("Finished")

def FilePlotter(Args, TDirectoryNames, Objects, Variables, NominalSamples, SystSamples, DataSamples, EventVar):
	logging.info("---------------------------------")
	logging.info("Beginning file plotting")
	logging.info("---------------------------------")
	logging.info("Working with filename: " + str(Args.PlotFilename))

	if not os.path.exists(Args.OutputDir):
		# Begin just by creating the output path:
		os.makedirs(Args.OutputDir)

	PlotFile = tfile(Args.PlotFilename)
	logging.debug("Successfully opened plotting file")

	# What we want to do is create a linux dir per TDir and create a plot per obj and var
	for TDirectoryName in TDirectoryNames:
		logging.info("Creating histograms for the TDirectoryName: " + TDirectoryName)

		cwd = Args.OutputDir + TDirectoryName + "/"
		if not os.path.exists(cwd):
			os.makedirs(cwd)

		for Obj in Objects:
			for Var in Variables:
				logging.info("------------------------------------------------------------")
				logging.info("      Gathering info for plot:\t\t" + str(Obj) + " " + str(Var))
				logging.info("------------------------------------------------------------")

				if EventVar:
					NameCheck = Obj
				else:
					NameCheck = Obj+"_"+Var

				NominalHists = GetNominalContributions(PlotFile, TDirectoryName, NameCheck, NominalSamples)
				logging.debug(NominalHists)

				DataHists = GetDataContributions(PlotFile, TDirectoryName, NameCheck, DataSamples)
				logging.debug(DataHists)

				SystBand = CreateSystematicBand(PlotFile, TDirectoryName, NameCheck, NominalSamples+SystSamples)
				logging.debug(SystBand)

				if not Args.NoBatch:
					ExportPlot(TDirectoryName, NominalHists, DataHists, SystBand, True)
				else:
					ExportPlot(TDirectoryName, NominalHists, DataHists, SystBand)
	PlotFile.Close()

def ExportPlot(TDirName, NominalHists, DataHists, SystBand, BatchMode = False):
	if BatchMode:
		ROOT.gROOT.SetBatch()

	c = ROOT.TCanvas("c2", "",800,700)
	c.cd()
	gStyle.SetOptStat(0)
	gStyle.SetHatchesLineWidth(1)
	gStyle.SetErrorX(0)

	logging.info("Defining colour palette ... ")
	colours = {}
	colours["Singletop"]= TColor(3000,54./255, 121./255,191./255)
	colours["ttbar"]=       TColor(3001,123./255, 178./255, 116./255)
	colours["Wjets"]= TColor(3002,130./255, 95./255, 135./255)
	colours["Zjets"]= TColor(3003,252./255, 176./255, 8./255)
	colours["Diboson"]=TColor(3007,168./255, 164./255, 150./255)

	# # Begin work on top pad
	pad1 = CreateTopPad("pad1")
	pad1.cd()
	pad1.SetLogy()
	SetTicks(pad1)

	# Create TLegend
	Legend = CreateLegend()

	# Get systematic band but dont draw
	TopUncertBand = DrawUncHist(SystBand, colours, Legend)

	# Draw MC
	SMHist, MCStack = DrawSMHists(NominalHists, colours, Legend)
	# Then draw the MC stack
	MCStack.Draw("HIST L")
	logging.info("Drawn MC stack")

	SMHist = ApplySystematicBand(SMHist, TopUncertBand)
	# This gives us the error on the SM Hist
	SMHist.Draw("E3 SAME")
	# This gives us the line as well as the hashed fill
	# SMHist.Draw("C SAME")
	logging.info("Drawn SM total with systematic band")

	# Draw Data
	DataHist = DrawDataHists(DataHists, colours, Legend)
	DataHist.Draw("E0 SAME")
	logging.info("Drawn data")

	# Draw Legend
	logging.info("Drawing legend ... ")
	Legend.Draw()

	logging.info("Adding more aesthetic changes ...")
	# Draw ATLAS Text
	ATLASText = DrawATLAS()

	# Begin work on DataMC Ratio
	c.cd()
	pad2 = CreateBottomPad("pad2")
	SetTicks(pad2)
	
	DataMCHist, UncBand = DrawDataMCPad(DataHist, SMHist, pad2)
	DataMCHist.Draw("ep")
	UncBand.Draw("E3 SAME")
	logging.info("Drawn dataMC ratio")

	c.Update()
	logging.info("Finished this plot!")

	PlotName = "Final_"+DataHist.GetName().split("_nominal")[0]
	c.SaveAs("Plots/" + TDirName + "/" + PlotName+".pdf")
	if not BatchMode:
		raw_input()

def DrawUncHist(UncBand, Colours, TLegend, TopPad=True):
	logging.info("Drawing systematic band ... ")

	# Draw line
	UncBand.SetLineWidth(3)
	UncBand.SetLineStyle(0)
	UncBand.SetLineColor(4)
	return UncBand

def SortByIntegral(Hist):
	return Hist.Integral()

def DrawSMHists(HistList, Colours, TLegend, SystematicBand=False):
	logging.info("Drawing MC ... ")

	Iter = 0
	MCHists = []
	for Hist in HistList:
		logging.debug(str(Hist) + " " + str(type(Hist)))
		if "ttbar" in Hist.GetName():
			SampleName = Hist.GetName().split("_")[-3]
		else:
			SampleName = Hist.GetName().split("_")[-2]
		logging.info("...... histogram: " + str(SampleName))

		# # For each SM sample
		Hist.SetFillStyle(1001)
		Hist.SetFillColor(Colours[SampleName].GetNumber())
		Hist.SetLineColor(Colours[SampleName].GetNumber())
		Hist.SetLineWidth(0)

		if Iter == 0:
			# Create a combined SM hist
			SMHist = Hist.Clone()
		else:
			SMHist.Add(Hist.Clone())

		TLegend.AddEntry(Hist,SampleName)
		# Add to MC Stack
		MCHists.append(Hist)
		Iter = Iter + 1
	
	# Sort by integral
	MCHists.sort(key=SortByIntegral)

	# Add to MCStack
	MCStack = ROOT.THStack("MCStack","")
	for Hist in MCHists:
		MCStack.Add(Hist)

	YMax = GetYAxisRange(SMHist)
	logging.debug("Y Max is " + str(YMax))
	MCStack.SetMaximum(YMax*1000.0)
	MCStack.SetMinimum(0.1)

	logging.debug("Drawing full sum of SM contributions ... ")

	# First label axes
	SMHist.GetYaxis().SetTitle("No. of Events")
	# Adjust y labelling
	SMHist.GetYaxis().SetTitleSize(20)
	SMHist.GetYaxis().SetTitleFont(43)
	SMHist.GetYaxis().SetTitleOffset(1.55)

	# Remove title
	SMHist.SetNameTitle("", "")

	# Draw combined SM sample
	SMHist.SetLineWidth(2)
	SMHist.SetLineStyle(0)
	SMHist.SetLineColor(28)
	SMHist.SetFillStyle(3345)
	SMHist.SetFillColor(1)
	TLegend.AddEntry(SMHist,"SM Total")
	return SMHist, MCStack

def DrawDataHists(HistList, Colours, TLegend):
	logging.info("Drawing data ... ")

	Iter = 0
	for Hist in HistList:
		logging.debug(str(Hist) + " " + str(type(Hist)))
		SampleName = Hist.GetName().split("_")[-2]
		logging.info("...... Drawing hist: " + str(SampleName))
		
		if Iter == 0:
			# Create a combined data hist
			ExistHist = Hist.Clone()
		else:
			ExistHist.Add(Hist.Clone())
		Iter = Iter + 1

	# Draw line
	ExistHist.SetMarkerStyle(ROOT.kFullCircle)
	ExistHist.SetLineWidth(2)
	ExistHist.SetLineStyle(0)
	ExistHist.SetLineColor(1)
	TLegend.AddEntry(ExistHist, "Data")
	return ExistHist

def RenameVarObj(HistName):
	Var = HistName.split("_")[5]
	Obj = HistName.split("_")[4]

	if Obj == "jet1":
		xObj = "j_{1}"
	elif Obj == "jet2":
		xObj = "j_{2}"
	elif Obj == "el1":
		xObj = "e^{-}_{1}"
	elif Obj == "mu1":
		xObj = "#mu_{1}"
	else:
		xObj = Obj
	xlabel = xObj + " - "

	if Var == "pt" or Var == "pt1":
		xVar = "p_{T}"
	elif Var == "eta" or Var == "eta1":
		xVar = "#eta"
	elif Var == "cl_eta1":
		xVar = "Calo. #eta"
	elif Var == "phi" or Var == "phi1":
		xVar = "#phi"
	elif Var == "m":
		xVar = "m_{lj}"
	else:
		xVar = Var
	xlabel = xlabel + xVar

	logging.debug(xlabel)
	return xlabel

def DrawDataMCPad(DataHist, SMHist, TPad):
	logging.info("Drawing dataMC pad ... ")
	TPad.cd()

	# Get the ratio
	DataMCHist = DataHist.Clone()
	BottomUncBand = SMHist.Clone()
	OutputHistUncBand2 = SMHist.Clone()

	# Divide by itself as this will be the bottom unc band
	BottomUncBand.Divide(BottomUncBand)

	for Bin in range(0,OutputHistUncBand2.GetNcells()):
		# Set the bin error to zero
		OutputHistUncBand2.SetBinError(Bin, 0.0)

	DataMCHist.Divide(OutputHistUncBand2)

	DataMCHist.SetMinimum(0.5)
	DataMCHist.SetMaximum(1.5)

	# Gunna get labels for x axis
	xlabel = RenameVarObj(DataHist.GetName())
	DataMCHist.GetXaxis().SetTitle(xlabel)

	DataMCHist.SetNameTitle("", "")

	# X axis ratio plot settings
	DataMCHist.GetXaxis().SetTitleSize(20)
	DataMCHist.GetXaxis().SetTitleFont(43)
	DataMCHist.GetXaxis().SetTitleOffset(4.)
	DataMCHist.GetXaxis().SetLabelFont(43)
	DataMCHist.GetXaxis().SetLabelSize(15)

	# Y axis ratio plot settings
	DataMCHist.GetYaxis().SetTitle("Data/MC")
	DataMCHist.GetYaxis().SetNdivisions(505)
	DataMCHist.GetYaxis().SetTitleSize(20)
	DataMCHist.GetYaxis().SetTitleFont(43)
	DataMCHist.GetYaxis().SetTitleOffset(1.55)
	DataMCHist.GetYaxis().SetLabelFont(43)
	DataMCHist.GetYaxis().SetLabelSize(15)
	return DataMCHist, BottomUncBand

def CreateTopPad(Name):
	Name = ROOT.TPad(Name,"distribution",0.01,0.30,0.99,0.99,0)
	Name.SetTopMargin(0.01)
	Name.SetBottomMargin(0.01)
	Name.Draw()
	return Name

def CreateBottomPad(Name):
	Name = ROOT.TPad(Name,"dataMC",0.01,0.01,0.99,0.29,0)
	Name.SetTopMargin(0.04)
	Name.SetBottomMargin(0.25)
	Name.SetGridy()
	Name.Draw()
	return Name

def SetTicks(TPad):
	tx = TPad.GetTickx()
	ty = TPad.GetTicky()
	tx = 1
	ty = 1
	TPad.SetTicks(tx,ty)

def DrawATLAS():
	latexObject = ROOT.TLatex()
	latexObject.SetTextFont(42)
	latexObject.SetTextAlign(11)
	latexObject.SetTextColor(1)

	# ATLAS Text                                                                                                                           
	latexObject.SetTextSize(0.064)
	latexObject.DrawLatexNDC(0.15, 0.90, "#scale[1.2]{#bf{#it{ATLAS}} Internal}")

	# CoM / Lumi Text                                                                                                                      
	latexObject.SetTextSize(0.060)
	comEnergy = 13
	luminosity = 139100
	latexObject.DrawLatexNDC(0.15, 0.90 - 0.06, str(comEnergy)+" TeV "+str(luminosity/1.e3) + " fb^{-1}")
	return latexObject

def CreateLegend():
	legend = ROOT.TLegend(0.54,0.72,0.77,0.96)                                                                              
	legend.SetBorderSize(0)                                                                                                 
	legend.SetTextSize(0.04)
	legend.SetNColumns(2)
	legend.SetColumnSeparation(0.25)
	return legend

def GetYAxisRange(Histogram):
	LargestBinValue = 0
	for Bin in range(0,Histogram.GetNcells()):
		if Histogram.GetBinContent(Bin) > LargestBinValue:
			LargestBinValue = Histogram.GetBinContent(Bin)
	return round_to_1(LargestBinValue)

def round_to_1(x):
	return round(x, -int(floor(log10(abs(x)))))

def ApplySystematicBand(Hist, SystematicBand):
	logging.debug("Setting systematic band bin values as errors ... ")
	for Bin in range(0,Hist.GetNcells()):
		logging.debug("------------------------------------------------------------")
		logging.debug("Previous bin error:\t" + str(Hist.GetBinError(Bin)))
		logging.debug("Bin content:\t\t" + str(Hist.GetBinContent(Bin)))
		logging.debug("New uncert value:\t" + str(SystematicBand.GetPointY(Bin)))

		new_uncert = SystematicBand.GetPointY(Bin)
		if Hist.GetBinContent(Bin) != 0.0:
			mc_stat_unc = Hist.GetBinError(Bin)/Hist.GetBinContent(Bin)
			rel_error = new_uncert / Hist.GetBinContent(Bin)
		else:
			mc_stat_unc = 0.0
			rel_error = 0.0

		logging.debug("MC Stat uncert:\t" + str(mc_stat_unc))

		new_rel_error = sqrt(mc_stat_unc*mc_stat_unc + rel_error*rel_error)
		new_error = Hist.GetBinContent(Bin)*new_rel_error
		logging.debug("New bin error:\t" + str(new_error))

		Hist.SetBinError(Bin, new_error)
	return Hist

def GetBinEventYield(RootFile, TDirectory, ObjVar, Sample, Bin, Nominal=True):
	if Nominal:
		NomHist = RootFile.Get(TDirectory + "/"+Sample+"/"+"h_"+TDirectory+"_"+ObjVar+"_nominal")
		NomGr = ROOT.TGraphAsymmErrors(NomHist)
		EventYield = NomGr.GetPointY(Bin)
	return EventYield

def CreateSystematicBand(RootFile, TDirectoryName, ObjVar, AllSamples):
	# Need to loop over each TDirectoryname and take its systematic band
	logging.info("Getting the combined uncert TGraph now ... ")
	Iter = 0
	for Sample in range(0, len(AllSamples)):
		logging.debug("Taking the total uncertainty TGraph from " + AllSamples[Sample])
		TDir = RootFile.Get(TDirectoryName).GetDirectory(AllSamples[Sample])
		# Check if the dir is empty
		if len(TDir.GetListOfKeys()) == 0:
			logging.debug("Directory empty, skimming uncertainty")
			continue

		SystGr = TDir.Get("Gr_" + TDirectoryName + "_" + ObjVar + "_tot_uncert")

		# Check if it returns a the graph
		if SystGr == None:
			logging.debug("Problem returning uncertainty graph")
			continue

		if Iter == 0:
			for i in range(0, SystGr.GetN()):
				# Set the values to sqaured value
				# logging.debug("Uncertainty in " + str(i)+"th bin is " + str(SystGr.GetPointY(i)))
				BinEventYield = GetBinEventYield(RootFile, TDirectoryName, ObjVar, AllSamples[Sample], i)
				SystGr.SetPointY(i, pow(SystGr.GetPointY(i)*BinEventYield,2))
				FinalBand = SystGr
		else:
			for i in range(0, SystGr.GetN()):
				# Add in quadrature
				# logging.debug("Uncertainty in " + str(i)+"th bin is " + str(SystGr.GetPointY(i)))
				BinEventYield = GetBinEventYield(RootFile, TDirectoryName, ObjVar, AllSamples[Sample], i)
				FinalBand.SetPointY(i, FinalBand.GetPointY(i)+pow(SystGr.GetPointY(i)*BinEventYield,2))
		Iter = Iter + 1

	logging.debug("Sqaure-rooting the summed TGraphs ... ")
	for i in range(0, FinalBand.GetN()):
		# Add in quadrature
		FinalBand.SetPointY(i, sqrt(FinalBand.GetPointY(i)))
	logging.info("Finished getting the combined uncert band!")
	return FinalBand

def GetDataContributions(RootFile, TDirectoryName, ObjVar, DataSamples):
	logging.info("Getting all the data histograms now ... ")
	DataHists = []
	for Sample in DataSamples:
		hist_name = TDirectoryName + "/" + Sample + "/h_" + TDirectoryName + "_" + ObjVar + "_nominal"
		logging.debug(hist_name)
		DataSampleHist = RootFile.Get(hist_name)
		logging.debug(str(DataSampleHist) + str(type(DataSampleHist)))
		if DataSampleHist != None:
			logging.info("Got data histogram from sample:\t" + str(Sample))
			DataHists.append(DataSampleHist.Clone(DataSampleHist.GetName() + "_" + Sample))
	logging.info("Finished getting all data contributions!")
	return DataHists

def GetNominalContributions(RootFile, TDirectoryName, ObjVar, NominalSamples):
	logging.info("Getting all the nominal histograms now ... ")
	NomHists = []
	for Sample in NominalSamples:
		hist_name = TDirectoryName + "/" + Sample + "/h_" + TDirectoryName + "_" + ObjVar + "_nominal"
		NomSampleHist = RootFile.Get(hist_name)
		logging.debug(str(NomSampleHist) + str(type(NomSampleHist)))
		if NomSampleHist != None:
			logging.info("Got nominal histogram from sample:\t" + str(Sample))
			NomHists.append(NomSampleHist.Clone(NomSampleHist.GetName() + "_" + ("_").join(Sample.split("_")[1:])))
	logging.info("Finsihed getting all nominal contributions!")
	return NomHists

def CreateSystFile(Args, TDirectoryNames, Objects, Variables, NominalSamples, SystSamples, DataSamples, EventVar):
	logging.info("---------------------------------------------------")
	logging.info("Beginning file creation needed for plotting")
	logging.info("---------------------------------------------------")
	logging.info('Input path is: ' + str(Args.InputPath))

	# Open up the output file
	OutFile = tfile(Args.OutputFile,"UPDATE")
	logging.info("Outputting everything to: " + str(Args.OutputFile))

	AllInputDirectories = NominalSamples + SystSamples + DataSamples

	AllInputDirsClone = list(AllInputDirectories)
	for TDirectoryName in TDirectoryNames:
		# Using a certain TDirectory, calculate all the systematics
		logging.info("Calculating systematic uncertainties for the directory " + TDirectoryName)
		OutFile.mkdir(TDirectoryName)
		for Sample in AllInputDirectories:
			# Take a samples directory, get the nominal file from that directory and compare it for each systematic
			logging.info("Working in directory: " + str(Sample))
			OutFile.cd(TDirectoryName)
			gDirectory.mkdir(Sample)
			base = Args.InputPath + Sample + "/"

			# Get the nominal file
			if "data" in Sample:
				nominal_file = Sample+"_data_combination.root"
			elif Sample == "FTAG2_ttbar_PhPy8_hdamp3mtop":
				nominal_file = Sample+"_weight_mc_rad_UP_combination.root"
			else:
				nominal_file = Sample+"_nominal_combination.root"
			NomFile = tfile(base+nominal_file)
			if NomFile == None: continue
			logging.debug("Successfully opened nominal file: " + nominal_file)

			# Get a list of the systematic files
			if not "data" in Sample and Sample != "FTAG2_ttbar_PhPy8_hdamp3mtop":
				SystFiles = [f for f in os.listdir(base) if "FTAG2_" in f]
				SystFiles.remove(nominal_file)
				# if Sample == "FTAG2_ttbar_PhPy8":
				# 	bad_file = "FTAG2_ttbar_PhPy8_weight_mc_shower_np_131_combination.root"
				# 	SystFiles.remove(bad_file)
			else:
				SystFiles = []

			CompletedMethod = CalculateSystematics(NomFile, SystFiles, TDirectoryName, base, OutFile, Sample, Objects, Variables, EventVar)
			if not CompletedMethod:
				logging.info("Removing sample with missing histograms "+str(Sample))
				AllInputDirsClone.remove(Sample)

			NomFile.Close()
		OutFile.cd(TDirectoryName)
		OutFile.Write('', ROOT.TObject.kOverwrite)

	logging.info("-----------------------------------------")
	logging.info("Now calculating each total uncertainty ...")
	logging.info("-----------------------------------------")
	for TDirectoryName in TDirectoryNames:
		logging.info("Creating syst band for the TDirectoryName: " + TDirectoryName)
		for Sample in AllInputDirsClone:
			if "data" in Sample: continue
			for Obj in Objects:
				logging.info("And object " + Obj)
				for Var in Variables:
					logging.info("and variable " + Var)
					if Sample in SystSamples:
						logging.debug("Calculating for a systematic sample and not tree based syst ... ")
						CalculateSampleUncertainty(OutFile, TDirectoryName, Obj, Var, Sample, EventVar, True, Args.InputPath)
					else:
						logging.debug("Calcating a combined tree based systematic uncertainty")
						CalculateSampleUncertainty(OutFile, TDirectoryName, Obj, Var, Sample, EventVar, False)
	OutFile.Close()

def CalculateSampleUncertainty(RootFile, TDirectoryName, Object, Variable, Sample, EventVar, CalculateSystSamp=False, FilePath = False):
	# First create total uncertainty for nom vs tree systematics
	# To do this we will take qudrature sum of all uncerts
	if EventVar:
		NameCheck = Object
	else:
		NameCheck = Object+"_"+Variable
	logging.debug("NameCheck is:\t" + NameCheck)

	logging.debug("Going to begin creating the total uncertainty")
	if CalculateSystSamp == False:
		logging.info("Getting qudrature sum of all tree uncertainties for the sample: \t" + str(Sample))
		CurrentTDir = RootFile.Get(TDirectoryName + "/" + Sample)
		Iter = 0
		for Key in CurrentTDir.GetListOfKeys():
			if (NameCheck) in Key.GetName() and "nominal" not in Key.GetName():
				logging.debug("Found graph " + str(Key.GetName()))

				SystGr = CurrentTDir.Get(Key.GetName()).Clone()
				if Iter == 0:
					for i in range(0, SystGr.GetN()):
						# Set the values to sqaured value
						SystGr.SetPointY(i, SystGr.GetPointY(i)*SystGr.GetPointY(i))
						UncGr = SystGr
				else:
					for i in range(0, SystGr.GetN()):
						# Add in quadrature
						UncGr.SetPointY(i, UncGr.GetPointY(i)+(SystGr.GetPointY(i)*SystGr.GetPointY(i)))
				Iter = Iter + 1

		logging.debug("Sqaure-rooting the summed histograms ... ")
		for i in range(0, UncGr.GetN()):
			# Add in quadrature
			if UncGr.GetPointY(i) < 0.0:
				new_val = sqrt(UncGr.GetPointY(i) * -1.0)
			else:
				new_val = sqrt(UncGr.GetPointY(i))
			UncGr.SetPointY(i, new_val)

		# Change the name for outputting
		GrName = ("_").join(UncGr.GetName().split("_")[:3]+[NameCheck,"tot_uncert"])
		UncGr.SetTitle(GrName)
		UncGr.SetName(GrName)
		logging.info("Final combined graph: " + str(UncGr.GetName()))

		logging.info("Writing to output file ... ")
		RootFile.cd(TDirectoryName + "/" + Sample)
		UncGr.Write('', ROOT.TObject.kOverwrite)
	else:
		logging.info("Calculating uncertainty based upon which syst sample it is " + Sample)
		# if "FTAG2_ttbar_PhPy8_hdamp3mtop" in Sample:
		# 	logging.info("Comparing to hdamp3mtop sample ... ")
		# 	NomGr = RootFile.Get(TDirectoryName + "/"+Sample+"/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
		# 	SystGr = RootFile.Get(TDirectoryName + "/"+Sample+"/"+"h_"+TDirectoryName+"_"+NameCheck+"_weight_mc_rad_UP")
		# 	UncHist = SystematicSampleWrapper(NomGr, SystGr)
		if "FTAG2_ttbar_Sherpa221" in Sample:
			logging.info("Comparing to nominal AF2 ttbar sample ... ")
			NomGr = RootFile.Get(TDirectoryName + "/FTAG2_ttbar_PhPy8/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			SystGr = RootFile.Get(TDirectoryName + "/"+Sample+"/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			UncHist = SystematicSampleWrapper(NomGr, SystGr)
		elif "ttbar" in Sample:
			logging.info("Comparing to nominal AF2 ttbar sample ... ")
			NomGr = RootFile.Get(TDirectoryName + "/FTAG2_ttbar_PhPy8_AF2/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			if NomGr == None: return None
			SystGr = RootFile.Get(TDirectoryName + "/"+Sample+"/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			UncHist = SystematicSampleWrapper(NomGr, SystGr)
		elif "Singletop" in Sample:
			logging.info("Comparing to nominal (AF2) singletop sample ... ")
			NomGr = RootFile.Get(TDirectoryName + "/FTAG2_Singletop_PowPy8_AF2/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			SystGr = RootFile.Get(TDirectoryName + "/"+Sample+"/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			UncHist = SystematicSampleWrapper(NomGr, SystGr)
		elif "Zjets" in Sample:
			logging.info("Comparing to nominal ZJets sample ... ")
			NomGr = RootFile.Get(TDirectoryName + "/FTAG2_Zjets_Sherpa221/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			SystGr = RootFile.Get(TDirectoryName + "/"+Sample+"/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			UncHist = SystematicSampleWrapper(NomGr, SystGr)
		elif "Diboson" in Sample:
			logging.info("Comparing to nominal Diboson sample ... ")
			NomGr = RootFile.Get(TDirectoryName + "/FTAG2_Diboson_Sherpa222/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			SystGr = RootFile.Get(TDirectoryName + "/"+Sample+"/"+"h_"+TDirectoryName+"_"+NameCheck+"_nominal")
			UncHist = SystematicSampleWrapper(NomGr, SystGr)
		else:
			logging.info("You should write some code if you've done down this road ... ")

		logging.debug("Exporting uncert hist to syst file as tot_uncert hist ... ")
		RootFile.cd(TDirectoryName + "/" + Sample)
		UncHist.Write('', ROOT.TObject.kOverwrite)

def CalculateSystematics(NominalFile, SystematicFiles, TDirectory, PathToFiles, OutputFile, Sample, Objs, Vars, EventVar):
	NomTDir = NominalFile.Get(TDirectory)

	NominalHists = GetGoodListOfHistograms(Objs, Vars, NomTDir, TDirectory, Sample, EventVar)
	# Recreate the Objs list from the histograms found

	if len(NominalHists) != 0:
		# Ensure we have results to work with
		Objs = list(dict.fromkeys(Objs))
		for Obj in Objs:
			# Can now group them and add them
			logging.info("Organising histograms for the object " + Obj)
			for Var in Vars:
				logging.info("Getting the histograms for the variable " + Var)

				# Now search for all the matches in using this
				if "data" in Sample and EventVar:
					SearchString = '^h_'+TDirectory+'_('+Obj+')_(data)$'
				elif "data" in Sample:
					SearchString = '^h_'+TDirectory+'_('+Obj+')_('+Var+')_(data)$'
				elif EventVar:
					SearchString = '^h_'+TDirectory+'_('+Obj+')_([blc]*)$'
				else:
					SearchString = '^h_'+TDirectory+'_('+Obj+')_('+Var+')_([blc]*)$'
				r = re.compile(r'' + SearchString)

				HistMatches = filter(r.match, NominalHists)
				if len(HistMatches) == 0:
					logging.debug("No histograms found for this obj and variable, skipping!")
					continue
				logging.debug("Regex search string: " + SearchString)
				logging.debug("Found matches after regex search:")
				logging.debug(HistMatches)

				logging.debug("About to create nominal added histogram")
				NomHist = GetAddedHistogram(HistMatches, NomTDir)
				OutputFile.cd(TDirectory + "/" + Sample)
				NomHist.Write('', ROOT.TObject.kOverwrite)
				if not "data" in Sample:
					for File in SystematicFiles:
						# Now compare every systematic to the nominal
						SystFile = tfile(PathToFiles+File)
						SystematicName = File.split(".root")[0].split(Sample+"_")[1].split("_combination")[0]
						if SystFile == None: continue
						logging.debug("Successfully opened systematic file: " + File)
						logging.info("Including systematic:\t" + SystematicName)
						if len(SystFile.GetListOfKeys()) == 0: continue
						SystTDir = SystFile.Get(TDirectory)

						# Get the systematic hist and calculate uncert w.r.t to nominal
						logging.debug("About to create systematic added histogram")
						SystHist = GetAddedHistogram(HistMatches, SystTDir, SystematicName)
						logging.debug("About to create uncertainty graph")
						UncGr = GetUncertaintyGr(NomHist, SystHist)

						OutputFile.cd(TDirectory + "/" + Sample)
						UncGr.Write('', ROOT.TObject.kOverwrite)
						SystFile.Close()
		return True
	else:
		logging.info("No histograms found in file, bad file!")
		return False

def SystematicSampleWrapper(Nom_Gr, Syst_Gr):
	logging.debug("Getting nominal sample graph ... ")
	logging.debug(str(Nom_Gr)+" "+str(type(Nom_Gr)))
	logging.debug("Getting nominal alt. sample graph ... ")
	logging.debug(str(Syst_Gr)+" "+str(type(Syst_Gr)))

	logging.debug("Create corresponding uncert graph ... ")
	Unc_Gr = GetUncertaintyGr(Nom_Gr, Syst_Gr)
	HistName = ("_").join(Unc_Gr.GetName().split("_")[:-1]) + "_tot_uncert"
	Unc_Gr.SetName(HistName)

	logging.debug(str(Unc_Gr)+" "+str(type(Unc_Gr)))
	logging.info("Final combined graph: " + str(Unc_Gr.GetName()))
	return Unc_Gr

def GetUncertaintyGr(NominalHist, SystematicHist):
	# Create TGraphs required
	Gr_Nom = ROOT.TGraphAsymmErrors(NominalHist)
	Gr_Syst = ROOT.TGraphAsymmErrors(SystematicHist)
	Gr_Diff = ROOT.TGraphAsymmErrors(Gr_Nom.GetN())

	for i in range(0, Gr_Nom.GetN()):
		# Use the nom and syst to calculate the rel syst.
		if Gr_Syst.GetPointY(i) == 0.0 and Gr_Nom.GetPointY(i) == 0.0:
			Gr_Diff.SetPoint(i, Gr_Nom.GetPointX(i), 0.0)
		elif Gr_Nom.GetPointY(i) == 0.0:
			Gr_Diff.SetPoint(i, Gr_Nom.GetPointX(i), 1.0)
		else:
			Gr_Diff.SetPoint(i, Gr_Nom.GetPointX(i), ((Gr_Syst.GetPointY(i) - Gr_Nom.GetPointY(i)) / Gr_Nom.GetPointY(i)) )
	Gr_Diff.SetName("Gr_"+("_").join(SystematicHist.GetName().split("_")[1:]))
	return Gr_Diff

def GetAddedHistogram(InputHists, TDirectory, SystName="nominal"):
	# Now that we have all the histogram names for an object
	# as well as a particular variable we can get them all,
	# add them all to get the complete variable  distribution
	for Index in range (0,len(InputHists)):
		shortened_name = "_".join(InputHists[Index].split("_")[0:-1])
		if Index == 0:
			# Create clone of hist to add to
			logging.debug("Creating a histogram with name " + shortened_name + " from " + str(InputHists[Index]))
			Hist = TDirectory.Get(InputHists[Index]).Clone(shortened_name + "_" + SystName)
		else:
			# Just add it to the exising hist
			logging.debug("Adding " + str(InputHists[Index]) + " to the histogram " + shortened_name)
			HistToAdd = TDirectory.Get(InputHists[Index]).Clone()
			Hist.Add(HistToAdd)
	Hist.SetTitle("")
	return Hist

def GetGoodListOfHistograms(Objects, Variables, TDirectory, TDirectoryName, Sample, EventVar):
	NominalHists = []
	for key in TDirectory.GetListOfKeys():
		# Loop over the TDirectory and pick out the histograms we want
		KeyName = key.GetName()
		for Obj in Objects:
			for Var in Variables:
				if "data" in Sample and EventVar:
					SearchString = '^h_'+TDirectoryName+'_('+''.join(Obj)+')_(data)$'
				elif "data" in Sample:
					SearchString = '^h_'+TDirectoryName+'_('+''.join(Obj)+')_('+''.join(Var)+')_(data)$'
				elif EventVar:
					SearchString = '^h_'+TDirectoryName+'_('+''.join(Obj)+')_([blc]*)$'
				else:
					SearchString = '^h_'+TDirectoryName+'_('+''.join(Obj)+')_('+''.join(Var)+')_([blc]*)$'
				pattern = re.compile(r'' + SearchString)
				match = pattern.search(KeyName)

				if match:
					full_histname = match.group(0)
					obj = match.group(1)
				else:
					continue

				# Now have a list of good histograms
				logging.debug("Taking nominal histogram " + full_histname)
				NominalHists.append(full_histname)
	logging.debug("These are the good (nominal) histograms found that we will be using:")
	logging.debug(NominalHists)
	return NominalHists

def LogDebuggingHistInfo(Histogram, HistogramType):
	# HistogramTypes: "nom", "syst"
	if HistogramType == "nom":
		HistogramType = "nominal"
	elif HistogramType == "syst":
		HistogramType = "systematic"
	else:
		HistogramType = "unknown"
	logging.debug("----------------------------------------------------------------------------------------------")
	logging.debug("Debugging info regarding the "+HistogramType+" histogram:")
	logging.debug("\nName:\t" + str(Histogram.GetName()) + "\nType:\t" + str(type(Histogram)) + "\nObject:\t" + str(Histogram))
	logging.debug("----------------------------------------------------------------------------------------------")

def byteify(input):
  if isinstance(input, dict):
    return {byteify(key): byteify(value)
      for key, value in input.iteritems()}
  elif isinstance(input, list):
    return [byteify(element) for element in input]
  elif isinstance(input, unicode):
    return input.encode('utf-8')
  else:
    return input

def tfile(path, mode='READ'):
	if not os.path.exists(path):
		if not mode=="UPDATE":
			raise RuntimeError("{} not found!".format(path))
	tf = ROOT.TFile.Open(path, mode)
	# if tf.IsZombie():
	if tf is None:
		logging.debug("Your file is broken, you wanna fix that, its as shit like spurs")
		# raise RuntimeError("Unable to open {}!".format(path))
	return tf

def get_args():
	args = argparse.ArgumentParser(description='')
	# Arguments related to the creating the input file for the plotter
	args.add_argument('--CreateFile', action="store_true", help="If running for the first time, turn on.")
	args.add_argument('--InputPath', type=str, default="/atlas/shatlas/FTAGCalibrations/code/AlvaroCode/TTbar-b-calib-final-selection-r21/Results/r21.2.130_combined/histograms/")
	args.add_argument('--Plots', type=str, default=os.getcwd()+"/Configs/plots.json")
	args.add_argument('--Samples', type=str, default=os.getcwd()+"/Configs/samples.json")
	args.add_argument('--OutputFile', type=str, default=os.getcwd()+"/BTagHistSysts.root")
	
	# Arguments related to the plotter part of the code
	args.add_argument('--PlotFile', action="store_true", help="Run after file created")
	args.add_argument('--PlotFilename', type=str, default=os.getcwd()+"/BTagHistSysts.root", help="Name of the file to plot")
	args.add_argument('--OutputDir', type=str, default=os.getcwd()+"/Plots/", help="Path for the plots")
	args.add_argument('--NoBatch', action="store_true", help="Turn off batch mode so you see plots")
	return args.parse_args()

if __name__ == '__main__':
	main()
