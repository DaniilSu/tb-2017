#!/usr/bin/python

import sys
import os
import numpy as np
from math import *
from ROOT import *
from array import array
from scipy.interpolate import interp1d
from decimal import *
import matplotlib.pyplot as plt

def my2DLineFit( nSP, zVals, uVals, zValsErr, uValsErr ):
	pars = array('d',2 * [0])
	parsErr = array('d',2 * [0])
	# print nSP, zVals, uVals
	gr = TGraphErrors(nSP,zVals,uVals,zValsErr,uValsErr)
	gr.SetMarkerStyle(25)
	gr.SetMarkerColor(1)
	# gr.Draw()
	fit = TF1("fit", "pol1", 0, 70)
	gr.Fit(fit,"QNR")
	fit.GetParameters(pars)
	parsErr = fit.GetParErrors()
	# raw_input("Press Enter to continue...")
	return [pars, parsErr]

def vshapeParametric( x, par) :
	if isinstance(x, float) or isinstance(x, int):
		result = par[0] * (x - par[1]) ** 2 + par[2]
	else :
		xx = x[0]
		result = par[0] * (xx - par[1]) ** 2 + par[2]
	return result

def vshapeParametricMod(x, par) :
	if isinstance(x, float) or isinstance(x, int):
		result = par[0] + par[1] * x + par[2] * x ** 2 + par[3] * x ** 3 + par[4] * x ** 4 + par[5] * x ** 5 + par[6] * x ** 6
	else:
		xx = x[0]
		result = par[0] + par[1] * xx + par[2] * xx ** 2 + par[3] * xx ** 3 + par[4] * xx ** 4 + par[5] * xx ** 5 + par[6] * xx ** 6
	return result

def inverseVshapeParametric(dt, par, y_track) :
	result = 2 * [0]
	if dt < par[2]:
		result[0] = par[1]
		result[1] = par[1]
	else:
		result[0] = par[1] + sqrt((dt - par[2])/par[0])
		result[1] = par[1] - sqrt((dt - par[2])/par[0])
	if abs(result[0] - y_track) < abs(result[1] - y_track) :
		return result[0]
	else :
		return result[1]

def inverseVshapeParametricMod(dt, par1, par2, minimum1, minimum2, y_track) :
	result = 2 * [0]
	if dt < minimum1 and dt < minimum2:
		result[0] =	par1[0] + par1[1] * minimum1 + par1[2] * minimum1 ** 2 + par1[3] * minimum1 ** 3 + par1[4] * minimum1 ** 4 + par1[5] * minimum1 ** 5 + par1[6] * minimum1 ** 6
		result[1] = par2[0] + par2[1] * minimum2 + par2[2] * minimum2 ** 2 + par2[3] * minimum2 ** 3 + par2[4] * minimum2 ** 4 + par2[5] * minimum2 ** 5 + par2[6] * minimum2 ** 6
	elif dt < minimum1 and dt >= minimum2:
		result[0] =	par1[0] + par1[1] * minimum1 + par1[2] * minimum1 ** 2 + par1[3] * minimum1 ** 3 + par1[4] * minimum1 ** 4 + par1[5] * minimum1 ** 5 + par1[6] * minimum1 ** 6
		result[1] = par2[0] + par2[1] * dt + par2[2] * dt ** 2 + par2[3] * dt ** 3 + par2[4] * dt ** 4 + par2[5] * dt ** 5 + par2[6] * dt ** 6
	elif dt >= minimum1 and dt < minimum2:
		result[0] = par1[0] + par1[1] * dt + par1[2] * dt ** 2 + par1[3] * dt ** 3 + par1[4] * dt ** 4 + par1[5] * dt ** 5 + par1[6] * dt ** 6
		result[1] = par2[0] + par2[1] * minimum2 + par2[2] * minimum2 ** 2 + par2[3] * minimum2 ** 3 + par2[4] * minimum2 ** 4 + par2[5] * minimum2 ** 5 + par2[6] * minimum2 ** 6
	else:
		result[0] = par1[0] + par1[1] * dt + par1[2] * dt ** 2 + par1[3] * dt ** 3 + par1[4] * dt ** 4 + par1[5] * dt ** 5 + par1[6] * dt ** 6
		result[1] = par2[0] + par2[1] * dt + par2[2] * dt ** 2 + par2[3] * dt ** 3 + par2[4] * dt ** 4 + par2[5] * dt ** 5 + par2[6] * dt ** 6
	if abs(result[0] - y_track) < abs(result[1] - y_track) :
		return result[0]
	else :
		return result[1]

def inverseVshapeParametricInter(dt, f1, f2, minimum1, minimum2, maximum1, maximum2, y_track) :
	result = 2 * [0]
	if dt < minimum1 and dt < minimum2:
		result[0] = f1(minimum1)
		result[1] = f2(minimum2)
	elif dt < minimum1 and dt >= minimum2:
		result[0] = f1(minimum1)
		result[1] = f2(dt)
	elif dt >= minimum1 and dt < minimum2:
		result[0] = f1(dt)
		result[1] = f2(minimum2)
	#elif dt > maximum1 and dt > maximum2:
		#result[0] = f1(maximum1)
		#result[1] = f2(maximum2)
	elif dt <= maximum1 and dt > maximum2:
		result[0] = f1(dt)
		#result[1] = f2(maximum2)
		return result[0]
	elif dt > maximum1 and dt <= maximum2:
		#result[0] = f1(maximum1)
		result[1] = f2(dt)
		return result[1]
	else:
		result[0] = f1(dt)
		result[1] = f2(dt)
	if abs(result[0] - y_track) < abs(result[1] - y_track) :
		return result[0]
	else:
		return result[1]

def insertInBetween(list_old):
	list_new = []
	for i,j in enumerate(list_old):
		list_new.append(j)
		if (i+2 < len(list_old)):
			list_new.append(np.mean(list_old[i:i+2]))
		elif (i+2 == len(list_old)):
			list_new.append((list_old[i]+list_old[i+1])/2.)
		else:
			continue
	return list_new

def interpolatePoints(coord,dtime,errorcoord,errordtime,apex):
	indexDiv = next(i for i,j in enumerate(coord) if j>=apex)
	#dtime = [round(i,4) for i in dtime]
	#dtime = np.around(dtime, decimals=4)
	getcontext().prec = 4
	dtime = [float(Decimal("%.4f" % e)) for e in dtime]
	f_left = interp1d(sorted(dtime[indexDiv:]),sorted(coord[indexDiv:]),kind='slinear')
	f_right = interp1d(sorted(dtime[:indexDiv]),sorted(coord[:indexDiv],reverse=True),kind='slinear')
#	dtime_precise_left = insertInBetween(insertInBetween(insertInBetween(sorted(dtime[indexDiv:]))))
#	dtime_precise_right = insertInBetween(insertInBetween(insertInBetween(sorted(dtime[:indexDiv]))))
	print dtime[indexDiv:]
#	plt.plot(sorted(dtime[indexDiv:]), sorted(coord[indexDiv:]), 'o')
#	plt.plot(dtime_precise_left, f_left(dtime_precise_left), '-')
#	plt.legend(['data', 'linear'], loc='best')
#	plt.xlabel('drift time (ns)')
#	plt.ylabel('Y (cm)')
#	plt.title('Left branch of the V-shape')
#	plt.show()
#	plt.plot(sorted(dtime[:indexDiv]), sorted(coord[:indexDiv],reverse=True), 'o')
#	plt.plot(dtime_precise_right, f_right(dtime_precise_right), '-')
#	plt.legend(['data', 'linear'], loc='best')
#	plt.xlabel('drift time (ns)')
#	plt.ylabel('Y (cm)')
#	plt.title('Right branch of the V-shape')
#	plt.show()
	top_dtime=sorted(dtime[indexDiv:])
	bot_dtime=sorted(dtime[:indexDiv])
	return f_left, f_right, dtime[indexDiv], dtime[indexDiv-1], top_dtime[-1], bot_dtime[-1]

def findInverseFunc(coord,dtime,errorcoord,errordtime,apex) :
	indexDiv = next(i for i,j in enumerate(coord) if j>=apex)
	gre1 = TGraphErrors(len(coord[indexDiv:]),dtime[indexDiv:].flatten('C'),coord[indexDiv:].flatten('C'),errordtime[indexDiv:].flatten('C'),errorcoord[indexDiv:].flatten('C'))
	gre2 = TGraphErrors(len(coord[:indexDiv-1]),dtime[:indexDiv-1].flatten('C'),coord[:indexDiv-1].flatten('C'),errordtime[:indexDiv-1].flatten('C'),errorcoord[:indexDiv-1].flatten('C'))
	par1 = array('d',7 * [0])
	par2 = array('d',7 * [0])
	pf1 = TF1("pf1","pol6",min(dtime[indexDiv:]),max(dtime[indexDiv:]))
	pf2 = TF1("pf2","pol6",min(dtime[:indexDiv-1]),max(dtime[:indexDiv-1]))
	gre1.Fit(pf1,"WWQNR")
	pf1.GetParameters(par1)
	gre1.Fit(pf1,"QNR")
#	a = TCanvas()
#	gre1.SetMarkerColor(kBlue)
#	gre1.SetMarkerStyle(21)
#	gre1.GetXaxis().SetTitle("t_{drift} (ns)")
#	gre1.GetYaxis().SetTitle("Y (cm)")
#	gre1.SetTitle("Left branch of the V-shape")
#	gre1.Draw("AP")
#	a.SaveAs("top_branch_old.png")
	gre2.Fit(pf2,"WWQNR")
	pf2.GetParameters(par2)
	gre2.Fit(pf2,"QNR")
#	b = TCanvas()
#	gre2.SetMarkerColor(kRed)
#	gre2.SetMarkerStyle(22)
#	gre2.GetXaxis().SetTitle("t_{drift} (ns)")
#	gre2.GetYaxis().SetTitle("Y (cm)")
#	gre2.SetTitle("Right branch of the V-shape")
#	gre2.Draw("AP")
#	b.SaveAs("bottom_branch_old.png")
#	raw_input("Press any key to continue\n")
	pf1.GetParameters(par1)
	pf2.GetParameters(par2)
	return par1, par2, dtime[indexDiv], dtime[indexDiv-1]

def chooseTheBest(nt, tt, tt08, vshapeVal, sigma=None) :
	best = 0
	valid = []
	dist_old = abs(conversion*(tt[0]-tt08[0])-vshapeVal)
	if not sigma:
		for i in xrange(nt-1) :
			dist = abs(conversion*(tt[i+1]-tt08[0])-vshapeVal)
			if (dist < dist_old) :
				best = i + 1
				dist_old = dist
		valid.append(best)
		return valid
	else:
		for i in range(0,nt) :
			dist = abs(conversion*(tt[i]-tt08[0])-vshapeVal)
			if (dist < sigma):
				best = i
				valid.append(best)
		return valid

def FindRMSofSlices(input_histo, onX, firstbin, lastbin, cut, option, arr=None):
	outerAxis = input_histo.GetYaxis() if onX else input_histo.GetXaxis()
	innerAxis = input_histo.GetXaxis() if onX else input_histo.GetYaxis()

	nbins  = outerAxis.GetNbins()
	if (firstbin < 0): firstbin = 0
	if (lastbin < 0 or lastbin > nbins + 1): lastbin = nbins + 1
	if (lastbin < firstbin):
		firstbin = 0
		lastbin = nbins + 1
	opt = TString(option)
	opt.ToLower()
	ngroup = 1
	if (opt.Contains("g2")):
		ngroup = 2
		opt.ReplaceAll("g2","")
	if (opt.Contains("g3")):
		ngroup = 3
		opt.ReplaceAll("g3","")
	if (opt.Contains("g4")):
		ngroup = 4
		opt.ReplaceAll("g4","")
	if (opt.Contains("g5")):
		ngroup = 5
		opt.ReplaceAll("g5","")
	
	# implement option S sliding merge for each bin using in conjunction with a given Gn
	nstep = ngroup
	if (opt.Contains("s")):  nstep = 1
	
	npar = 1
	
	if (arr):
		arr.SetOwner()
		arr.Expand(npar)
		
	#Create one histogram for each function parameter
	name   = '\0' * 2000
	title  = '\0' * 2000
	bins = outerAxis.GetXbins()
	name = '{}_RMS'.format(input_histo.GetName())
	title = 'Resolution of the straw tube (RMS per slice)'
	if (bins.fN == 0):
		hlist = TH1D(name,title, nbins, outerAxis.GetXmin(), outerAxis.GetXmax())
	else:
		hlist = TH1D(name,title, nbins, bins.fArray)
	hlist.GetXaxis().SetTitle(outerAxis.GetTitle())
		
	#Loop on all bins in Y, generate a projection along X
	# in case of sliding merge nstep=1, i.e. do slices starting for every bin
	# now do not slices case with overflow (makes more sense)
	for bin in range(firstbin,lastbin-ngroup):
		if (onX):
			hp = input_histo.ProjectionX("_temp",bin,bin+ngroup-1,"e")
		else:
			hp = input_histo.ProjectionY("_temp",bin,bin+ngroup-1,"e")
		if (hp == 0): continue
		nentries = hp.GetEntries()
 		if (nentries == 0 or nentries < cut):
 			del hp
 			continue
		binOn = bin + ngroup/2
		hlist.Fill(outerAxis.GetBinCenter(binOn),hp.GetRMS())
		hlist.SetBinError(binOn,hp.GetRMSError())
 		del hp
	if (arr):
		arr = hlist
		return arr
	else:
		return hlist

def FitSlicesFindBinMax(vshape, onX, firstbin, lastbin, cut, option) :
	outerAxis = vshape.GetYaxis() if onX else vshape.GetXaxis()
	innerAxis = vshape.GetXaxis() if onX else vshape.GetYaxis()

	nbins  = outerAxis.GetNbins()
	if (firstbin < 0):
		firstbin = 0
	if (lastbin < 0 or lastbin > nbins + 1) :
		lastbin = nbins + 1
	if (lastbin < firstbin) :
		firstbin = 0
		lastbin = nbins + 1

	opt = TString(option)
	opt.ToLower()
	ngroup = 1
	if (opt.Contains("g2")) :
		ngroup = 2
		opt.ReplaceAll("g2","")
	if (opt.Contains("g3")) :
		ngroup = 3
		opt.ReplaceAll("g3","")
	if (opt.Contains("g4")) :
		ngroup = 4
		opt.ReplaceAll("g4","")
	if (opt.Contains("g5")) :
		ngroup = 5
		opt.ReplaceAll("g5","")

	# implement option S sliding merge for each bin using in conjunction with a given Gn
	nstep = ngroup
	if (opt.Contains("s")):
		nstep = 1
	
	# default is to fit with a gaussian
	# f1 = TF1("f1","gaus")
	# f2 = TF1("f2","gaus")
	# f3 = TF1("double_gaus", "gaus(0) + gaus(3)",innerAxis.GetXmin()+3.5,innerAxis.GetXmax()-3.5)
	# f1.SetRange(innerAxis.GetXmin()+3.5,0.)
	# f2.SetRange(0.,innerAxis.GetXmax()-3.5)

	# npar = f1.GetNpar() + f2.GetNpar()
	# if (npar <= 0):
	#	return

	# parsavef1 = array('d',f1.GetNpar() * [0])
	# parsavef2 = array('d',f2.GetNpar() * [0])
	# parsavef3 = array('d',npar * [0])
	# f1.GetParameters(parsavef1)
	# f2.GetParameters(parsavef2)
	# f3.GetParameters(parsavef3)
	
	# if (arr) :
	# 	arr.SetOwner()
	#	arr.Expand(npar + 1)
		
		
	# Create one histogram for each function parameter

	#	hlist = npar * [0]
	#	bins = outerAxis.GetXbins()
	#	for ipar in xrange(npar) :
	#		name = "%s_%d" % (vshape.GetName(),ipar)
	#		title = "Fitted value of par[%d]=%s" % (ipar,f3.GetParName(ipar))
	#		# del gDirectory.FindObject(name)
	#		if (bins.fN == 0) :
	#			hlist[ipar] = TH1D(name,title, nbins, outerAxis.GetXmin(), outerAxis.GetXmax())
	#		else :
	#			hlist[ipar] = TH1D(name,title, nbins,bins.fArray)
	#
	#		hlist[ipar].GetXaxis().SetTitle(outerAxis.GetTitle())
	#		if (arr) :
	#			arr.AddAt(hlist[ipar],ipar)

	#	name = "%s_chi2" % vshape.GetName()
	#	# del gDirectory.FindObject(name)
	#	hchi2 = 0
	#	if (bins.fN == 0) :
	#		hchi2 = TH1D(name,"chisquare", nbins, outerAxis.GetXmin(), outerAxis.GetXmax())
	#	else :
	#		hchi2 = TH1D(name,"chisquare", nbins, bins.fArray)
	#
	#	hchi2.GetXaxis().SetTitle(outerAxis.GetTitle())
	#	if (arr) :
	#		arr.AddAt(hchi2,npar)

	# Loop on all bins in Y, generate a projection along X
	# in case of sliding merge nstep=1, i.e. do slices starting for every bin
	# now do not slices case with overflow (makes more sense)
	bin = firstbin
	coord = []
	dtime = []
	errordtime = []

	while (bin+ngroup-1<=lastbin) :
		if (onX) :
  			hp = vshape.ProjectionX("_temp",bin,bin+ngroup-1,"e")
  		else :
  			hp = vshape.ProjectionY("_temp",bin,bin+ngroup-1,"e")
		if (hp == 0) :
			bin += nstep
			continue
  		nentries = hp.GetEntries()
  		kurt = hp.GetKurtosis()
  		if (nentries == 0 or nentries < cut or hp.GetMaximum() < cut) :
  			bin += nstep
  			del hp
  			continue

		if kurt < 0:
			bin += nstep
			del hp
			continue

#		hp.Rebin(3)
		binmax = hp.GetMaximumBin()
		hp.GetXaxis().SetRange(binmax-3,binmax+3)
		xmax = hp.GetMean()
		errorBinmax = hp.GetMeanError()
#		errorBinmax = hp.GetBinError(binmax)
#		xmax = hp.GetXaxis().GetBinCenter(binmax)
		dtime.append(xmax)
		coord.append(outerAxis.GetBinCenter(bin))
		errordtime.append(errorBinmax)

		#		hp_first = hp.Clone()
		#		hp_second = hp.Clone()
		#		hp_first.SetAxisRange(hp_first.GetMinimum(),0.0)
		#		hp_second.SetAxisRange(0.0,hp_second.GetMaximum())
		#		binmax_first = hp_first.GetMaximumBin()
		#		binmax_second = hp_second.GetMaximumBin()
		#		xmax_first = hp_first.GetXaxis().GetBinCenter(binmax_first)
		#		xmax_second = hp_second.GetXaxis().GetBinCenter(binmax_second)
		#		hp_first.SetMinimum(3)
		#		hp_second.SetMinimum(3)
		#		hp_first.Draw()
		#		hp_second.Draw()
		#		raw_input("Press any key to continue\n")

		# f1.SetParameters(parsavef1)
  		# f2.SetParameters(parsavef2)
		#  		f1.SetParameter(0, 100.)
		#  		f2.SetParameter(0, 100.)
		#  		f1.SetParameter(1, xmax_first)
		#  		f2.SetParameter(1, xmax_second)
		#  		f1.SetParameter(2, 0.05)
		#  		f2.SetParameter(2, 0.05)
		#		f1.SetParLimits(1,innerAxis.GetXmin()+3.5,0.0)
		#		f2.SetParLimits(1,-0.0,innerAxis.GetXmax()-3.5)
		#		f1.SetParLimits(2,0.0,0.1)
		#		f2.SetParLimits(2,0.0,0.1)
		#		f1.SetParLimits(0,cut,1000)
		#		f2.SetParLimits(0,cut,1000)
		#
		#		r1 = hp.Fit(f1,"S")
		#		r2 = hp.Fit(f2,"S+")
		#		gPad.Update()
		#
		#		if (f1.GetNumberFitPoints() > f1.GetNpar() and f2.GetNumberFitPoints() > f2.GetNpar()):
		#
		#			# parameters of first gaussian
		#			f3.SetParameter(0,r1.Parameter(0))
		#			f3.SetParameter(1,r1.Parameter(1))
		#			f3.SetParameter(2,r1.Parameter(2))
		#			# parameters of second gaussian
		#			f3.SetParameter(3,r2.Parameter(0))
		#			f3.SetParameter(4,r2.Parameter(1))
		#			f3.SetParameter(5,r2.Parameter(2))
		#			hp.Fit(f3,opt.Data())
		#			npfits = f3.GetNumberFitPoints()
		#			if (npfits > npar and npfits >= cut) :
		#				binOn = bin + ngroup/2
		#				for ipar in xrange(npar) :
		##					if ipar == 1:
		##						hlist[ipar].Fill(outerAxis.GetBinCenter(binOn),hp_first.GetMean())
		##						hlist[ipar].SetBinError(binOn,hp_first.GetMeanError())
		##					if ipar == 4:
		##						hlist[ipar].Fill(outerAxis.GetBinCenter(binOn),hp_second.GetMean())
		##						hlist[ipar].SetBinError(binOn,hp_second.GetMeanError())
		##					if ipar == 2:
		##						hlist[ipar].Fill(outerAxis.GetBinCenter(binOn),hp_first.GetRMS())
		##						hlist[ipar].SetBinError(binOn,hp_first.GetRMSError())
		##					if ipar == 5:
		##						hlist[ipar].Fill(outerAxis.GetBinCenter(binOn),hp_second.GetRMS())
		##						hlist[ipar].SetBinError(binOn,hp_second.GetRMSError())
		##					else:
		#					hlist[ipar].Fill(outerAxis.GetBinCenter(binOn),f3.GetParameter(ipar))
		#					hlist[ipar].SetBinError(binOn,f3.GetParError(ipar))
		#
		#				hchi2.SetBinContent(binOn,f3.GetChisquare()/(npfits-npar))
		#
		del hp
		bin += nstep

	return [coord, dtime, errordtime]
		#		print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		#
		#	if arr :
		#		return arr


# Usage
if len(sys.argv) < 4:
	sys.stderr.write('Usage:\t' + str(sys.argv[0]) + ' mamba_root vme_rooot numberOfRun (optional): bin_size (um) statistics_fraction (0..1.0) optional(time_window (ns)) ...\n')
	sys.exit(1)

#################################################

# bin_size flag

binSizeFlag = True

# toy MC flag

enableToyMC = True

if enableToyMC:
	enableUniform = True

# Declare histograms

maxdiff = 10
maxDiff = 2000
max1ms   = 3600000
max800ns = max1ms*1250
Clocks   =  TH2I("Clockt", "VME PC clock vs MAMBA time (1 ms)"   , 300, 0, max1ms, 300, 0, max1ms)
Clockt   =  TH2I("Clocks", "VME TDC clock vs MAMBA time (800ns)" , 300, 0, max800ns, 300, 0, max800ns)
T_ms     =  TH1I("T_ms", "MAMBA time (1 ms)" , 36000, 0, max1ms)
T_us     =  TH1I("T_us", "MAMBA time (800ns)", 36000, 0, max800ns)
TDCstamp =  TH1I("TDCstamp", "VME TDC clock (800ns)", 36000, 0, max800ns)
Pclock   =  TH1I("Pclock", "VME PC clock (1 ms)", 36000, 0, max1ms)
TSdiff   =  TH1I("TSdiff", "dif = MAMBA - VME PC clock (1 ms)", 201, -maxdiff, maxdiff)
TSdif2   =  TH2I("TSdif2", "dif vs MAMBA time (1 ms)", 300, 0, max1ms, 401, -maxdiff, maxdiff)
TSDiff   =  TH1I("TSDiff", "Dif = MAMBA - VME  TDC clock (800ns)", 201, -maxDiff, maxDiff)
TSDif2   =  TH2I("TSDif2", "Dif vs MAMBA time (800ns)", 300, 0, max800ns, 401, -maxDiff, maxDiff)

XYdut  =  TH2F("XY_dut"           , "XY_dut"            , 200, -10.0, 10.0, 200, -10.0, 10.0)
XYstr  =  TH2F("XY_dut_withstraw" , "XY_dut_withstraw"  , 200, -10.0, 10.0, 200, -10.0, 10.0)
XYtim  =  TH2F("XY_dut_withtiming", "XY_dut_withtiming" , 200, -10.0, 10.0, 200, -10.0, 10.0)
XYboth =  TH2F("XY_dut_withboth"  , "XY_dut_withboth"   , 200, -10.0, 10.0, 200, -10.0, 10.0)

XYdut_long  =  TH2F("XY_dut_long"           , "XY_dut_long"            , 200, -10.0, 10.0, 200, -10.0, 10.0)
XYstr_long  =  TH2F("XY_dut_withstraw_long" , "XY_dut_withstraw_long"  , 200, -10.0, 10.0, 200, -10.0, 10.0)
XYtim_long  =  TH2F("XY_dut_withtiming_long", "XY_dut_withtiming_long" , 200, -10.0, 10.0, 200, -10.0, 10.0)
XYboth_long =  TH2F("XY_dut_withboth_long"  , "XY_dut_withboth_long"   , 200, -10.0, 10.0, 200, -10.0, 10.0)

# VME TDC hits:
maxTDC_S = 1300.0
maxTDC_L =  500.0
minTDC   = -200.0
maxTDC_raw = 1600.0
minTDC_raw = -200.0

# VME TDC hits for short straw:

TDC_tim_trk =  TH1F("TDC_timing_scint_trk"    , "TDC_timing_scint_trk"    , 600, minTDC_raw, maxTDC_raw ) # ns
TDC_tim_trk_extra = TH1F("TDC_timing_scint_trk_extra"    , "TDC_timing_scint_trk_extra"    , 600, minTDC_raw, maxTDC_raw ) # ns

TDC_tim_raw =  TH1F("TDC_timing_scint_raw"    , "TDC_timing_scint_raw"    , 600, minTDC_raw, maxTDC_raw ) # ns
TDC_str_raw =  TH1F("TDC_straw_raw"           , "TDC_straw_raw"           , 600, minTDC_raw, maxTDC_raw ) # ns
TDC_dif_raw =  TH1F("TDC_straw_minus_tim_raw" , "TDC_straw_minus_tim_raw" , 500, -100.0,  1400.0 ) # ns
TDC_str_trk =  TH1F("TDC_straw_trk"           , "TDC_straw_trk"           , 600, minTDC_raw, maxTDC_raw ) # ns
TDC_dif_trk =  TH1F("TDC_straw_minus_tim_trk" , "TDC_straw_minus_tim_trk" , 500, -100.0,  1400.0 ) # ns

# VME TDC hits for long straw:

TDC_tim_raw_long =  TH1F("TDC_timing_scint_raw_long"    , "TDC_timing_scint_raw_long"    , 600, minTDC_raw, maxTDC_raw ) # ns
TDC_str_raw_long =  TH1F("TDC_straw_raw_long"           , "TDC_straw_raw_long"           , 600, minTDC_raw, maxTDC_raw ) # ns
TDC_dif_raw_long =  TH1F("TDC_straw_minus_tim_raw_long" , "TDC_straw_minus_tim_raw_long" , 500, -100.0,  700.0 ) # ns
TDC_str_trk_long =  TH1F("TDC_straw_trk_long"           , "TDC_straw_trk_long"           , 600, minTDC_raw, maxTDC_raw ) # ns
TDC_dif_trk_long =  TH1F("TDC_straw_minus_tim_trk_long" , "TDC_straw_minus_tim_trk_long" , 500, -100.0,  700.0 ) # ns

if binSizeFlag:
	VshapeX =  TH2F("VshapeX" , "VshapeX" , 833, -5.0, 5.0, 500, minTDC,  maxTDC_S )
	VshapeY =  TH2F("VshapeY" , "VshapeY" , 200, -1.2, 1.2, 500, minTDC,  maxTDC_S )
	VshapeY_SiPM_up = TH2F("VshapeY_uSiPMs" , "VshapeY_uSiPMs" , 200, -1.2, 1.2, 500, minTDC,  maxTDC_S )
	VshapeY_SiPM_down = TH2F("VshapeY_dSiPMs" , "VshapeY_dSiPMs" , 200, -1.2, 1.2, 500, minTDC,  maxTDC_S )
else:
	VshapeX =  TH2F("VshapeX" , "VshapeX" , int(10.0/(0.0001*float(sys.argv[4]))), -5.0, 5.0, 500, minTDC,  maxTDC_S )
	VshapeY =  TH2F("VshapeY" , "VshapeY" , int(2.4/(0.0001*float(sys.argv[4]))), -1.2, 1.2, 500, minTDC,  maxTDC_S )
	VshapeY_SiPM_up = TH2F("VshapeY_uSiPMs" , "VshapeY_uSiPMs" , int(2.4/(0.0001*float(sys.argv[4]))), -1.2, 1.2, 500, minTDC,  maxTDC_S )
	VshapeY_SiPM_down = TH2F("VshapeY_dSiPMs" , "VshapeY_dSiPMs" , int(2.4/(0.0001*float(sys.argv[4]))), -1.2, 1.2, 500, minTDC,  maxTDC_S )

if enableToyMC:
	VshapeY_toy = TH2F("VshapeY_toy", "VshapeY_toyMC", 200, -1.2, 1.2, 500, minTDC,  maxTDC_S)
	VshapeY_long_toy = TH2F("VshapeY_long_toy" , "VshapeY_long_toyMC" , 130, -0.7, 0.6, 233, minTDC,  maxTDC_L)
	sigmaYtrak = 0.006
	sigma1time = 5. # 3.6
	sigma2time = 5. # 35.

VshapeX_long =  TH2F("VshapeX_long" , "VshapeX_long" , 1000, -5.0, 5.0, 233, minTDC,  maxTDC_L )
VshapeY_long =  TH2F("VshapeY_long" , "VshapeY_long" , 130, -0.7, 0.6, 233, minTDC,  maxTDC_L )
VshapeY_long_SiPM_up = TH2F("VshapeY_long_uSiPMs" , "VshapeY_long_uSiPMs" , 130, -0.7, 0.6, 233, minTDC,  maxTDC_L )
VshapeY_long_SiPM_down = TH2F("VshapeY_long_dSiPMs" , "VshapeY_long_dSiPMs" , 130, -0.7, 0.6, 233, minTDC,  maxTDC_L )

passed_tracks = TH1F("Short_straw_passed_tracks", "Short_straw_passed_tracks", 250, -1.5, 1.5)
passed_tracks_long = TH1F("Long_straw_passed_tracks", "Long_straw_passed_tracks", 300, -1.5, 1.5)
tot_tracks = TH1F("Short_straw_all_tracks", "Short_straw_all_tracks", 250, -1.5, 1.5)
tot_tracks_long = TH1F("Long_straw_all_tracks", "Long_straw_all_tracks", 300, -1.5, 1.5)

if not binSizeFlag:
	resolutionS = TH2F("Resolution_short", "Resolution_short", 250, -1.5, 1.5, 50, -0.3, 0.3)
	if enableToyMC:
		resolutionS_toy = TH2F("Resolution_short_toy", "Resolution_short_toyMC", 250, -1.5, 1.5, 50, -0.3, 0.3)
else:
	resolutionS = TH2F("Resolution_short", "Resolution_short", int(3.0/(0.0001*float(sys.argv[4]))), -1.5, 1.5, 1000, -0.3, 0.3)
	if enableToyMC:
		resolutionS_toy = TH2F("Resolution_short_toy", "Resolution_short_toyMC", int(3.0/(0.0001*float(sys.argv[4]))), -1.5, 1.5, 1000, -0.3, 0.3)
resolutionL = TH2F("Resolution_long", "Resolution_long", 300, -1.5, 1.5, 60, -0.3, 0.3)
if enableToyMC:
	resolutionL_toy = TH2F("Resolution_long_toy", "Resolution_long_toyMC", 300, -1.5, 1.5, 60, -0.3, 0.3)

residualY_short = TH1F("varY"+"_Short_Straw", "Y_sigma_on"+"_Short_Straw", 40, 0., 0.01)
residualY_long = TH1F("varY"+"_Long_Straw", "Y_sigma_on"+"_Long_Straw", 40, 0., 0.01)
residualX_short = TH1F("varX"+"_Short_Straw", "X_sigma_on"+"_Short_Straw", 40, 0., 0.01)
residualX_long = TH1F("varX"+"_Long_Straw", "X_sigma_on"+"_Long_Straw", 40, 0., 0.01)

TDC_hit_cor = TH1F("Difference_between_2_hits_in_T0","Difference_between_2_hits_in_T0",400,-1100,100)
TDC_hit_cor_SiPM_up = TH1F("Difference_between_hits_in_T0_and_uSiPMs","Difference_between_hits_in_T0_and_uSiPMs",600,0,600)
TDC_hit_cor_SiPM_down = TH1F("Difference_between_hits_in_T0_and_dSiPMs","Difference_between_hits_in_T0_and_dSiPMs",600,0,600)
TDC_hit_cor_SiPM_up_down = TH1F("Difference_between_hits_in_uSiPMs_and_dSiPMs","Difference_between_hits_in_uSiPMs_and_dSiPMs",100,-40,40)

TDC_scin_2_hits = TH2F("Two_hits_in_T0","Two_hits_in_T0", 167, -50, 450, 600, minTDC_raw, maxTDC_raw)
TDC_scin_2_hits_SiPM_up = TH2F("Hits_in_T0_and_uSiPMs","Hits_in_T0_and_uSiPMs", 100, 0, 300, 100, 0, 300)
TDC_scin_2_hits_SiPM_down = TH2F("Hits_in_T0_and_dSiPMs","Hits_in_T0_and_dSiPMs", 100, 0, 300, 100, 0, 300)
TDC_scin_2_hits_SiPM_up_down = TH2F("Hits_in_uSiPMs_and_dSiPMs","Hits_in_uSiPMs_and_dSiPMs", 100, 0, 300, 100, 0, 300)

passed_tracks.Sumw2()
passed_tracks_long.Sumw2()
tot_tracks.Sumw2()
tot_tracks_long.Sumw2()

# noise reduction flag

noiseRedFlag = False

# batch flag

batchFlag = True

# create sync file

syncFileFlag = False

# use sync file

loadSyncFileFlag = True

# reference run flag

referenceRunFlag = False

# make plots of tracks

plotTracksFlag = False

# advanced V-shape fitting flag

advVshapeFit = True

# sigma test Flag

sigmaTest = False

# declare the histograms for residual calculations

nDets = 8

resU_Det = nDets * [0] #trick to declare an list filled with zeroes

resUvsU_Det = nDets * [0]
resUvsV_Det = nDets * [0]

resUvsSlpX_Det = nDets * [0]
resUvsSlpY_Det = nDets * [0]

for iDet in range(nDets):
	resU_Det[iDet] = TH1F("resU"+"_Det"+str(iDet), "U_residual_on"+"_Det"+str(iDet), 40, -0.1, 0.1)

	resUvsU_Det[iDet] = TH2F("resUvsU"+"_Det"+str(iDet), "U_residual_vs_U_on"+"_Det"+str(iDet), 40, -5, 5, 40, -0.1, 0.1)
	resUvsV_Det[iDet] = TH2F("resUvsV"+"_Det"+str(iDet), "U_residual_vs_V_on"+"_Det"+str(iDet), 40, -5, 5, 40, -0.1, 0.1)

	resUvsSlpX_Det[iDet] = TH2F("resUvsSlpX"+"_Det"+str(iDet), "U_residual_vs_SlpX_on"+"_Det"+str(iDet), 40, -.5, .5, 40, -0.1, 0.1)
	resUvsSlpY_Det[iDet] = TH2F("resUvsSlpY"+"_Det"+str(iDet), "U_residual_vs_SlpY_on"+"_Det"+str(iDet), 40, -.5, .5, 40, -0.1, 0.1)

vshf = "vshape_pars_"
if not binSizeFlag:
	vshf = vshf + "run_" + sys.argv[3] + ".txt"
else:
	vshf = vshf + "run_" + sys.argv[3] + "_bin_" + sys.argv[4] + "_" + sys.argv[5] + ".txt"

if not advVshapeFit:
	if referenceRunFlag:
		refFile = "vshape_pars_run_328.txt"
		ff = open(refFile,"r")
		ff.readline()
		parRefS = map(float, ff.readline().split())
		ff.readline()
		parRefL = map(float, ff.readline().split())
		ff.close()


# Output root file

if syncFileFlag and not loadSyncFileFlag:
	# out_file_name = os.environ['TBD'] + "synchronized_trees/Synchronized_data_run_" + sys.argv[3] + ".root"
	out_file_name = "Synchronized_data_run_" + sys.argv[3] + ".root"
	if not os.path.exists("rootfile/"):
		os.makedirs("rootfile/")
	out_file = TFile("rootfile/" + out_file_name, "RECREATE")
	if out_file.IsOpen():
		print 'The output file was opened successfully!'

if not binSizeFlag:
	plots_dir = "plots_with_tracks_run_" + sys.argv[3]
else:
	plots_dir = "plots_with_tracks_run_" + sys.argv[3] + "_bin_" + sys.argv[4] + "_" + sys.argv[5]
if os.path.exists(plots_dir + "/"):
	print "Directory exists!"
else:
	os.makedirs(plots_dir + "/")

if not os.path.exists("outputfiles/"):
	os.makedirs("outputfiles/")

if os.path.isfile("outputfiles/" + vshf) :
	print "Vshape file exists! Applying noise reduction algorythm..."
	noiseRedFlag = True
else :
	print "Creating vshape file..."


if referenceRunFlag:
	if not binSizeFlag:
		resol_file_name = "resolution_run_" + sys.argv[3] + ".root"
	else:
		resol_file_name = "resolution_run_" + sys.argv[3] + "_bin_" + sys.argv[4] + "_" + sys.argv[5] + ".root"
else:
	if not binSizeFlag:
		resol_file_name = "resolution_separate_run_" + sys.argv[3] + ".root"
	else:
		resol_file_name = "resolution_separate_run_" + sys.argv[3] + "_bin_" + sys.argv[4] + "_" + sys.argv[5] + ".root"
resol_file = TFile("outputfiles/" + resol_file_name, "RECREATE")
if resol_file.IsOpen():
	print "The output file with resolution studies was opened successfully!"

if noiseRedFlag:
	f = open( "outputfiles/" + vshf,"r")
	if not advVshapeFit:
		f.readline()
		parS = map(float, (f.readline()).split())
		f.readline()
		parL = map(float, (f.readline()).split())
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		bordS = map(float, (f.readline()).split())
		f.readline()
		bordL = map(float, (f.readline()).split())
	else:
		f.readline()
		parS = map(float, (f.readline()).split())
		f.readline()
		parL = map(float, (f.readline()).split())
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		bordS = map(float, (f.readline()).split())
		f.readline()
		bordL = map(float, (f.readline()).split())
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		apexS = map(float, (f.readline()).split())
		f.readline()
		apexL = map(float, (f.readline()).split())
	f.close()

	if advVshapeFit:
		if not binSizeFlag:
			f2S = open("outputfiles/" + "points_short_" + sys.argv[3] + ".dat", "r")
		else:
			f2S = open("outputfiles/" + "points_short_" + sys.argv[3] + "_bin_" + sys.argv[4] + "_" + sys.argv[5] + ".dat", "r")
		coordS,dtimeS,errorcoordS,errordtimeS = np.loadtxt(f2S, skiprows=1, unpack=True)
		f2S.close()
		f2L = open("outputfiles/" + "points_long_" + sys.argv[3] + ".dat", "r")
		coordL,dtimeL,errorcoordL,errordtimeL = np.loadtxt(f2L, skiprows=1, unpack=True)
		f2L.close()
#		parSinvTop, parSinvBot, SinvTopMin, SinvBotMin = findInverseFunc(coordS,dtimeS,errorcoordS,errordtimeS,apexS)
#		parLinvTop, parLinvBot, LinvTopMin, LinvBotMin = findInverseFunc(coordL,dtimeL,errorcoordL,errordtimeL,apexL)
		fSinterTop, fSinterBot, SinterTopMin, SinterBotMin, SinterTopMax, SinterBotMax = interpolatePoints(coordS,dtimeS,errorcoordS,errordtimeS,apexS)
		fLinterTop, fLinterBot, LinterTopMin, LinterBotMin, LinterTopMax, LinterBotMax = interpolatePoints(coordL,dtimeL,errorcoordL,errordtimeL,apexL)

	plots_dir = plots_dir + "_noise_red/"
	if os.path.exists(plots_dir):
		print "Directory exists!"
	else:
		os.makedirs(plots_dir)
else :
	plots_dir = plots_dir + "/"

def improveNoise(dirname, vshf, tracks_RMS, parRef=None):
	global VshapeY
	global VshapeY_long
	global referenceRunFlag
	parS = array('d',3 * [0])
	parL = array('d',3 * [0])
	parSerr = array('d',3 * [0])
	parLerr = array('d',3 * [0])
	[coordS,dtimeS,errordtimeS] = FitSlicesFindBinMax(VshapeY, False, VshapeY.GetXaxis().GetFirst(), VshapeY.GetXaxis().GetLast(),10,"QNR") #8
	[coordL,dtimeL,errordtimeL] = FitSlicesFindBinMax(VshapeY_long, False, VshapeY_long.GetXaxis().GetFirst(),VshapeY_long.GetXaxis().GetLast(),12,"QNR") #10
	errorcoordS = array('d',len(coordS)*[tracks_RMS[0]])
	errorcoordL = array('d',len(coordL)*[tracks_RMS[1]])
	pfS = TF1("pfS", "[0]*(x-[1])**2+[2]", coordS[0], coordS[-1])
	pfL = TF1("pfL", "[0]*(x-[1])**2+[2]", coordL[0], coordL[-1])
	pfS.SetParLimits(0,300,1000)
	pfS.SetParLimits(1,0.0,0.25)
	pfS.SetParLimits(2,-100,100)
	pfS.SetParameters(600.,0.1,-30.)
	pfS.SetParNames("a","x0","c")
	pfL.SetParLimits(0,300,1000)
	pfL.SetParLimits(1,-0.2,0.1)
	pfL.SetParLimits(2,-100,100)
	pfL.SetParameters(600.,0.01,-30.)
	pfL.SetParNames("a","x0","c")
	if referenceRunFlag:
		pfS.FixParameter(0,parRef[0])
		pfL.FixParameter(0,parRef[1])
	par_1S = TGraphErrors(len(coordS),np.array(coordS).flatten('C'),np.array(dtimeS).flatten('C'),errorcoordS,np.array(errordtimeS).flatten('C'))
	par_1S.Fit(pfS,"QNR")
	# a = TCanvas()
	# par_1S.SetMarkerColor(kBlue)
	# par_1S.SetMarkerStyle(21)
	# par_1S.Draw("AP")
	par_1L = TGraphErrors(len(coordL),np.array(coordL).flatten('C'),np.array(dtimeL).flatten('C'),errorcoordL,np.array(errordtimeL).flatten('C'))
	par_1L.Fit(pfL,"QNR")
	# b = TCanvas()
	# par_1L.SetMarkerColor(kBlue)
	# par_1L.SetMarkerStyle(21)
	# par_1L.Draw("AP")
	# raw_input("Press any key to continue\n")
	pfS.GetParameters(parS)
	parSerr = pfS.GetParErrors()
	pfL.GetParameters(parL)
	parLerr = pfL.GetParErrors()
	outPars = open(dirname + vshf, "w")
	outPars.write( "Short straw parameters y = a * (x - x0)^2 + c \n")
	outPars.write( str(parS[0]) + "   " + str(parS[1]) + "   " + str(parS[2]) + "\n")
	outPars.write( "Long straw parameters y = a * (x - x0)^2 + c \n")
	outPars.write( str(parL[0]) + "   " + str(parL[1]) + "   " + str(parL[2]) + "\n")
	outPars.write( "Short straw chi^2 and ndf of the fit \n")
	outPars.write( str(pfS.GetChisquare()) + "   " + str(pfS.GetNDF()) + "\n")
	outPars.write( "Long straw chi^2 and ndf of the fit \n")
	outPars.write( str(pfL.GetChisquare()) + "   " + str(pfL.GetNDF()) + "\n")
	outPars.write( "Short straw tube borders \n" )
	outPars.write( str(coordS[0]) + "   " + str(coordS[-1]) + "\n")
	outPars.write( "Long straw tube borders \n")
	outPars.write( str(coordL[0]) + "   " + str(coordL[-1]) + "\n")
	outPars.write( "Short straw errors of parameters a, x0, c \n")
	outPars.write( str(parSerr[0]) + "   " + str(parSerr[1]) + "   " + str(parSerr[2]) + "\n")
	outPars.write( "Long straw errors of parameters a, x0, c \n")
	outPars.write( str(parLerr[0]) + "   " + str(parLerr[1]) + "   " + str(parLerr[2]) + "\n")
	outPars.write( "Short straw tube drift time difference at edges and error \n")
	outPars.write( str(dtimeS[0] - dtimeS[-1]) + "   " + str(sqrt(errordtimeS[0]**2 + errordtimeS[-1]**2)) + "\n")
	outPars.write( "Long straw tube drift time difference at edges and error \n")
	outPars.write( str(dtimeL[0] - dtimeL[-1]) + "   " + str(sqrt(errordtimeL[0]**2 + errordtimeL[-1]**2)) + "\n")
	outPars.close()

def improveNoiseMod(dirname, vshf, nRun, tracks_RMS, binSize=None, stat_fraction=1.0):
	global VshapeY
	global VshapeY_long
	global binSizeFlag
	parS = array('d',7 * [0])
	parL = array('d',7 * [0])
	parSerr = array('d',7 * [0])
	parLerr = array('d',7 * [0])
	[coordS,dtimeS,errordtimeS] = FitSlicesFindBinMax(VshapeY, False, VshapeY.GetXaxis().GetFirst(), VshapeY.GetXaxis().GetLast(),12,"QNR")
	[coordL,dtimeL,errordtimeL] = FitSlicesFindBinMax(VshapeY_long, False, VshapeY_long.GetXaxis().GetFirst(),VshapeY_long.GetXaxis().GetLast(),5,"QNR") # 14
	errorcoordS = array('d',len(coordS)*[tracks_RMS[0]])
	errorcoordL = array('d',len(coordL)*[tracks_RMS[1]])
	print "Straw tube borders: %f .. %f" % (coordS[0], coordS[-1])
	pfS = TF1("pfS", "pol6", coordS[0], coordS[-1])
	pfL = TF1("pfL", "pol6", coordL[0], coordL[-1])
	par_1S = TGraphErrors(len(coordS),np.array(coordS).flatten('C'),np.array(dtimeS).flatten('C'),errorcoordS,np.array(errordtimeS).flatten('C'))
	par_1S.Fit(pfS,"WWQNR")
	par_1S.Fit(pfS,"QNR")
	#a = TCanvas()
	#par_1S.SetMarkerColor(kBlue)
	#par_1S.SetMarkerStyle(21)
	#par_1S.Draw("AP")
	par_1L = TGraphErrors(len(coordL),np.array(coordL).flatten('C'),np.array(dtimeL).flatten('C'),errorcoordL,np.array(errordtimeL).flatten('C'))
	par_1L.Fit(pfL,"WWQNR")
	par_1L.Fit(pfL,"QNR")
	#b = TCanvas()
	#par_1L.SetMarkerColor(kBlue)
	#par_1L.SetMarkerStyle(21)
	#par_1L.Draw("AP")
	# raw_input("Press any key to continue\n")
	pfS.GetParameters(parS)
	parSerr = pfS.GetParErrors()
	pfL.GetParameters(parL)
	parLerr = pfL.GetParErrors()
	apexS = pfS.GetMinimumX(coordS[0], coordS[-1])
	apexL = pfL.GetMinimumX(coordL[0], coordL[-1])
	outPars = open(dirname + vshf, "w")
	outPars.write( "Short straw parameters a0 .. a6 \n")
	outPars.write( str(parS[0]) + "  " + str(parS[1]) + "  " + str(parS[2]) + "  " + str(parS[3]) + "  " + str(parS[4]) + "  " + str(parS[5]) + "  " + str(parS[6]) + "\n")
	outPars.write( "Long straw parameters a0 .. a6 \n")
	outPars.write( str(parL[0]) + "  " + str(parL[1]) + "  " + str(parL[2]) + "  " + str(parL[3]) + "  " + str(parL[4]) + "  " + str(parL[5]) + "  " + str(parL[6]) + "\n")
	outPars.write( "Short straw chi^2 and ndf of the fit \n")
	outPars.write( str(pfS.GetChisquare()) + "   " + str(pfS.GetNDF()) + "\n")
	outPars.write( "Long straw chi^2 and ndf of the fit \n")
	outPars.write( str(pfL.GetChisquare()) + "   " + str(pfL.GetNDF()) + "\n")
	outPars.write( "Short straw tube borders \n" )
	outPars.write( str(coordS[0]) + "   " + str(coordS[-1]) + "\n")
	outPars.write( "Long straw tube borders \n")
	outPars.write( str(coordL[0]) + "   " + str(coordL[-1]) + "\n")
	outPars.write( "Short straw errors of parameters a0 .. a6 \n")
	outPars.write( str(parSerr[0]) + "  " + str(parSerr[1]) + "  " + str(parSerr[2]) + "  " + str(parSerr[3]) + "  " + str(parSerr[4]) + "  " + str(parSerr[5]) + "  " + str(parSerr[6]) + "\n")
	outPars.write( "Long straw errors of parameters a0 .. a6 \n")
	outPars.write( str(parLerr[0]) + "  " + str(parLerr[1]) + "  " + str(parLerr[2]) + "  " + str(parLerr[3]) + "  " + str(parLerr[4]) + "  " + str(parLerr[5]) + "  " + str(parLerr[6]) + "\n")
	outPars.write( "Short straw tube drift time difference at edges and error \n")
	outPars.write( str(dtimeS[0] - dtimeS[-1]) + "   " + str(sqrt(errordtimeS[0]**2 + errordtimeS[-1]**2)) + "\n")
	outPars.write( "Long straw tube drift time difference at edges and error \n")
	outPars.write( str(dtimeL[0] - dtimeL[-1]) + "   " + str(sqrt(errordtimeL[0]**2 + errordtimeL[-1]**2)) + "\n")
	outPars.write( "Short straw apex \n")
	outPars.write( str(apexS) + "\n")
	outPars.write( "Long straw apex \n")
	outPars.write( str(apexL) + "\n")
	outPars.close()
	if not binSizeFlag:
		vshapePointsS = open(dirname + "points_short_" + nRun + "_" + stat_fraction + ".dat", "w")
	else:
		vshapePointsS = open(dirname + "points_short_" + nRun + "_bin_" + binSize + "_" + stat_fraction + ".dat", "w")
	vshapePointsS.write("Y (cm)   dTime (ns)   Error Y (cm)   Error dTime (ns) \n")
	for i in xrange(len(coordS)):
		vshapePointsS.write("%f %f %f %f\n" % (coordS[i], dtimeS[i], errorcoordS[i], errordtimeS[i]))
	vshapePointsS.close()
	vshapePointsL = open(dirname + "points_long_" + nRun + ".dat", "w")
	vshapePointsL.write("Y (cm)   dTime (ns)   Error Y (cm)   Error dTime (ns) \n")
	for i in xrange(len(coordL)):
		vshapePointsL.write("%f %f %f %f\n" % (coordL[i], dtimeL[i], errorcoordL[i], errordtimeL[i]))
	vshapePointsL.close()

# Read the filenames form sys.argv and open the tree
mamba_tree = TChain("tree")
vme_tree = TChain("ADC1")
if loadSyncFileFlag and not syncFileFlag:
	for (i,filename) in enumerate(sys.argv):
		if i == 0:
			continue
		elif i == 1:
			mamba_tree.Add(filename+"/tree")
			vme_tree.Add(filename+"/ADC1")
		elif i == 2:
			continue
else:
	for (i,filename) in enumerate(sys.argv):
		if i == 0:
			continue
		elif i == 1:
			mamba_tree.Add(filename)
		elif i == 2:
			vme_tree.Add(filename)

#################################################

# Read the total number of events in the TChain
n_mamba_events = mamba_tree.GetEntries()
print "Number of MAMBA Entries:\t", n_mamba_events, "\n"
n_vme_events = vme_tree.GetEntries()
print "Number of VME Entries:\t", n_vme_events, "\n"

if syncFileFlag and not loadSyncFileFlag:
	mamba_tree.SetBranchStatus("*",0)
	mamba_tree.SetBranchStatus("nevent",1)
	mamba_tree.SetBranchStatus("timestamp",1)
	mamba_tree.SetBranchStatus("nDet",1)
	mamba_tree.SetBranchStatus("ntrk",1)
	mamba_tree.SetBranchStatus("chi2",1)
	mamba_tree.SetBranchStatus("ndof",1)
	mamba_tree.SetBranchStatus("fitStatus",1)
	mamba_tree.SetBranchStatus("ItpX",1)
	mamba_tree.SetBranchStatus("ItpY",1)
	mamba_tree.SetBranchStatus("SlpX",1)
	mamba_tree.SetBranchStatus("SlpY",1)
	mamba_tree.SetBranchStatus("XCov00",1)
	mamba_tree.SetBranchStatus("XCov01",1)
	mamba_tree.SetBranchStatus("XCov10",1)
	mamba_tree.SetBranchStatus("XCov11",1)
	mamba_tree.SetBranchStatus("YCov00",1)
	mamba_tree.SetBranchStatus("YCov01",1)
	mamba_tree.SetBranchStatus("YCov10",1)
	mamba_tree.SetBranchStatus("YCov11",1)
	mamba_tree.SetBranchStatus("XrecoT0",1)
	mamba_tree.SetBranchStatus("YrecoT0",1)
	mamba_tree.SetBranchStatus("XrecoT1",1)
	mamba_tree.SetBranchStatus("YrecoT1",1)
	mamba_tree.SetBranchStatus("XrecoT2",1)
	mamba_tree.SetBranchStatus("YrecoT2",1)
	mamba_tree.SetBranchStatus("XrecoT3",1)
	mamba_tree.SetBranchStatus("YrecoT3",1)
	mamba_tree.SetBranchStatus("XfitT0",1)
	mamba_tree.SetBranchStatus("YfitT0",1)
	mamba_tree.SetBranchStatus("XfitT1",1)
	mamba_tree.SetBranchStatus("YfitT1",1)
	mamba_tree.SetBranchStatus("XfitT2",1)
	mamba_tree.SetBranchStatus("YfitT2",1)
	mamba_tree.SetBranchStatus("XfitT3",1)
	mamba_tree.SetBranchStatus("YfitT3",1)
	mamba_tree.SetBranchStatus("nIntersect",1)
	mamba_tree.SetBranchStatus("intTrkID",1)
	mamba_tree.SetBranchStatus("intLayer",1)
	mamba_tree.SetBranchStatus("intInside",1)
	mamba_tree.SetBranchStatus("intXPos",1)
	mamba_tree.SetBranchStatus("intYPos",1)
	mamba_tree.SetBranchStatus("intZPos",1)
	mamba_tree.SetBranchStatus("intUPos",1)
	mamba_tree.SetBranchStatus("intVPos",1)
	mamba_tree.SetBranchStatus("nSP",1)
	mamba_tree.SetBranchStatus("spLayer",1)
	mamba_tree.SetBranchStatus("spXPos",1)
	mamba_tree.SetBranchStatus("spYPos",1)
	mamba_tree.SetBranchStatus("spZPos",1)
	mamba_tree.SetBranchStatus("spXErr",1)
	mamba_tree.SetBranchStatus("spYErr",1)
	mamba_tree.SetBranchStatus("spZErr",1)
	mamba_tree.SetBranchStatus("spUPos",1)
	mamba_tree.SetBranchStatus("spVPos",1)
	mamba_tree.SetBranchStatus("spIsOnTrk",1)

	vme_tree.SetBranchStatus("*",0)
	vme_tree.SetBranchStatus("PCclock",1)
	vme_tree.SetBranchStatus("pattern",1)
	vme_tree.SetBranchStatus("TDCtrig",1)
	vme_tree.SetBranchStatus("t",1)
	vme_tree.SetBranchStatus("nt00",1)
	vme_tree.SetBranchStatus("tt00",1)
	vme_tree.SetBranchStatus("nt01",1)
	vme_tree.SetBranchStatus("tt01",1)
	vme_tree.SetBranchStatus("nt02",1)
	vme_tree.SetBranchStatus("tt02",1)
	vme_tree.SetBranchStatus("nt03",1)
	vme_tree.SetBranchStatus("tt03",1)
	vme_tree.SetBranchStatus("nt04",1)
	vme_tree.SetBranchStatus("tt04",1)
	vme_tree.SetBranchStatus("nt05",1)
	vme_tree.SetBranchStatus("tt05",1)
	vme_tree.SetBranchStatus("nt06",1)
	vme_tree.SetBranchStatus("tt06",1)
	vme_tree.SetBranchStatus("nt08",1)
	vme_tree.SetBranchStatus("tt08",1)

	sync_mamba_tree = mamba_tree.CloneTree(0)
	sync_vme_tree = vme_tree.CloneTree(0)

# z position of detectors

zPos_x = [10., 30.0, 40.0, 60.0] # cm
zPos_y = [20., 0., 50., 70.] # cm
dzPos_x = [-0.0, -0.0, -0.0, -0.0]
dzPos_y = [-0., -0., -0., -0.]
# fullZlength = 70.0 # cm
straw_zPos = -24.0 # cm
straw_zPos_long = -10.0 # cm
X_zvals = array('d', 6 * [0])
Y_zvals = array('d', 6 * [0])
Xtrak = array('d', 6 * [0])
Ytrak = array('d', 6 * [0])
Xtrak_test = array('d', 6 * [0])
Ytrak_test = array('d', 6 * [0])
Xmeas = 4 * [0]
Ymeas = 4 * [0]
Yerrs = [ 0.005, 0.005, 0.005, 0.005] # estimated
Xerrs = [ 0.005, 0.005, 0.005, 0.005] # estimated
Zerrs = [0.1, 0.1, 0.1, 0.1]

if sigmaTest:
	distSigma = float(sys.argv[6]) # ns
else:
	distSigma = None

for i in range(6):
	if i < 4:
		X_zvals[i] = (zPos_x[i] + dzPos_x[i])
		Y_zvals[i] = (zPos_y[i] + dzPos_y[i])
	elif i == 4:
		X_zvals[i] = straw_zPos
		Y_zvals[i] = straw_zPos
	else:
		X_zvals[i] = straw_zPos_long
		Y_zvals[i] = straw_zPos_long


# TDC conversion for hits: 1024 is 25ns

conversion = 25.0/1024
max27bits = 134217728
ns_per_cycle = 800*max27bits
count_processed = 0
cycles = 0

# improvement factor for PCclock and the time slope for it

improvefactor = 1000
tslope = 320000 - 1.2 # default 1.2
slope = 400005000
tsdiffCut = 50 #set to very large because does not seem to be stable
tsDiffCut = 500

vme_missed_event = 0
mamba_tree.GetEntry(0)

packID = mamba_tree.nevent
vmeID = packID - 1
vme_tree.GetEntry(vmeID)
pattern = vme_tree.pattern

vmeTmp = 0
packID_old = packID
packID_diff = 0
if not loadSyncFileFlag:
	if (pattern != 0) :
		print "1"
		for k in xrange(50000 - vmeID - 1) :
			vme_tree.GetEntry(k + vmeID + 1)
			pattern = vme_tree.pattern
			vmeTmp = k + vmeID + 1
			if (pattern == 0) :
				break
		vmeID = vmeTmp
		print 'vmeID : %d' % vmeID

vmeTmp = vmeID

# create aliases for variables of mamba_tree

# clock

timestamp = mamba_tree.timestamp

# tracks

k_x = 0
k_y = 0
good_track = 0
very_good_track = 0
very_good_track_long = 0
ultra_good_track = 0

# create aliases for variables of vme_tree

PCclock = vme_tree.PCclock
TDCtrig = vme_tree.TDCtrig

timestamp_0 = timestamp
TDCstamp_0 = TDCtrig
PCclock_0 = PCclock
tsDiff = 0
tsDiff_old = 0
mambaTmp = 0
if not loadSyncFileFlag:
	for i in xrange(9):
		mamba_tree.GetEntry(i)
		timestamp = mamba_tree.timestamp
		packID = mamba_tree.nevent
		vmeID = packID - 1
		vme_tree.GetEntry(vmeID)
		pattern = vme_tree.pattern
		TDCtrig = vme_tree.TDCtrig
		if (pattern == 0) :
			if i == 0:
				TDCstamp_i_1 = TDCstamp_0
				time_us_i_1 = 0
			timestamp_iEvt = (timestamp - timestamp_0)
			time_us_i   = (timestamp_iEvt * improvefactor)/tslope
			TDCstamp_i  = (TDCtrig - TDCstamp_0)
			print 'i:  %d , MAMBA_time:   %f , diff:   %f , VME_TDC:   %f , diff:   %f' % (i, time_us_i, time_us_i - time_us_i_1, TDCstamp_i, TDCstamp_i - TDCstamp_i_1)
			TDCstamp_i_1 = TDCstamp_i
			time_us_i_1 = time_us_i
mamba_tree.GetEntry(0)
packID = mamba_tree.nevent
if not loadSyncFileFlag:
	for i in xrange(9) :
		tsDiff_old = tsDiff
		packID_old = packID
		mamba_tree.GetEntry(i+1)
		timestamp = mamba_tree.timestamp
		packID = mamba_tree.nevent
		packID_diff = packID - packID_old
		vmeID = vmeTmp + packID_diff - vme_missed_event # packID starts at 1, for 1=0
		vme_tree.GetEntry(vmeID)
		pattern = vme_tree.pattern
		if (pattern != 0) :
			print "1"
			for k in xrange(n_vme_events - vmeID - 1) :
				vme_tree.GetEntry(k + vmeID + 1)
				pattern = vme_tree.pattern
				vmeTmp = k + vmeID + 1
				if (pattern == 0):
					break

			vmeID = vmeTmp
		
		TDCtrig = vme_tree.TDCtrig
		vmeTmp = vmeID + vme_missed_event
		timestamp_iEvt = (timestamp - timestamp_0)
		time_us_i   = (timestamp_iEvt * improvefactor)/tslope
		TDCstamp_i  = (TDCtrig - TDCstamp_0)  # assuming one cycle is 2^27 = 134217728
		tsDiff      =  time_us_i - TDCstamp_i
		print 'tsDiff_old : %f' % tsDiff_old
		print 'tsDiff : %f' % tsDiff
		if (tsDiff< -tsDiffCut):
			vme_missed_event += 1
			print 'vme_missed_event upped to %d' % vme_missed_event
			print 'i     vmeID    packID       TDCtrig  time_us_i  TDCstamp_i  tsDiff'
			print '%d    %d       %d           %f       %f         %f          %f' % (i, vmeID, packID, TDCtrig, time_us_i, TDCstamp_i, tsDiff)

		if (tsDiff> tsDiffCut) :
			vme_missed_event -= 1
			print 'vme_missed_event downed to %d' % vme_missed_event
			print 'i     vmeID    packID       TDCtrig  time_us_i  TDCstamp_i  tsDiff'
			print '%d    %d       %d           %f       %f         %f          %f' % (i, vmeID, packID, TDCtrig, time_us_i, TDCstamp_i, tsDiff)

		if (abs(tsDiff - tsDiff_old) < 50) :
			mambaTmp = i+1
			break
		if i == 8 :
			mambaTmp = i+1

nnn = 0
ntics = 10
countdown = 10
last_ntics_timestamp = []
last_ntics_PCclock = []
last_ntics_TDCstamp = []
sum_last_ntics_timestamp = 0
sum_last_ntics_PCclock = 0
sum_last_ntics_TDCstamp = 0

PCclock = vme_tree.PCclock

print 'First event (i = %d and vmeID = %d) :' % (mambaTmp, vmeID)
print 'MAMBA Timestamp 0 is ', timestamp
print 'VME PCclock 0 is ', PCclock
print 'TDCstamp 0 is ', TDCtrig
print 'Improve factor is ', improvefactor
print 'Using tslope', tslope
print 'Using slope', slope

timestamp_0 = timestamp
TDCstamp_0 = TDCtrig
PCclock_0 = PCclock

n_mamba_events = mamba_tree.GetEntries()
n_vme_events = vme_tree.GetEntries()

N_tracks_tot = 0.
N_tracks_straw = 0.
N_passed_2nd_hits = 0
N_tracks_straw_long = 0.
r = TRandom()
resU_test = array('d', nDets * [0])
XCov = np.zeros((2, 2))
YCov = np.zeros((2, 2))

r = TRandom3()

for iEvt in xrange(n_mamba_events-mambaTmp):

	if abs(vme_missed_event) > 100000:
		break
	
	coin = r.Rndm()
	if coin > float(sys.argv[5]): continue
	
	if (iEvt + mambaTmp + 1) % 100 == 0:
		sys.stdout.write("Processing event:\t" + str(iEvt + mambaTmp +1)+ "\r" )
		sys.stdout.flush()

	TDCtrig_old = TDCtrig
	packID_old = packID
	mamba_tree.GetEntry(iEvt + mambaTmp)
	packID = mamba_tree.nevent
	packID_diff = packID - packID_old

	# clock

	vmeID = vmeTmp + packID_diff - vme_missed_event
	vme_tree.GetEntry(vmeID)
	if not loadSyncFileFlag:
		pattern = vme_tree.pattern
		if (pattern != 0) :
			for k in xrange(n_vme_events - vmeID - 1) :
				vme_tree.GetEntry(k + vmeID + 1)
				pattern = vme_tree.pattern
				vmeTmp = k + vmeID + 1
				if (pattern == 0):
					break

			vmeID = vmeTmp

		vmeTmp = vmeID + vme_missed_event

		timestamp = mamba_tree.timestamp
		PCclock = vme_tree.PCclock
		TDCtrig = vme_tree.TDCtrig

		PCclock_i = PCclock - PCclock_0
		timestamp_iEvt = timestamp - timestamp_0
		time_ms_i = timestamp_iEvt * improvefactor/slope
		time_us_i = timestamp_iEvt * improvefactor/tslope # us means scaling to 800 ns

		if (TDCtrig < TDCtrig_old):
			cycles += 1
			print 'i = %d, TDCtrig = %f, TDCtrig_old = %f, VME cycles = %d' % (iEvt + mambaTmp, TDCtrig, TDCtrig_old, cycles)

		TDCstamp_i  = (TDCtrig + (cycles*max27bits) - TDCstamp_0)  # assuming one cycle is 2^27 = 134217728
		tsdiff      =  time_ms_i - PCclock_i
		tsDiff      =  time_us_i - TDCstamp_i

		if abs(tsDiff) > tsDiffCut:
			cycles_from_ms_clock = int(PCclock_i * 1000000 / ns_per_cycle)
			if cycles_from_ms_clock > cycles:
				print 'i = %d, cycles = %d adjusted to %d' % (iEvt + mambaTmp, cycles, cycles_from_ms_clock)
				cycles = cycles_from_ms_clock
				TDCstamp_i  =  TDCtrig + cycles * max27bits - TDCstamp_0
				tsdiff      =  time_ms_i - PCclock_i
				tsDiff      =  time_us_i - TDCstamp_i

		# this requires a good calibration of PC clock scale vs mamba cliock scale (slope)

		if tsdiff< -tsdiffCut :
			vme_missed_event += 1
			print 'vme_missed_event upped to %d' % vme_missed_event
			print 'i     vmeID    packID    time_ms_i    pcclock_i    tsdiff    slope     tslope     TDCtrig    time_us_i    TDCstamp_i      tsDiff'
			print '%d    %d       %d        %f         %f         %f      %f     %f      %f       %f         %f          %f' % (iEvt, vmeID, packID, time_ms_i, PCclock_i, tsdiff, slope, tslope, TDCtrig, time_us_i, TDCstamp_i, tsDiff)
			countdown -= 1
			continue

		if tsdiff> tsdiffCut :
			vme_missed_event -= 1
			print 'vme_missed_event downed to %d' % vme_missed_event
			print 'i     vmeID    packID      time_ms_i     pcclock_i     tsdiff      slope      tslope       TDCtrig     time_us_i     TDCstamp_i     tsDiff'
			print '%d    %d       %d        %f         %f         %f      %f      %f       %f       %f         %f          %f' % (iEvt, vmeID, packID, time_ms_i, PCclock_i, tsdiff, slope, tslope, TDCtrig, time_us_i, TDCstamp_i, tsDiff)
			countdown -= 1
			continue

		if ( tsDiff< -tsDiffCut):
			vme_missed_event += 1
			print 'vme_missed_event upped to %d' % vme_missed_event
			print 'i     vmeID    packID    time_ms_i    pcclock_i    tsdiff    slope     tslope     TDCtrig    time_us_i    TDCstamp_i      tsDiff'
			print '%d    %d       %d        %f         %f         %f      %f     %f      %f       %f         %f          %f' % (iEvt, vmeID, packID, time_ms_i, PCclock_i, tsdiff, slope, tslope, TDCtrig, time_us_i, TDCstamp_i, tsDiff)
			countdown -= 1
			continue

		if ( tsDiff> tsDiffCut) :
			vme_missed_event -= 1
			print 'vme_missed_event downed to %d' % vme_missed_event
			print 'i     vmeID    packID      time_ms_i     pcclock_i     tsdiff      slope      tslope       TDCtrig     time_us_i     TDCstamp_i     tsDiff'
			print '%d    %d       %d        %f         %f         %f      %f      %f       %f       %f         %f          %f' % (iEvt, vmeID, packID, time_ms_i, PCclock_i, tsdiff, slope, tslope, TDCtrig, time_us_i, TDCstamp_i, tsDiff)
			countdown -= 1
			continue

		# make running slope
		# make running tslope
		if len(last_ntics_timestamp) < ntics:
			last_ntics_timestamp.append(timestamp_iEvt)
			sum_last_ntics_timestamp += timestamp_iEvt
			last_ntics_TDCstamp.append(TDCstamp_i)
			last_ntics_PCclock.append(PCclock_i)
			sum_last_ntics_TDCstamp += TDCstamp_i
			sum_last_ntics_PCclock += PCclock_i
		else:
			sum_last_ntics_timestamp -= last_ntics_timestamp[0]
			last_ntics_timestamp.remove(last_ntics_timestamp[0])
			last_ntics_timestamp.append(timestamp_iEvt)
			sum_last_ntics_timestamp += timestamp_iEvt
			sum_last_ntics_TDCstamp -= last_ntics_TDCstamp[0]
			last_ntics_TDCstamp.remove(last_ntics_TDCstamp[0])
			sum_last_ntics_PCclock -= last_ntics_PCclock[0]
			last_ntics_PCclock.remove(last_ntics_PCclock[0])
			last_ntics_TDCstamp.append(TDCstamp_i)
			sum_last_ntics_TDCstamp += TDCstamp_i
			last_ntics_PCclock.append(PCclock_i)
			sum_last_ntics_PCclock += PCclock_i
			tslope = (improvefactor * sum_last_ntics_timestamp * 1.0 / sum_last_ntics_TDCstamp)
			slope = (improvefactor * sum_last_ntics_timestamp / sum_last_ntics_PCclock)
			# print 'Adjusted ms slope = ', slope

		Clocks.Fill(time_ms_i,PCclock_i)
		Clockt.Fill(time_us_i,TDCstamp_i)
		T_ms.Fill(time_ms_i)
		T_us.Fill(time_us_i)
		Pclock.Fill(PCclock_i)
		TDCstamp.Fill(TDCstamp_i)
		TSdiff.Fill(tsdiff)
		TSdif2.Fill(time_ms_i , tsdiff )
		TSDiff.Fill(tsDiff)
		TSDif2.Fill(time_us_i , tsDiff )

	else :
		vmeID = iEvt + mambaTmp
		vme_tree.GetEntry(vmeID)


	if syncFileFlag and not loadSyncFileFlag:
		sync_mamba_tree.Fill()
		sync_vme_tree.Fill()

	#tracks

	ntrk = mamba_tree.ntrk

	nt00 = vme_tree.nt00
	nt01 = vme_tree.nt01
	nt03 = vme_tree.nt03
	nt04 = vme_tree.nt04
	nt05 = vme_tree.nt05
	nt06 = vme_tree.nt06

	nt08 = vme_tree.nt08

	tt00 = vme_tree.tt00
	tt01 = vme_tree.tt01
	tt03 = vme_tree.tt03
	tt04 = vme_tree.tt04
	tt05 = vme_tree.tt05
	tt06 = vme_tree.tt06

	tt08 = vme_tree.tt08

	# spacepoints
	nSP = mamba_tree.nSP

	k = 0
	l = 0

	'''if nSP >= 6:
		
		spLayer = mamba_tree.spLayer
		spXPos = mamba_tree.spXPos
		spYPos = mamba_tree.spYPos
		spZPos = mamba_tree.spZPos
		spXErr = mamba_tree.spXErr
		spYErr = mamba_tree.spYErr
		spZErr = mamba_tree.spZErr
		spUPos = mamba_tree.spUPos
		spVPos = mamba_tree.spVPos
		spIsOnTrk = mamba_tree.spIsOnTrk

		# intersects
		nIntersect = mamba_tree.nIntersect
		intLayer = mamba_tree.intLayer
		intTrkID = mamba_tree.intTrkID
		intXPos = mamba_tree.intXPos
		intYPos = mamba_tree.intYPos
		intUPos = mamba_tree.intUPos
		intVPos = mamba_tree.intVPos
		intInside = mamba_tree.intInside
		
		spLayer = np.array(spLayer).tolist()
		nSP0 = spLayer.count(0)
		nSP1 = spLayer.count(1)
		nSP2 = spLayer.count(2)
		nSP3 = spLayer.count(3)
		nSP4 = spLayer.count(4)
		nSP5 = spLayer.count(5)
		nSP6 = spLayer.count(6)
		nSP7 = spLayer.count(7)
		
		spXPos = np.array(spXPos).tolist()
		spYPos = np.array(spYPos).tolist()
		
		X_ind = 4 * [False]
		Y_ind = 4 * [False]

		if nSP0 == 1 :
			Ymeas[k] = spYPos[spLayer.index(0)]
			Y_ind[k] = True
			k += 1
		if nSP1 == 1 :
			Ymeas[k] = spYPos[spLayer.index(1)]
			Y_ind[k] = True
			k += 1
		if nSP2 == 1 :
			Ymeas[k] = spYPos[spLayer.index(2)]
			Y_ind[k] = True
			k += 1
		if nSP3 == 1 :
			Ymeas[k] = spYPos[spLayer.index(3)]
			Y_ind[k] = True
			k += 1
		if nSP4 == 1 :
			Xmeas[l] = spXPos[spLayer.index(4)]
			X_ind[l] = True
			l += 1
		if nSP5 == 1 :
			Xmeas[l] = spXPos[spLayer.index(5)]
			X_ind[l] = True
			l += 1
		if nSP6 == 1 :
			Xmeas[l] = spXPos[spLayer.index(6)]
			X_ind[l] = True
			l += 1
		if nSP7 == 1 :
			Xmeas[l] = spXPos[spLayer.index(7)]
			X_ind[l] = True
			l += 1

		if k > 2 and l > 2 :
			
			good_track += 1
			if k > 3 and l > 3 : ultra_good_track += 1
			
			Zerrs_trak = np.array([Zerrs[el] for el in xrange(4) if X_ind[el]]).flatten('C')
			
			X_zvals_trak = np.array([X_zvals[el] for el in xrange(4) if X_ind[el]]).flatten('C')
			Xmeas_trak = np.array([Xmeas[el] for el in xrange(4) if X_ind[el]]).flatten('C')
			Xerrs_trak = np.array([Xerrs[el] for el in xrange(4) if X_ind[el]]).flatten('C')
			[resultX, resErrX] = my2DLineFit(l, X_zvals_trak, Xmeas_trak, Zerrs_trak, Xerrs_trak)
			
			Y_zvals_trak = np.array([Y_zvals[el] for el in xrange(4) if Y_ind[el]]).flatten('C')
			Ymeas_trak = np.array([Ymeas[el] for el in xrange(4) if Y_ind[el]]).flatten('C')
			Yerrs_trak = np.array([Yerrs[el] for el in xrange(4) if Y_ind[el]]).flatten('C')
			[resultY, resErrY] = my2DLineFit(k, Y_zvals_trak, Ymeas_trak, Zerrs_trak, Yerrs_trak)
			
			# print resultY
			for j in range(6):
				Xtrak[j] = resultX[1] * X_zvals[j] + resultX[0]
				Ytrak[j] = resultY[1] * Y_zvals[j] + resultY[0]
				# print Ytrak '''

	if ntrk == 1:

		#tracks
		good_track += 1
		chi2 = mamba_tree.chi2
		ndof = mamba_tree.ndof
		ItpX = mamba_tree.ItpX
		ItpY = mamba_tree.ItpY
		SlpX = mamba_tree.SlpX
		SlpY = mamba_tree.SlpY
		
		XCov[0][0] = mamba_tree.XCov00[0]
		XCov[0][1] = mamba_tree.XCov01[0]
		XCov[1][0] = mamba_tree.XCov10[0]
		XCov[1][1] = mamba_tree.XCov11[0]
		
		YCov[0][0] = mamba_tree.YCov00[0]
		YCov[0][1] = mamba_tree.YCov01[0]
		YCov[1][0] = mamba_tree.YCov10[0]
		YCov[1][1] = mamba_tree.YCov11[0]

		# spacepoints
		nSP = mamba_tree.nSP
		spLayer = mamba_tree.spLayer
		spXPos = mamba_tree.spXPos
		spYPos = mamba_tree.spYPos
		spUPos = mamba_tree.spUPos
		spVPos = mamba_tree.spVPos
		spIsOnTrk = mamba_tree.spIsOnTrk

		# intersects
		nIntersect = mamba_tree.nIntersect
		intLayer = mamba_tree.intLayer
		intTrkID = mamba_tree.intTrkID
		intXPos = mamba_tree.intXPos
		intYPos = mamba_tree.intYPos
		intUPos = mamba_tree.intUPos
		intVPos = mamba_tree.intVPos
		intInside = mamba_tree.intInside

		# print 'Number of single tracks : ', good_track

		k_x = ItpX[0]
		k_y = ItpY[0]

		# print 'k_x : %f, k_y : %f' % (k_x, k_y)

		for j in range(6):
			Xtrak[j] = SlpX[0] * X_zvals[j] + k_x
			Ytrak[j] = SlpY[0] * Y_zvals[j] + k_y
		
		varX_short = sqrt((X_zvals[4] ** 2) * XCov[1][1] + XCov[0][0] + X_zvals[4] * (XCov[1][0] + XCov[0][1]))
		varY_short = sqrt((Y_zvals[4] ** 2) * YCov[1][1] + YCov[0][0] + Y_zvals[4] * (YCov[1][0] + YCov[0][1]))
		varX_long = sqrt((X_zvals[5] ** 2) * XCov[1][1] + XCov[0][0] + X_zvals[5] * (XCov[1][0] + XCov[0][1]))
		varY_long = sqrt((Y_zvals[5] ** 2) * YCov[1][1] + YCov[0][0] + Y_zvals[5] * (YCov[1][0] + YCov[0][1]))
		
		residualY_short.Fill(varY_short)
		residualY_long.Fill(varY_long)
		residualX_short.Fill(varX_short)
		residualX_long.Fill(varX_long)

		# if good_track % 50000 == 0:
		if plotTracksFlag and good_track % 50000 == 0:
			nnn += 1
			XZtrakgraph =  TGraph(6,X_zvals,Xtrak)
			# XZtrakgraph.SetLineColor(2)
			# XZtrakgraph.SetLineStyle(1)
			YZtrakgraph =  TGraph(6,Y_zvals,Ytrak)
			YZtrakgraph.SetLineColor(2)
			YZtrakgraph.SetLineStyle(1)
			
			# Create the fitting function
   			fpol = TF1("fpol", "pol1", 0., 70.)
   			fpol.SetLineWidth(2)
   			XZtrakgraph.Fit(fpol, "QR")
			
			grintX = TGraphErrors(6)
			grintX.SetTitle("Fitted line with .95 conf. band")
			for i in xrange(6) :
				grintX.SetPoint(i, XZtrakgraph.GetX()[i], 0)
			# Compute the confidence intervals at the x points of the created graph
			(TVirtualFitter.GetFitter()).GetConfidenceIntervals(grintX)
			# Now the "grintX" graph contains function values as its y-coordinates
			# and confidence intervals as the errors on these coordinates
			# Draw the graph, the function and the confidence intervals
			
			
			c1 =  TCanvas("good_track_"+str(good_track),"good_track_"+str(good_track),600,400)
			c1.Divide(2,1)
			c1.cd(1)
			grintX.SetLineColor(kRed)
   			grintX.Draw("p")
   			XZtrakgraph.SetTitle("XZ projection of the track")
   			XZtrakgraph.GetXaxis().SetLimits(-50.0,100.0)
			XZtrakgraph.SetMinimum(-10.0)
			XZtrakgraph.SetMaximum(11.0)
			XZtrakgraph.GetYaxis().SetTitle("X (cm)")
			XZtrakgraph.GetXaxis().SetTitle("Z (cm)")
			XZtrakgraph.SetMarkerStyle(5)
			XZtrakgraph.SetMarkerSize(0.7)
   			XZtrakgraph.Draw("apsame")
			c1.cd(2)
			YZtrakgraph.Fit(fpol, "QR")
			grintY = TGraphErrors(6)
			grintY.SetTitle("Fitted line with .95 conf. band")
			for i in xrange(6) :
				grintY.SetPoint(i, YZtrakgraph.GetX()[i], 0)
			(TVirtualFitter.GetFitter()).GetConfidenceIntervals(grintY)
			grintY.SetLineColor(kRed)
   			grintY.Draw("p")
			YZtrakgraph.SetTitle("YZ projection of the track")
			YZtrakgraph.GetXaxis().SetLimits(-50.0,100.0)
			YZtrakgraph.SetMinimum(-10.0)
			YZtrakgraph.SetMaximum(11.0)
			YZtrakgraph.GetYaxis().SetTitle("Y (cm)")
			YZtrakgraph.GetXaxis().SetTitle("Z (cm)")
			YZtrakgraph.SetMarkerStyle(5)
			YZtrakgraph.SetMarkerSize(0.7)
			YZtrakgraph.Draw("PAsame")
			aname = "visual_" + `nnn` + ".png"
			c1.Modified()
			c1.SaveAs(plots_dir + aname)

		# loop on intersects
		for iIntersect in range(nIntersect):
			if intLayer[iIntersect] >= nDets:
				continue

			# loop on spacepoints
			for iSP in range(nSP):
				if spLayer[iSP] != intLayer[iIntersect]:
					continue

				# fill the histograms
				iDet = spLayer[iSP]
				resU = spUPos[iSP] - intUPos[iIntersect]

				resU_Det[iDet].Fill( resU )

				resUvsU_Det[iDet].Fill( intUPos[iIntersect], resU )
				resUvsV_Det[iDet].Fill( intVPos[iIntersect], resU )

				resUvsSlpX_Det[iDet].Fill( SlpX[ntrk-1], resU )
				resUvsSlpY_Det[iDet].Fill( SlpY[ntrk-1], resU )
			
				# resUvsSlpX_Det[iDet].Fill( resultX[1], resU )
				# resUvsSlpY_Det[iDet].Fill( resultY[1], resU )
		
		
		# DUT

		XYdut.Fill(Xtrak[4],Ytrak[4])
		XYdut_long.Fill(Xtrak[5],Ytrak[5])
		if (nt08==1) :
			XYtim.Fill(Xtrak[4],Ytrak[4])
			XYtim_long.Fill(Xtrak[5],Ytrak[5])
			TDC_tim_trk.Fill(conversion*tt08[0])
			N_tracks_tot += 1
			tot_tracks.Fill(Ytrak[4])
			tot_tracks_long.Fill(Ytrak[5])
		if nt08 == 2:
			TDC_hit_cor.Fill(conversion*(tt08[0] - tt08[1]))
			TDC_scin_2_hits.Fill(conversion*tt08[0],conversion*tt08[1])
			if abs(conversion*(tt08[0] - tt08[1])) < 23 :
				TDC_tim_trk_extra.Fill(conversion*(tt08[1]))
		if nt08 == 1 and nt03 == 1 and nt04 == 1:
			TDC_hit_cor_SiPM_up.Fill(conversion*(tt08[0] - (tt03[0] + tt04[0])/2))
			TDC_scin_2_hits_SiPM_up.Fill(conversion*tt08[0],conversion*(tt03[0] + tt04[0])/2)
		if nt08 == 1 and nt05 == 1 and nt06 == 1:
			TDC_hit_cor_SiPM_down.Fill(conversion*(tt08[0] - (tt05[0] + tt06[0])/2))
			TDC_scin_2_hits_SiPM_down.Fill(conversion*tt08[0],conversion*(tt05[0] + tt06[0])/2)
		if nt03 == 1 and nt04 == 1 and nt05 == 1 and nt06 == 1:
			TDC_hit_cor_SiPM_up_down.Fill(conversion*((tt03[0] + tt04[0])/2 - (tt05[0] + tt06[0])/2))
			TDC_scin_2_hits_SiPM_up_down.Fill(conversion*(tt03[0] + tt04[0])/2,conversion*(tt05[0] + tt06[0])/2)
		if noiseRedFlag :
			if (nt00>=1) :
				XYstr.Fill(Xtrak[4],Ytrak[4])
				for i in xrange(nt00) :
					TDC_str_trk.Fill(conversion*tt00[i])
						 
						 
			if ((nt00==1) and (nt08==1)) :
				XYboth.Fill(Xtrak[4],Ytrak[4])
				TDC_dif_trk.Fill(conversion*(tt00[0]-tt08[0]))
				VshapeX.Fill(Xtrak[4],conversion*(tt00[0]-tt08[0]))
				VshapeY.Fill(Ytrak[4],conversion*(tt00[0]-tt08[0]))
				very_good_track += 1
				N_tracks_straw += 1
				passed_tracks.Fill(Ytrak[4])
				if advVshapeFit:
#					resolutionS.Fill(Ytrak[4],inverseVshapeParametricMod(conversion*(tt00[0]-tt08[0]),parSinvTop,parSinvBot,SinvTopMin,SinvBotMin,Ytrak[4])-Ytrak[4])
					if conversion*(tt00[0]-tt08[0]) <= SinterTopMax or conversion*(tt00[0]-tt08[0]) <= SinterBotMax :
						if enableToyMC:
							if Ytrak[4] < bordS[1] and Ytrak[4] > bordS[0]:
								if enableUniform:
									Ytrak_uni = r.Uniform(bordS[0],bordS[1])
									Ytrak_toy = r.Gaus(Ytrak_uni,sigmaYtrak)
									time_toy = r.Gaus(vshapeParametricMod(Ytrak_uni,parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak_uni,parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
								else:
									Ytrak_toy = r.Gaus(Ytrak[4],sigmaYtrak)
									time_toy = r.Gaus(vshapeParametricMod(Ytrak[4],parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak[4],parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
								if time_toy <= SinterTopMax or time_toy <= SinterBotMax:
									VshapeY_toy.Fill(Ytrak_toy,time_toy)
									resolutionS_toy.Fill(Ytrak_toy,inverseVshapeParametricInter(time_toy,fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak_toy)-Ytrak_toy)
							resolutionS.Fill(Ytrak[4],inverseVshapeParametricInter(conversion*(tt00[0]-tt08[0]),fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak[4])-Ytrak[4])
						elif not enableToyMC:
							resolutionS.Fill(Ytrak[4],inverseVshapeParametricInter(conversion*(tt00[0]-tt08[0]),fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak[4])-Ytrak[4])
				else:
					resolutionS.Fill(Ytrak[4],inverseVshapeParametric(conversion*(tt00[0]-tt08[0]),parS,Ytrak[4])-Ytrak[4])
			   
			if ((nt00>1) and (nt08==1)) :
				if advVshapeFit:
					bestPoints = chooseTheBest(nt00,tt00,tt08,vshapeParametricMod(Ytrak[4],parS),distSigma)
				else:
					bestPoints = chooseTheBest(nt00,tt00,tt08,vshapeParametric(Ytrak[4],parS),distSigma)
				if len(bestPoints) < 1: continue
				if len(bestPoints) > 1: N_passed_2nd_hits += len(bestPoints) - 1
				for best in bestPoints:
					XYboth.Fill(Xtrak[4],Ytrak[4])
					TDC_dif_trk.Fill(conversion*(tt00[best]-tt08[0]))
					VshapeX.Fill(Xtrak[4],conversion*(tt00[best]-tt08[0]))
					VshapeY.Fill(Ytrak[4],conversion*(tt00[best]-tt08[0]))
					very_good_track += 1
					N_tracks_straw += 1
					passed_tracks.Fill(Ytrak[4])
					if advVshapeFit:
	#					resolutionS.Fill(Ytrak[4],inverseVshapeParametricMod(conversion*(tt00[best]-tt08[0]),parSinvTop,parSinvBot,SinvTopMin,SinvBotMin,Ytrak[4])-Ytrak[4])
						if conversion*(tt00[best]-tt08[0]) <= SinterTopMax or conversion*(tt00[best]-tt08[0]) <= SinterBotMax :
							if enableToyMC:
								if Ytrak[4] < bordS[1] and Ytrak[4] > bordS[0]:
									if enableUniform:
										Ytrak_uni = r.Uniform(bordS[0],bordS[1])
										Ytrak_toy = r.Gaus(Ytrak_uni,sigmaYtrak)
										time_toy = r.Gaus(vshapeParametricMod(Ytrak_uni,parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak_uni,parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
									else:
										Ytrak_toy = r.Gaus(Ytrak[4],sigmaYtrak)
										time_toy = r.Gaus(vshapeParametricMod(Ytrak[4],parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak[4],parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
									if time_toy <= SinterTopMax or time_toy <= SinterBotMax:
										VshapeY_toy.Fill(Ytrak_toy,time_toy)
										resolutionS_toy.Fill(Ytrak_toy,inverseVshapeParametricInter(time_toy,fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak_toy)-Ytrak_toy)
								resolutionS.Fill(Ytrak[4],inverseVshapeParametricInter(conversion*(tt00[best]-tt08[0]),fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak[4])-Ytrak[4])
							elif not enableToyMC:
								resolutionS.Fill(Ytrak[4],inverseVshapeParametricInter(conversion*(tt00[best]-tt08[0]),fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak[4])-Ytrak[4])
					else:
						resolutionS.Fill(Ytrak[4],inverseVshapeParametric(conversion*(tt00[best]-tt08[0]),parS,Ytrak[4])-Ytrak[4])

			if (nt01>=1) :
				XYstr_long.Fill(Xtrak[5],Ytrak[5])
				for i in xrange(nt01) :
					TDC_str_trk_long.Fill(conversion*tt01[i])


			if ((nt01==1) and (nt08==1)) :
				XYboth_long.Fill(Xtrak[5],Ytrak[5])
				TDC_dif_trk_long.Fill(conversion*(tt01[0]-tt08[0]))
				VshapeX_long.Fill(Xtrak[5],conversion*(tt01[0]-tt08[0]))
				VshapeY_long.Fill(Ytrak[5],conversion*(tt01[0]-tt08[0]))
				very_good_track_long += 1
				N_tracks_straw_long += 1
				passed_tracks_long.Fill(Ytrak[5])
				if advVshapeFit:
#					resolutionL.Fill(Ytrak[5],inverseVshapeParametricMod(conversion*(tt01[0]-tt08[0]),parLinvTop,parLinvBot,LinvTopMin,LinvBotMin,Ytrak[5])-Ytrak[5])
					if conversion*(tt01[0]-tt08[0]) <= LinterTopMax or conversion*(tt01[0]-tt08[0]) <= LinterBotMax :
						if enableToyMC:
							if Ytrak[5] < bordL[1] and Ytrak[5] > bordL[0]:
								if enableUniform:
									Ytrak_uni_long = r.Uniform(bordL[0],bordL[1])
									Ytrak_toy_long = r.Gaus(Ytrak_uni_long,sigmaYtrak)
									time_toy_long = r.Gaus(vshapeParametricMod(Ytrak_uni_long,parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak_uni_long,parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
								else:
									Ytrak_toy_long = r.Gaus(Ytrak[5],sigmaYtrak)
									time_toy_long = r.Gaus(vshapeParametricMod(Ytrak[5],parL),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak[5],parL)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
								if time_toy_long <= LinterTopMax or time_toy_long <= LinterBotMax:
									VshapeY_long_toy.Fill(Ytrak_toy_long,time_toy_long)
									resolutionL_toy.Fill(Ytrak_toy_long,inverseVshapeParametricInter(time_toy_long,fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak_toy_long)-Ytrak_toy_long)
							resolutionL.Fill(Ytrak[5],inverseVshapeParametricInter(conversion*(tt01[0]-tt08[0]),fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak[5])-Ytrak[5])
						elif not enableToyMC:
							resolutionL.Fill(Ytrak[5],inverseVshapeParametricInter(conversion*(tt01[0]-tt08[0]),fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak[5])-Ytrak[5])
				else:
					resolutionL.Fill(Ytrak[5],inverseVshapeParametric(conversion*(tt01[0]-tt08[0]),parL,Ytrak[5])-Ytrak[5])

			if ((nt01>1) and (nt08==1)) :
				if advVshapeFit:
					bestPoints = chooseTheBest(nt01,tt01,tt08,vshapeParametricMod(Ytrak[5],parL),distSigma)
				else:
					bestPoints = chooseTheBest(nt01,tt01,tt08,vshapeParametric(Ytrak[5],parL),distSigma)
				if len(bestPoints) < 1: continue
				for best in bestPoints:
					XYboth_long.Fill(Xtrak[5],Ytrak[5])
					TDC_dif_trk_long.Fill(conversion*(tt01[best]-tt08[0]))
					VshapeX_long.Fill(Xtrak[5],conversion*(tt01[best]-tt08[0]))
					VshapeY_long.Fill(Ytrak[5],conversion*(tt01[best]-tt08[0]))
					very_good_track_long += 1
					N_tracks_straw_long += 1
					passed_tracks_long.Fill(Ytrak[5])
					if advVshapeFit:
	#					resolutionL.Fill(Ytrak[5],inverseVshapeParametricMod(conversion*(tt01[best]-tt08[0]),parLinvTop,parLinvBot,LinvTopMin,LinvBotMin,Ytrak[5])-Ytrak[5])
						if conversion*(tt01[best]-tt08[0]) <= LinterTopMax or conversion*(tt01[best]-tt08[0]) <= LinterBotMax :
							if enableToyMC:
								if Ytrak[5] < bordL[1] and Ytrak[5] > bordL[0]:
									if enableUniform:
										Ytrak_uni_long = r.Uniform(bordL[0],bordL[1])
										Ytrak_toy_long = r.Gaus(Ytrak_uni_long,sigmaYtrak)
										time_toy_long = r.Gaus(vshapeParametricMod(Ytrak_uni_long,parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak_uni_long,parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
									else:
										Ytrak_toy_long = r.Gaus(Ytrak[5],sigmaYtrak)
										time_toy_long = r.Gaus(vshapeParametricMod(Ytrak[5],parL),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak[5],parL)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
									if time_toy_long <= LinterTopMax or time_toy_long <= LinterBotMax:
										VshapeY_long_toy.Fill(Ytrak_toy_long,time_toy_long)
										resolutionL_toy.Fill(Ytrak_toy_long,inverseVshapeParametricInter(time_toy_long,fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak_toy_long)-Ytrak_toy_long)
								resolutionL.Fill(Ytrak[5],inverseVshapeParametricInter(conversion*(tt01[best]-tt08[0]),fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak[5])-Ytrak[5])
							elif not enableToyMC:
								resolutionL.Fill(Ytrak[5],inverseVshapeParametricInter(conversion*(tt01[best]-tt08[0]),fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak[5])-Ytrak[5])
					else:
						resolutionL.Fill(Ytrak[5],inverseVshapeParametric(conversion*(tt01[best]-tt08[0]),parL,Ytrak[5])-Ytrak[5])
		else :
			if nt00 == 1:
				XYstr.Fill(Xtrak[4],Ytrak[4])
				for i in range(nt00):
					TDC_str_trk.Fill(conversion*tt00[i])

			if nt01 == 1:
				XYstr_long.Fill(Xtrak[5],Ytrak[5])
				TDC_str_trk_long.Fill(conversion*tt01[0])

			if (nt00 == 1 and nt08 == 1):
				XYboth.Fill(Xtrak[4],Ytrak[4])
				for i in range(nt00):
					TDC_dif_trk.Fill(conversion*(tt00[i]-tt08[0]))
					VshapeX.Fill(Xtrak[4],conversion*(tt00[i]-tt08[0]))
					VshapeY.Fill(Ytrak[4],conversion*(tt00[i]-tt08[0]))
					N_tracks_straw += 1
				very_good_track += 1
				passed_tracks.Fill(Ytrak[4])

			if (nt00 == 1 and nt03 == 1 and nt04 == 1):
				VshapeY_SiPM_up.Fill(Ytrak[4],conversion*(tt00[0]-(tt03[0] + tt04[0])/2))

			if (nt00 == 1 and nt05 == 1 and nt06 == 1):
				VshapeY_SiPM_down.Fill(Ytrak[4],conversion*(tt00[0]-(tt05[0] + tt06[0])/2))

			if (nt01 == 1 and nt08 == 1):
				XYboth_long.Fill(Xtrak[5],Ytrak[5])
				TDC_dif_trk_long.Fill(conversion*(tt01[0]-tt08[0]))
				VshapeX_long.Fill(Xtrak[5],conversion*(tt01[0]-tt08[0]))
				VshapeY_long.Fill(Ytrak[5],conversion*(tt01[0]-tt08[0]))
				very_good_track_long += 1
				N_tracks_straw_long += 1
				passed_tracks_long.Fill(Ytrak[5])
				
			if (nt01 == 1 and nt03 == 1 and nt04 == 1):
				VshapeY_long_SiPM_up.Fill(Ytrak[5],conversion*(tt01[0]-(tt03[0] + tt04[0])/2))
				
			if (nt01 == 1 and nt05 == 1 and nt06 == 1):
				VshapeY_long_SiPM_down.Fill(Ytrak[5],conversion*(tt01[0]-(tt05[0] + tt06[0])/2))

	if nt00 == 1:
		TDC_str_raw.Fill(conversion*tt00[0])

	if nt01 == 1:
		TDC_str_raw_long.Fill(conversion*tt01[0])

	if nt08 == 1:
		TDC_tim_raw.Fill(conversion*tt08[0])
		TDC_tim_raw_long.Fill(conversion*tt08[0])

	if (nt00 == 1 and nt08 == 1):
		for i in range(nt00):
			TDC_dif_raw.Fill(conversion*(tt00[i]-tt08[0]))

	if (nt01 == 1 and nt08 == 1):
		TDC_dif_raw_long.Fill(conversion*(tt01[0]-tt08[0]))

	count_processed +=1

print  "   Events processed: ", count_processed
print  "   Events with a single track: ", good_track
print  "   Events with a single track and drift time in a short tube: ", very_good_track
print  "   Events with a single track and drift time in a long tube: ", very_good_track_long
if not loadSyncFileFlag:
	print  "   timeStamp_i = ",   timestamp_iEvt
	print  "   time_ms_i = ",  time_ms_i
	print  "   time_us_i = ",  time_us_i
	print  "   pcclock_i   = ",  PCclock_i
	print  "   TDCstamp_i  = ",  TDCstamp_i
print  "   Efficiency of the short straw = ", N_tracks_straw / N_tracks_tot
print  "   Efficiency of the long straw = ", N_tracks_straw_long / N_tracks_tot
print  "   Probability to get a hit between a correct hit and the extrapolated value = ", N_passed_2nd_hits / N_tracks_straw

if not noiseRedFlag :
	residualY_short_RMS = residualY_short.GetMean()
	residualY_long_RMS = residualY_long.GetMean()
	tracks_RMS = [residualY_short_RMS,residualY_long_RMS]
	if not advVshapeFit:
		if referenceRunFlag:
			improveNoise("outputfiles/", vshf, tracks_RMS, [parRefS[0],parRefL[0]])
		else:
			improveNoise("outputfiles/", vshf, tracks_RMS)
	else:
		if not binSizeFlag:
			improveNoiseMod("outputfiles/", vshf, sys.argv[3], tracks_RMS)
		else:
			improveNoiseMod("outputfiles/", vshf, sys.argv[3], tracks_RMS, sys.argv[4], sys.argv[5])

# Calculate straw resolution

if noiseRedFlag :
	resolution_slices_toy = TObjArray()
	resolution_slices = TObjArray()
	resolution_slices_long = TObjArray()
	resolution_slices_long_toy = TObjArray()
	
	gf = TF1("gf", "gaus", -0.3, 0.3)
	gf.SetParLimits(0,5,5000)
	gf.SetParLimits(1,-0.2,0.2)
	gf.SetParLimits(2,0.001,0.1)
	gf.SetParameter(0,500.)
	gf.SetParameter(1,0.0001)
	gf.SetParameter(2,0.01)
	if enableToyMC:
		resolutionS.FitSlicesY(gf,resolutionS.GetXaxis().FindBin(bordS[0]),resolutionS.GetXaxis().FindBin(bordS[1]),16,"QNR",resolution_slices) # 16 18 22
		resolutionS_toy.FitSlicesY(gf,resolutionS_toy.GetXaxis().FindBin(bordS[0]),resolutionS_toy.GetXaxis().FindBin(bordS[1]),4,"QNR",resolution_slices_toy) # 16 18 22
	else:
		resolutionS.FitSlicesY(gf,resolutionS.GetXaxis().FindBin(bordS[0]),resolutionS.GetXaxis().FindBin(bordS[1]),16,"QNR",resolution_slices) # 16 18 22
	#resolution_slices = FindRMSofSlices(resolutionS,False,resolutionS.GetXaxis().FindBin(bordS[0]),resolutionS.GetXaxis().FindBin(bordS[1]),18,"QNR") # 16 18 22
	gf.SetParameter(0,500.)
	gf.SetParameter(1,0.0001)
	gf.SetParameter(2,0.01)
	if enableToyMC:
		resolutionL.FitSlicesY(gf,resolutionL.GetXaxis().FindBin(bordL[0]),resolutionL.GetXaxis().FindBin(bordL[1]),18,"QNR",resolution_slices_long) # 18 20
		resolutionL_toy.FitSlicesY(gf,resolutionL_toy.GetXaxis().FindBin(bordL[0]),resolutionL_toy.GetXaxis().FindBin(bordL[1]),4,"QNR",resolution_slices_long_toy) # 18 20
	else:
		resolutionL.FitSlicesY(gf,resolutionL.GetXaxis().FindBin(bordL[0]),resolutionL.GetXaxis().FindBin(bordL[1]),18,"QNR",resolution_slices_long) # 18 20
	#resolution_slices_long = FindRMSofSlices(resolutionL,False,resolutionL.GetXaxis().FindBin(bordL[0]),resolutionL.GetXaxis().FindBin(bordL[1]),18,"QNR") # 18 20

'''VshapeY_first = VshapeY.Clone()
VshapeY_second = VshapeY.Clone()
(VshapeY_first.GetXaxis()).SetRangeUser(-1.5, 0.)
(VshapeY_second.GetXaxis()).SetRangeUser(0.,1.5)

VshapeY_long_first = VshapeY_long.Clone()
VshapeY_long_second = VshapeY_long.Clone()
(VshapeY_long_first.GetXaxis()).SetRangeUser(-1.5, 0.)
(VshapeY_long_second.GetXaxis()).SetRangeUser(0.,1.5)

gf_first = TF1("gf_first", "gaus", -1.5, 0.3)
gf_second = TF1("gf_second", "gaus", -0.3, 1.5)

aSlicesS_first = TObjArray()
aSlicesS_second = TObjArray()
vshape_pars_histo = TObjArray()
 VshapeY_first.FitSlicesX(gf_first,(VshapeY_first.GetYaxis()).FindBin(-100),(VshapeY_first.GetYaxis()).FindBin(1200),10,"QNR",aSlicesS_first)
 VshapeY_second.FitSlicesX(gf_second,(VshapeY_second.GetYaxis()).FindBin(-100),(VshapeY_second.GetYaxis()).FindBin(1200),10,"QNR",aSlicesS_second)
VshapeY.GetXaxis().SetRangeUser(-1.5,1.5)
 VshapeY.FitSlicesX(gf_two,VshapeY.GetYaxis().FindBin(-100),VshapeY.GetYaxis().FindBin(1200),7,"R",vshape_pars_histo)
vshape_pars_histo = FitSlicesDoubleGaussian(VshapeY,True,VshapeY.GetYaxis().FindBin(-100),VshapeY.GetYaxis().FindBin(1200),10,"Rg2",vshape_pars_histo)

aSlicesL_first = TObjArray()
aSlicesL_second = TObjArray()
vshape_pars_histo_long = TObjArray()
 VshapeY_long_first.FitSlicesX(gf_first,(VshapeY_long_first.GetYaxis()).FindBin(-100),(VshapeY_long_first.GetYaxis()).FindBin(400),20,"QNR",aSlicesL_first)
 VshapeY_long_second.FitSlicesX(gf_second,(VshapeY_long_second.GetYaxis()).FindBin(-100),(VshapeY_long_second.GetYaxis()).FindBin(400),20,"QNR",aSlicesL_second)
VshapeY_long.GetXaxis().SetRangeUser(-1.5,1.5)
 VshapeY_long.FitSlicesX(gf_two,VshapeY_long.GetYaxis().FindBin(-100),VshapeY_long.GetYaxis().FindBin(400),7,"R",vshape_pars_histo_long)
vshape_pars_histo_long = FitSlicesDoubleGaussian(VshapeY_long,True,VshapeY_long.GetYaxis().FindBin(-100),VshapeY_long.GetYaxis().FindBin(400),10,"Rg2",vshape_pars_histo_long)'''

# plotting

if not batchFlag :

	gStyle.SetOptStat(111111)

	# timestamp plots

	C_cl =  TCanvas("C_cl", "C_cl",1280,800)
	C_cl.Divide(2,2)
	C_cl.cd(1)
	Pclock.Draw()
	C_cl.cd(2)
	TDCstamp.Draw()
	C_cl.cd(3)
	T_ms.Draw()
	C_cl.cd(4)
	T_us.Draw()
	C_cl.SaveAs(plots_dir + "cl.png")

	C_cldif =  TCanvas("C_cldif", "C_cldif",1280,800)
	C_cldif.Divide(2,3)
	C_cldif.cd(1)
	Clocks.Draw()
	C_cldif.cd(2)
	Clockt.Draw()
	C_cldif.cd(3)
	TSdiff.Draw()
	C_cldif.cd(4)
	TSDiff.Draw()
	C_cldif.cd(5)
	TSdif2.Draw()
	C_cldif.cd(6)
	TSDif2.Draw()
	C_cldif.SaveAs(plots_dir + "cldif.png")

	# DUT, TDC and Vshape plots
	statXdef = gStyle.GetStatX()
	statYdef = gStyle.GetStatY()
	gStyle.SetStatX(0.9)
	gStyle.SetStatY(0.9)
	C_XYdut =  TCanvas("C_XYdut", "C_XYdut",1280,800)
	C_XYdut.Divide(2,2)
	C_XYdut.cd(1)
	XYdut.SetXTitle("X (cm)")
	XYdut.SetYTitle("Y (cm)")
	XYdut.Draw("colz")
	C_XYdut.cd(2)
	XYstr.SetXTitle("X (cm)")
	XYstr.SetYTitle("Y (cm)")
	XYstr.Draw("colz")
	C_XYdut.cd(3)
	XYtim.SetXTitle("X (cm)")
	XYtim.SetYTitle("Y (cm)")
	XYtim.Draw("colz")
	C_XYdut.cd(4)
	XYboth.SetXTitle("X (cm)")
	XYboth.SetYTitle("Y (cm)")
	XYboth.Draw("colz")
	C_XYdut.SaveAs(plots_dir + "XYdut.png")

	C_XYdut_long =  TCanvas("C_XYdut_long", "C_XYdut_long",1280,800)
	C_XYdut_long.Divide(2,2)
	C_XYdut_long.cd(1)
	XYdut_long.SetXTitle("X (cm)")
	XYdut_long.SetYTitle("Y (cm)")
	XYdut_long.Draw("colz")
	C_XYdut_long.cd(2)
	XYstr_long.SetXTitle("X (cm)")
	XYstr_long.SetYTitle("Y (cm)")
	XYstr_long.Draw("colz")
	C_XYdut_long.cd(3)
	XYtim_long.SetXTitle("X (cm)")
	XYtim_long.SetYTitle("Y (cm)")
	XYtim_long.Draw("colz")
	C_XYdut_long.cd(4)
	XYboth_long.SetXTitle("X (cm)")
	XYboth_long.SetYTitle("Y (cm)")
	XYboth_long.Draw("colz")
	C_XYdut_long.SaveAs(plots_dir + "XYdut_long.png")
	gStyle.SetStatX(statXdef)
	gStyle.SetStatY(statYdef)

	C_tdc =  TCanvas("C_tdc", "C_tdc",1280,800)
	C_tdc.Divide(3,2)
	C_tdc.cd(1)
	TDC_str_raw.SetXTitle("TDC (ns)")
	TDC_str_raw.Draw()
	C_tdc.cd(2)
	TDC_tim_raw.SetXTitle("TDC (ns)")
	TDC_tim_raw.Draw()
	C_tdc.cd(3)
	TDC_dif_raw.SetXTitle("TDC (ns)")
	TDC_dif_raw.Draw()
	C_tdc.cd(4)
	TDC_str_trk.SetXTitle("TDC (ns)")
	TDC_str_trk.Draw()
	C_tdc.cd(5)
	TDC_tim_trk.SetXTitle("TDC (ns)")
	TDC_tim_trk.Draw()
	TDC_tim_trk_extra.SetLineColor(kRed)
	TDC_tim_trk_extra.Draw("same")
	C_tdc.cd(6)
	TDC_dif_trk.SetXTitle("TDC (ns)")
	TDC_dif_trk.Draw()
	C_tdc.SaveAs(plots_dir + "tdc.png")

	C_tdc_long =  TCanvas("C_tdc_long", "C_tdc_long",1280,800)
	C_tdc_long.Divide(3,2)
	C_tdc_long.cd(1)
	TDC_str_raw_long.SetXTitle("TDC (ns)")
	TDC_str_raw_long.Draw()
	C_tdc_long.cd(2)
	TDC_tim_raw_long.SetXTitle("TDC (ns)")
	TDC_tim_raw_long.Draw()
	C_tdc_long.cd(3)
	TDC_dif_raw_long.SetXTitle("TDC (ns)")
	TDC_dif_raw_long.Draw()
	C_tdc_long.cd(4)
	TDC_str_trk_long.SetXTitle("TDC (ns)")
	TDC_str_trk_long.Draw()
	C_tdc_long.cd(5)
	TDC_tim_trk.SetXTitle("TDC (ns)")
	TDC_tim_trk.Draw()
	TDC_tim_trk_extra.SetLineColor(kRed)
	TDC_tim_trk_extra.Draw("same")
	C_tdc_long.cd(6)
	TDC_dif_trk_long.SetXTitle("TDC (ns)")
	TDC_dif_trk_long.Draw()
	C_tdc_long.SaveAs(plots_dir + "tdc_long.png")

	C_t =  TCanvas("C_t", "C_t",1280,800)
	C_t.Divide(2,2)
	C_t.cd(1)
	VshapeX.SetXTitle("X (cm)")
	VshapeX.SetYTitle("t (ns)")
	VshapeX.Draw("colz")
	gPad.Update()
	stVshapeX = VshapeX.GetListOfFunctions().FindObject("stats").Clone("stVshapeX")
	stVshapeX.SetX1NDC(0.7)
	stVshapeX.SetX2NDC(0.9)
	stVshapeX.SetY1NDC(0.7)
	stVshapeX.SetY2NDC(0.9)
	gPad.Update()
	C_t.cd(2)
	VshapeY.SetXTitle("Y (cm)")
	VshapeY.SetYTitle("t (ns)")
	VshapeY.Draw("colz")
	gPad.Update()
	stVshapeY = VshapeY.GetListOfFunctions().FindObject("stats").Clone("stVshapeY")
	stVshapeY.SetX1NDC(0.4)
	stVshapeY.SetX2NDC(0.7)
	stVshapeY.SetY1NDC(0.7)
	stVshapeY.SetY2NDC(0.9)
	gPad.Update()
	if (noiseRedFlag) :
		if advVshapeFit:
			VshapeY_fit = TF1("VshapeY_fit", vshapeParametricMod, bordS[0], bordS[1], 7)
		else:
			VshapeY_fit = TF1("VshapeY_fit", vshapeParametric, bordS[0], bordS[1], 3)
		parS = np.array(parS).flatten('C')
		VshapeY_fit.SetParameters(parS)
		VshapeY_fit.Draw("same")
	C_t.cd(3)
	VshapeX_long.SetXTitle("X (cm)")
	VshapeX_long.SetYTitle("t (ns)")
	VshapeX_long.Draw("colz")
	gPad.Update()
	stVshapeX_long = VshapeX_long.GetListOfFunctions().FindObject("stats").Clone("stVshapeX_long")
	stVshapeX_long.SetX1NDC(0.7)
	stVshapeX_long.SetX2NDC(0.9)
	stVshapeX_long.SetY1NDC(0.7)
	stVshapeX_long.SetY2NDC(0.9)
	gPad.Update()
	C_t.cd(4)
	VshapeY_long.SetXTitle("Y (cm)")
	VshapeY_long.SetYTitle("t (ns)")
	VshapeY_long.Draw("colz")
	gPad.Update()
	stVshapeY_long = VshapeY_long.GetListOfFunctions().FindObject("stats").Clone("stVshapeY_long")
	stVshapeY_long.SetX1NDC(0.4)
	stVshapeY_long.SetX2NDC(0.7)
	stVshapeY_long.SetY1NDC(0.7)
	stVshapeY_long.SetY2NDC(0.9)
	gPad.Update()
	if (noiseRedFlag) :
		if advVshapeFit:
			VshapeY_fit_long = TF1("VshapeY_fit_long", vshapeParametricMod, bordL[0], bordL[1], 7)
		else:
			VshapeY_fit_long = TF1("VshapeY_fit_long", vshapeParametric, bordL[0], bordL[1], 3)
		parL = np.array(parL).flatten('C')
		VshapeY_fit_long.SetParameters(parL)
		VshapeY_fit_long.Draw("same")
	C_t.SaveAs(plots_dir + "vshapes.png")
	
	gStyle.SetStatX(0.6)
	gStyle.SetStatY(0.9)
	C_t_SiPMs = TCanvas("C_t_SiPMs", "C_t_SiPMs", 1280, 800)
	C_t_SiPMs.Divide(2,2)
	C_t_SiPMs.cd(1)
	VshapeY_SiPM_up.SetXTitle("Y (cm)")
	VshapeY_SiPM_up.SetYTitle("t (ns)")
	VshapeY_SiPM_up.Draw("colz")
	C_t_SiPMs.cd(2)
	VshapeY_SiPM_down.SetXTitle("Y (cm)")
	VshapeY_SiPM_down.SetYTitle("t (ns)")
	VshapeY_SiPM_down.Draw("colz")
	C_t_SiPMs.cd(3)
	VshapeY_long_SiPM_up.SetXTitle("Y (cm)")
	VshapeY_long_SiPM_up.SetYTitle("t (ns)")
	VshapeY_long_SiPM_up.Draw("colz")
	C_t_SiPMs.cd(4)
	VshapeY_long_SiPM_down.SetXTitle("Y (cm)")
	VshapeY_long_SiPM_down.SetYTitle("t (ns)")
	VshapeY_long_SiPM_down.Draw("colz")
	C_t_SiPMs.SaveAs(plots_dir + "vshapes_SiPMs.png")
	gStyle.SetStatX(statXdef)
	gStyle.SetStatY(statYdef)

	C_eff = TCanvas("C_eff", "C_eff", 1280, 800)
	C_eff.Divide(2,2)
	C_eff.cd(1)
	effY = TGraphAsymmErrors()
	effY.Divide(passed_tracks, tot_tracks, "cl=0.683 b(1,1) mode")
	effY.SetTitle("Efficiency_short_tube;Y (cm);#epsilon")
	# effY.SetYTitle("eff")
	effY.Draw()
	C_eff.cd(2)
	effY_zoom = effY.Clone()
	effY_zoom.GetYaxis().SetRangeUser(0.9,1.1)
	effY_zoom.Draw()
	C_eff.cd(3)
	effY_long = TGraphAsymmErrors()
	effY_long.Divide(passed_tracks_long, tot_tracks_long, "cl=0.683 b(1,1) mode")
	effY_long.SetTitle("Efficiency_long_tube;Y (cm);#epsilon")
	# effY_long.SetYTitle("eff")
	effY_long.Draw()
	C_eff.cd(4)
	effY_long_zoom = effY_long.Clone()
	effY_long_zoom.GetYaxis().SetRangeUser(0.9,1.1)
	effY_long_zoom.Draw()
	C_eff.SaveAs(plots_dir + "efficiency.png")

	if noiseRedFlag :
		C_resolSL = TCanvas("C_resolSL", "C_resolSL", 1280, 800)
		C_resolSL.Divide(5,2)
		C_resolSL.cd(1)
		if enableToyMC:
			resolutionS_toy.SetXTitle("Y (cm)")
			resolutionS_toy.SetYTitle("residual Y (cm)")
			resolutionS_toy.Draw("colz")
			gPad.Update()
			stResolutionS_toy = resolutionS_toy.GetListOfFunctions().FindObject("stats").Clone("stResolutionS")
			stResolutionS_toy.SetX1NDC(0.7)
			stResolutionS_toy.SetX2NDC(0.9)
			stResolutionS_toy.SetY1NDC(0.7)
			stResolutionS_toy.SetY2NDC(0.9)
			gPad.Update()
		else:
			resolutionS.SetXTitle("Y (cm)")
			resolutionS.SetYTitle("residual Y (cm)")
			resolutionS.Draw("colz")
			gPad.Update()
			stResolutionS = resolutionS.GetListOfFunctions().FindObject("stats").Clone("stResolutionS")
			stResolutionS.SetX1NDC(0.7)
			stResolutionS.SetX2NDC(0.9)
			stResolutionS.SetY1NDC(0.7)
			stResolutionS.SetY2NDC(0.9)
			gPad.Update()
		#VshapeY.SetXTitle("Y (cm)")
		#VshapeY.SetYTitle("t (ns)")
		#VshapeY.Draw("colz")
		C_resolSL.cd(2)
		resolution_slices[0].SetXTitle("Y (cm)")
		resolution_slices[0].Draw()
		#vshape_pars_histo[1].Draw()
		#meanS_first = aSlicesS_first[1].Clone()
		#meanS_first.Draw()
		C_resolSL.cd(3)
		resolution_slices[1].SetXTitle("Y (cm)")
		resolution_slices[1].SetYTitle("residual Y (cm)")
		resolution_slices[1].Draw()
		#vshape_pars_histo[2].Draw()
		#rmsS_first = aSlicesS_first[2].Clone()
		#rmsS_first.Draw()
		C_resolSL.cd(4)
		resolution_slices[2].SetMaximum(0.05)
		resolution_slices[2].SetMinimum(0.0)
		resolution_slices[2].SetXTitle("Y (cm)")
		resolution_slices[2].SetYTitle("resolution (cm)")
		resolution_slices[2].Draw()
		#vshape_pars_histo[4].Draw()
		#meanS_second = aSlicesS_second[1].Clone()
		#meanS_second.Draw()
		C_resolSL.cd(5)
		resolution_slices[3].SetXTitle("Y (cm)")
		resolution_slices[3].Draw()
		#vshape_pars_histo[5].Draw()
		#rmsS_second = aSlicesS_second[2].Clone()
		#rmsS_second.Draw()
		#raw_input("Press any key to continue\n")
		C_resolSL.cd(6)
		if enableToyMC:
			resolutionL_toy.SetXTitle("Y (cm)")
			resolutionL_toy.SetYTitle("residual Y (cm)")
			resolutionL_toy.Draw("colz")
			gPad.Update()
			stResolutionL_toy = resolutionL_toy.GetListOfFunctions().FindObject("stats").Clone("stResolutionL")
			stResolutionL_toy.SetX1NDC(0.7)
			stResolutionL_toy.SetX2NDC(0.9)
			stResolutionL_toy.SetY1NDC(0.7)
			stResolutionL_toy.SetY2NDC(0.9)
			gPad.Update()
		else:
			resolutionL.SetXTitle("Y (cm)")
			resolutionL.SetYTitle("residual Y (cm)")
			resolutionL.Draw("colz")
			gPad.Update()
			stResolutionL = resolutionL.GetListOfFunctions().FindObject("stats").Clone("stResolutionL")
			stResolutionL.SetX1NDC(0.7)
			stResolutionL.SetX2NDC(0.9)
			stResolutionL.SetY1NDC(0.7)
			stResolutionL.SetY2NDC(0.9)
			gPad.Update()
		#VshapeY_long.SetXTitle("Y (cm)")
		#VshapeY_long.SetYTitle("t (ns)")
		#VshapeY_long.Draw("colz")
		C_resolSL.cd(7)
		resolution_slices_long[0].SetXTitle("Y (cm)")
		resolution_slices_long[0].Draw()
		#vshape_pars_histo_long[1].Draw()
		#meanL_first = aSlicesL_first[0].Clone()
		#meanL_first.Draw()
		C_resolSL.cd(8)
		resolution_slices_long[1].SetXTitle("Y (cm)")
		resolution_slices_long[1].SetYTitle("residual Y (cm)")
		resolution_slices_long[1].Draw()
		#vshape_pars_histo_long[2].Draw()
		#rmsL_first = aSlicesL_first[2].Clone()
		#rmsL_first.Draw()
		C_resolSL.cd(9)
		resolution_slices_long[2].SetMaximum(0.05)
		resolution_slices_long[2].SetMinimum(0.0)
		resolution_slices_long[2].SetXTitle("Y (cm)")
		resolution_slices_long[2].SetYTitle("resolution (cm)")
		resolution_slices_long[2].Draw()
		#vshape_pars_histo_long[4].Draw()
		#meanL_second = aSlicesL_second[0].Clone()
		#meanL_second.Draw()
		C_resolSL.cd(10)
		resolution_slices_long[3].SetXTitle("Y (cm)")
		resolution_slices_long[3].Draw()
		#vshape_pars_histo_long[5].Draw()
		#rmsL_second = aSlicesL_second[2].Clone()
		#rmsL_second.Draw()
		C_resolSL.SaveAs(plots_dir + "resolutionS+L.png")

		if enableToyMC:
			C_toyMC =  TCanvas("C_toyMC", "C_toyMC",1280,800)
			C_toyMC.Divide(1,2)
			C_toyMC.cd(1)
			VshapeY_toy.SetXTitle("Y (cm)")
			VshapeY_toy.SetYTitle("t (ns)")
			VshapeY_toy.Draw("colz")
			gPad.Update()
			stVshapeY_toy = VshapeY_toy.GetListOfFunctions().FindObject("stats").Clone("stVshapeY_toy")
			stVshapeY_toy.SetX1NDC(0.4)
			stVshapeY_toy.SetX2NDC(0.7)
			stVshapeY_toy.SetY1NDC(0.7)
			stVshapeY_toy.SetY2NDC(0.9)
			gPad.Update()
			if advVshapeFit:
				VshapeY_fit = TF1("VshapeY_fit", vshapeParametricMod, bordS[0], bordS[1], 7)
			else:
				VshapeY_fit = TF1("VshapeY_fit", vshapeParametric, bordS[0], bordS[1], 3)
			parS = np.array(parS).flatten('C')
			VshapeY_fit.SetParameters(parS)
			VshapeY_fit.Draw("same")
			C_toyMC.cd(2)
			VshapeY_long_toy.SetXTitle("Y (cm)")
			VshapeY_long_toy.SetYTitle("t (ns)")
			VshapeY_long_toy.Draw("colz")
			gPad.Update()
			stVshapeY_long_toy = VshapeY_long_toy.GetListOfFunctions().FindObject("stats").Clone("stVshapeY")
			stVshapeY_long_toy.SetX1NDC(0.4)
			stVshapeY_long_toy.SetX2NDC(0.7)
			stVshapeY_long_toy.SetY1NDC(0.7)
			stVshapeY_long_toy.SetY2NDC(0.9)
			gPad.Update()
			if advVshapeFit:
				VshapeY_long_fit = TF1("VshapeY_long_fit", vshapeParametricMod, bordL[0], bordL[1], 7)
			else:
				VshapeY_long_fit = TF1("VshapeY_long_fit", vshapeParametric, bordL[0], bordL[1], 3)
			parL = np.array(parL).flatten('C')
			VshapeY_long_fit.SetParameters(parL)
			VshapeY_long_fit.Draw("same")
			C_toyMC.SaveAs(plots_dir + "vshapes_toyMC.png")

	'''raw_input("Press any key to continue\n")

	C_resolL = TCanvas("C_resolL", "C_resolL", 2000, 1000)
	C_resolL.Divide(3,2)
	C_resolL.cd(1)
	VshapeY_long_first.SetXTitle("Y (cm)")
	VshapeY_long_first.SetYTitle("t (ns)")
	VshapeY_long_first.Draw("colz")
	C_resolL.cd(2)
	meanL_first = aSlicesL_first[0].Clone()
	meanL_first.Draw()
	C_resolL.cd(3)
	rmsL_first = aSlicesL_first[2].Clone()
	rmsL_first.Draw()
	C_resolL.cd(4)
	VshapeY_long_second.SetXTitle("Y (cm)")
	VshapeY_long_second.SetYTitle("t (ns)")
	VshapeY_long_second.Draw("colz")
	C_resolL.cd(5)
	meanL_second = aSlicesL_second[0].Clone()
	meanL_second.Draw()
	C_resolL.cd(6)
	rmsL_second = aSlicesL_second[2].Clone()
	rmsL_second.Draw()
	C_resolL.SaveAs(plots_dir + "resolutionL.png")'''


	# canvas for residuals

	c_resU_Det = TCanvas("U_residual_on_detectors", "U_residual_on_detectors", 1280,800)

	c_resUvsU_Det = TCanvas("U_residual_vs_U_on_detectors", "U_residual_vs_U_on_detectors", 1280,800)
	c_resUvsV_Det = TCanvas("U_residual_vs_V_on_detectors", "U_residual_vs_V_on_detectors", 1280,800)


	c_resUvsSlpX_Det = TCanvas("U_residual_vs_SlpX_on_detectors", "U_residual_vs_SlpX_on_detectors", 1280,800)
	c_resUvsSlpY_Det = TCanvas("U_residual_vs_SlpY_on_detectors", "U_residual_vs_SlpY_on_detectors", 1280,800)

	# canvas fill residuals
	c_resU_Det.Divide(4,2)

	c_resUvsU_Det.Divide(4,2)
	c_resUvsV_Det.Divide(4,2)

	c_resUvsSlpX_Det.Divide(4,2)
	c_resUvsSlpY_Det.Divide(4,2)

	for iDet in range(nDets):
		c_resU_Det.cd(iDet+1)
		gStyle.SetStatX(0.9)
		gStyle.SetStatY(0.9)
		resU_Det[iDet].Draw()
		resU_Det[iDet].Fit("gaus","R")
		
		c_resUvsU_Det.cd(iDet+1)
		resUvsU_Det[iDet].Draw("colz")

		c_resUvsV_Det.cd(iDet+1)
		resUvsV_Det[iDet].Draw("colz")

		c_resUvsSlpX_Det.cd(iDet+1)
		resUvsSlpX_Det[iDet].Draw("colz")

		c_resUvsSlpY_Det.cd(iDet+1)
		resUvsSlpY_Det[iDet].Draw("colz")
	
	gStyle.SetStatX(statXdef)
	gStyle.SetStatY(statYdef)

	# canvas draw
	c_resU_Det.Draw()
	c_resU_Det.SaveAs(plots_dir + "U_res.png")

	c_resUvsU_Det.Draw()
	c_resUvsV_Det.Draw()
	c_resUvsU_Det.SaveAs(plots_dir + "UvsU_det_res.png")
	c_resUvsV_Det.SaveAs(plots_dir + "UvsV_det_res.png")

	c_resUvsSlpX_Det.Draw()
	c_resUvsSlpY_Det.Draw()
	c_resUvsSlpX_Det.SaveAs(plots_dir + "UvsSlpX_det_res.png")
	c_resUvsSlpY_Det.SaveAs(plots_dir + "UvsSlpY_det_res.png")

	c_residuals_DUT = TCanvas("XY_sigma_in_DUTs", "XY_sigma_in_DUTs", 1280, 800)
	c_residuals_DUT.Divide(2,2)
	c_residuals_DUT.cd(1)
	residualX_short.SetXTitle("X (cm)")
	residualX_short.Draw()
	# residualX_short.Fit("gaus","R")
	c_residuals_DUT.cd(2)
	residualY_short.SetXTitle("Y (cm)")
	residualY_short.Draw()
	# residualY_short.Fit("gaus","R")
	c_residuals_DUT.cd(3)
	residualX_long.SetXTitle("X (cm)")
	residualX_long.Draw()
	# residualX_long.Fit("gaus","R")
	c_residuals_DUT.cd(4)
	residualY_long.SetXTitle("Y (cm)")
	residualY_long.Draw()
	# residualY_long.Fit("gaus","R")
	c_residuals_DUT.SaveAs(plots_dir + "sigmas_DUT.png")

	c_tdc_hit_cor = TCanvas("Hits_correlation","Hits_correlation", 1280, 800)
	c_tdc_hit_cor.Divide(4,2)
	c_tdc_hit_cor.cd(1)
	TDC_hit_cor.SetXTitle("t (ns)")
	TDC_hit_cor.Draw()
	c_tdc_hit_cor.cd(2)
	TDC_hit_cor_SiPM_up.SetXTitle("t (ns)")
	TDC_hit_cor_SiPM_up.Draw()
	c_tdc_hit_cor.cd(3)
	TDC_hit_cor_SiPM_down.SetXTitle("t (ns)")
	TDC_hit_cor_SiPM_down.Draw()
	c_tdc_hit_cor.cd(4)
	TDC_hit_cor_SiPM_up_down.SetXTitle("t (ns)")
	TDC_hit_cor_SiPM_up_down.Draw()
	c_tdc_hit_cor.cd(5)
	TDC_scin_2_hits.SetXTitle("First hit, t (ns)")
	TDC_scin_2_hits.SetYTitle("Second hit, t (ns)")
	TDC_scin_2_hits.Draw("colz")
	gPad.Update()
	stTDC_scin_2_hits = TDC_scin_2_hits.GetListOfFunctions().FindObject("stats").Clone("stTDC_scin_2_hits")
	stTDC_scin_2_hits.SetX1NDC(0.7)
	stTDC_scin_2_hits.SetX2NDC(0.9)
	stTDC_scin_2_hits.SetY1NDC(0.7)
	stTDC_scin_2_hits.SetY2NDC(0.9)
	gPad.Update()
	c_tdc_hit_cor.cd(6)
	TDC_scin_2_hits_SiPM_up.SetXTitle("T0 hit, t (ns)")
	TDC_scin_2_hits_SiPM_up.SetYTitle("Upstream SiPMs hit, t (ns)")
	TDC_scin_2_hits_SiPM_up.Draw("colz")
	gPad.Update()
	stTDC_scin_2_hits_SiPM_up = TDC_scin_2_hits_SiPM_up.GetListOfFunctions().FindObject("stats").Clone("stTDC_scin_2_hits_SiPM_up")
	stTDC_scin_2_hits_SiPM_up.SetX1NDC(0.1)
	stTDC_scin_2_hits_SiPM_up.SetX2NDC(0.4)
	stTDC_scin_2_hits_SiPM_up.SetY1NDC(0.7)
	stTDC_scin_2_hits_SiPM_up.SetY2NDC(0.9)
	gPad.Update()
	c_tdc_hit_cor.cd(7)
	TDC_scin_2_hits_SiPM_down.SetXTitle("T0 hit, t (ns)")
	TDC_scin_2_hits_SiPM_down.SetYTitle("Downstream SiPMs hit, t (ns)")
	TDC_scin_2_hits_SiPM_down.Draw("colz")
	gPad.Update()
	stTDC_scin_2_hits_SiPM_down = TDC_scin_2_hits_SiPM_down.GetListOfFunctions().FindObject("stats").Clone("stTDC_scin_2_hits_SiPM_down")
	stTDC_scin_2_hits_SiPM_down.SetX1NDC(0.1)
	stTDC_scin_2_hits_SiPM_down.SetX2NDC(0.4)
	stTDC_scin_2_hits_SiPM_down.SetY1NDC(0.7)
	stTDC_scin_2_hits_SiPM_down.SetY2NDC(0.9)
	gPad.Update()
	c_tdc_hit_cor.cd(8)
	TDC_scin_2_hits_SiPM_up_down.SetXTitle("Upstream SiPMs hit, t (ns)")
	TDC_scin_2_hits_SiPM_up_down.SetYTitle("Downstream SiPMs hit, t (ns)")
	TDC_scin_2_hits_SiPM_up_down.Draw("colz")
	gPad.Update()
	stTDC_scin_2_hits_SiPM_up_down = TDC_scin_2_hits_SiPM_up_down.GetListOfFunctions().FindObject("stats").Clone("stTDC_scin_2_hits_SiPM_up_down")
	stTDC_scin_2_hits_SiPM_up_down.SetX1NDC(0.1)
	stTDC_scin_2_hits_SiPM_up_down.SetX2NDC(0.4)
	stTDC_scin_2_hits_SiPM_up_down.SetY1NDC(0.7)
	stTDC_scin_2_hits_SiPM_up_down.SetY2NDC(0.9)
	gPad.Update()
	c_tdc_hit_cor.SaveAs(plots_dir + "hits_cor.png")

#	raw_input("Press any key to continue\n")


if syncFileFlag and not loadSyncFileFlag:
	sync_mamba_tree.AddFriend(sync_vme_tree)
	sync_mamba_tree.Write()
	sync_vme_tree.Write()

	out_file.Close()


VshapeY.SetXTitle("Y (cm)")
VshapeY.SetYTitle("t (ns)")
VshapeY.Write()
if enableToyMC:
	VshapeY_toy.SetXTitle("Y (cm)")
	VshapeY_toy.SetYTitle("t (ns)")
	VshapeY_toy.Write()
TDC_dif_trk.SetXTitle("TDC (ns)")
TDC_dif_trk.Write()
VshapeY_long.SetXTitle("Y (cm)")
VshapeY_long.SetYTitle("t (ns)")
VshapeY_long.Write()
if enableToyMC:
	VshapeY_long_toy.SetXTitle("Y (cm)")
	VshapeY_long_toy.SetYTitle("t (ns)")
	VshapeY_long_toy.Write()
TDC_dif_trk_long.SetXTitle("TDC (ns)")
TDC_dif_trk_long.Write()
if enableToyMC and noiseRedFlag:
	resolutionS_toy.SetXTitle("Y (cm)")
	resolutionS_toy.SetYTitle("residual Y (cm)")
	resolutionS_toy.Write()
if noiseRedFlag:
	resolutionS.SetXTitle("Y (cm)")
	resolutionS.SetYTitle("residual Y (cm)")
	resolutionS.Write()
	resolution_slices[0].SetXTitle("Y (cm)")
	resolution_slices[0].Write()
	resolution_slices[1].SetXTitle("Y (cm)")
	resolution_slices[1].SetYTitle("residual Y (cm)")
	resolution_slices[1].Write()
	resolution_slices[2].SetXTitle("Y (cm)")
	resolution_slices[2].SetYTitle("sigma (cm)")
	resolution_slices[2].Write()
	if enableToyMC:
		resolution_slices_toy[0].SetXTitle("Y (cm)")
		resolution_slices_toy[0].Write()
		resolution_slices_toy[1].SetXTitle("Y (cm)")
		resolution_slices_toy[1].SetYTitle("residual Y (cm)")
		resolution_slices_toy[1].Write()
		resolution_slices_toy[2].SetXTitle("Y (cm)")
		resolution_slices_toy[2].SetYTitle("sigma (cm)")
		resolution_slices_toy[2].Write()
	if enableToyMC:
		resolutionL_toy.SetXTitle("Y (cm)")
		resolutionL_toy.SetYTitle("residual Y (cm)")
		resolutionL_toy.Write()
if noiseRedFlag:
	resolutionL.SetXTitle("Y (cm)")
	resolutionL.SetYTitle("residual Y (cm)")
	resolutionL.Write()
	resolution_slices_long[0].SetXTitle("Y (cm)")
	resolution_slices_long[0].Write()
	resolution_slices_long[1].SetXTitle("Y (cm)")
	resolution_slices_long[1].SetYTitle("residual Y (cm)")
	resolution_slices_long[1].Write()
	resolution_slices_long[2].SetXTitle("Y (cm)")
	resolution_slices_long[2].SetYTitle("sigma (cm)")
	resolution_slices_long[2].Write()
	if enableToyMC:
		resolution_slices_long_toy[0].SetXTitle("Y (cm)")
		resolution_slices_long_toy[0].Write()
		resolution_slices_long_toy[1].SetXTitle("Y (cm)")
		resolution_slices_long_toy[1].SetYTitle("residual Y (cm)")
		resolution_slices_long_toy[1].Write()
		resolution_slices_long_toy[2].SetXTitle("Y (cm)")
		resolution_slices_long_toy[2].SetYTitle("sigma (cm)")
		resolution_slices_long_toy[2].Write()


resol_file.Close()

# end of the script

