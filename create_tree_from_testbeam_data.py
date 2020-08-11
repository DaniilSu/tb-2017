import sys
import os
import numpy as np
from math import *
from ROOT import *

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
	result = 2 * [0] # f1 = t(y) top; f2 = t(y) bottom
	if dt < minimum1 and dt < minimum2:
		result[0] = f1(minimum1)
		result[1] = f2(minimum2)
	elif dt < minimum1 and dt >= minimum2:
		#result[0] = f1(minimum1)
		result[1] = f2(dt)
		return result[1]
	elif dt >= minimum1 and dt < minimum2:
		result[0] = f1(dt)
		#result[1] = f2(minimum2)
		return result[0]
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
		
def interpolatePoints(coord,dtime,errorcoord,errordtime,apex):
	indexDiv = next(i for i,j in enumerate(coord) if j>=apex)
	#dtime = [round(i,4) for i in dtime]
	#dtime = np.around(dtime, decimals=4)
	getcontext().prec = 4
	dtime = [float(Decimal("%.4f" % e)) for e in dtime]
	f_right = interp1d(sorted(dtime[indexDiv:]),sorted(coord[indexDiv:]),kind='slinear') # will be top when inverted
	f_left = interp1d(sorted(dtime[:indexDiv]),sorted(coord[:indexDiv],reverse=True),kind='slinear') # will become bottom
	#dtime_precise_left = insertInBetween(insertInBetween(insertInBetween(sorted(dtime[indexDiv:]))))
	#dtime_precise_right = insertInBetween(insertInBetween(insertInBetween(sorted(dtime[:indexDiv]))))
	#print dtime[indexDiv:]
	#plt.plot(sorted(dtime[indexDiv:]), sorted(coord[indexDiv:]), 'o')
	#plt.plot(dtime_precise_left, f_left(dtime_precise_left), '-')
	#plt.legend(['data', 'linear'], loc='best')
	#plt.xlabel('drift time (ns)')
	#plt.ylabel('Y (cm)')
	#plt.title('Left branch of the V-shape')
	#plt.show()
	#plt.plot(sorted(dtime[:indexDiv]), sorted(coord[:indexDiv],reverse=True), 'o')
	#plt.plot(dtime_precise_right, f_right(dtime_precise_right), '-')
	#plt.legend(['data', 'linear'], loc='best')
	#plt.xlabel('drift time (ns)')
	#plt.ylabel('Y (cm)')
	#plt.title('Right branch of the V-shape')
	#plt.show()
	top_dtime=sorted(dtime[indexDiv:])
	bot_dtime=sorted(dtime[:indexDiv])
	return f_right, f_left, top_dtime[0], bot_dtime[0], top_dtime[-1], bot_dtime[-1]

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
	#a = TCanvas()
	#gre1.SetMarkerColor(kBlue)
	#gre1.SetMarkerStyle(21)
	#gre1.GetXaxis().SetTitle("t_{drift} (ns)")
	#gre1.GetYaxis().SetTitle("Y (cm)")
	#gre1.SetTitle("Left branch of the V-shape")
	#gre1.Draw("AP")
	#a.SaveAs("top_branch_old.pdf")
	gre2.Fit(pf2,"WWQNR")
	pf2.GetParameters(par2)
	gre2.Fit(pf2,"QNR")
	#b = TCanvas()
	#gre2.SetMarkerColor(kRed)
	#gre2.SetMarkerStyle(22)
	#gre2.GetXaxis().SetTitle("t_{drift} (ns)")
	#gre2.GetYaxis().SetTitle("Y (cm)")
	#gre2.SetTitle("Right branch of the V-shape")
	#gre2.Draw("AP")
	#b.SaveAs("bottom_branch_old.pdf")
	#raw_input("Press any key to continue\n")
	pf1.GetParameters(par1)
	pf2.GetParameters(par2)
	return par1, par2, dtime[indexDiv], dtime[indexDiv-1]

def chooseTheBest(nt, tt, tt08, vshapeVal, sigma=None) :
	best = 0
	trigger = -1
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

# Usage
if len(sys.argv) < 4:
	sys.stderr.write('Usage:\t' + str(sys.argv[0]) + ' mamba_root vme_rooot (optional) numberOfRun (optional): bin_size (um) statistics_fraction (0..1.0) input_dir...\n')
	sys.exit(1)

#################################################

# bin_size flag

binSizeFlag = True

# toy MC flag

enableToyMC = True

if enableToyMC:
	enableUniform = True
	sigmaYtrak = 0.006
	sigma1time = 5. # 3.6
	sigma2time = 5. # 35.

# advanced V-shape fitting flag

advVshapeFit = True

# Switch between 6-deg polynomial and interpolation

if advVshapeFit:
	method = "Inter" # Inter or Pol6
	
vshf = "vshape_pars_"
if not binSizeFlag:
	vshf = vshf + "run_" + sys.argv[3] + ".txt"
else:
	vshf = vshf + "run_" + sys.argv[3] + "_bin_" + sys.argv[4] + "_" + sys.argv[5] + ".txt"

if not os.path.exists("outputfiles/"):
	os.makedirs("outputfiles/")

if not binSizeFlag:
	resol_file_name = "coordinates_and_drift_time_" + sys.argv[3] + ".root"
else:
	resol_file_name = "coordinates_and_drift_time_" + sys.argv[3] + "_bin_" + sys.argv[4] + "_" + sys.argv[5] + ".root"
resol_file = TFile("outputfiles/" + resol_file_name, "RECREATE")
if resol_file.IsOpen():
	print("The output file with coordinates and drift times was opened successfully!")

f = open(sys.argv[6]+vshf,"r")
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
		f2S = open(sys.argv[6] + "points_short_" + sys.argv[3] + ".dat", "r")
	else:
		f2S = open(sys.argv[6] + "points_short_" + sys.argv[3] + "_bin_" + sys.argv[4] + "_" + sys.argv[5] + ".dat", "r")
	coordS,dtimeS,errorcoordS,errordtimeS = np.loadtxt(f2S, skiprows=1, unpack=True)
	f2S.close()
	f2L = open(sys.argv[6] + "points_long_" + sys.argv[3] + ".dat", "r")
	coordL,dtimeL,errorcoordL,errordtimeL = np.loadtxt(f2L, skiprows=1, unpack=True)
	f2L.close()
	if method == "Pol6":
		parSinvTop, parSinvBot, SinvTopMin, SinvBotMin = findInverseFunc(coordS,dtimeS,errorcoordS,errordtimeS,apexS)
		parLinvTop, parLinvBot, LinvTopMin, LinvBotMin = findInverseFunc(coordL,dtimeL,errorcoordL,errordtimeL,apexL)
	if method == "Inter":
		fSinterTop, fSinterBot, SinterTopMin, SinterBotMin, SinterTopMax, SinterBotMax = interpolatePoints(coordS,dtimeS,errorcoordS,errordtimeS,apexS)
		fLinterTop, fLinterBot, LinterTopMin, LinterBotMin, LinterTopMax, LinterBotMax = interpolatePoints(coordL,dtimeL,errorcoordL,errordtimeL,apexL)

testbeam_data_file = TFile(sys.argv[1],"read")

mamba_tree = TTree()
vme_tree = TTree()

testbeam_data_file.GetObject("tree",mamba_tree)
testbeam_data_file.GetObject("ADC1",vme_tree)

n_mamba_events = mamba_tree.GetEntries()
print("Number of MAMBA Entries:\t", n_mamba_events, "\n")
n_vme_events = vme_tree.GetEntries()
print("Number of VME Entries:\t", n_vme_events, "\n")

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

# Declaring the new Tree

newTree = TTree("short_tube")
newTree_long = TTree("long_tube")

Y_track_short = Double_t()
Y_track_short_err = Double_t()
Y_short = Double_t()

Y_track_long = Double_t()
Y_track_long_err = Double_t()
Y_long = Double_t()

if enableToyMC:
	Y_toy_long = Double_t()
	Y_toy_short = Double_t()
	Y_track_toy_short = Double_t()
	Y_track_toy_long = Double_t()
	drift_toy_short = Double_t()
	drift_toy_long = Double_t()

drift_short = Double_t()
drift_long = Double_t()

newTree.Branch("Y_track_short",Y_track_short,"Y_track_short/D")
newTree.Branch("Y_track_short_err",Y_track_short_err,"Y_track_short_err/D")
newTree.Branch("Y_short",Y_short,"Y_short/D")

newTree_long.Branch("Y_track_long",Y_track_long,"Y_track_long/D")
newTree_long.Branch("Y_track_long_err",Y_track_long_err,"Y_track_long_err/D")
newTree_long.Branch("Y_long",Y_long,"Y_long/D")

newTree.Branch("drift_time_short",drift_short,"drift_short/D")
newTree_long.Branch("drift_time_long",drift_long,"drift_long/D")

if enableToyMC:
	newTree.Branch("Y_toyMC_short",Y_toy_short,"Y_toy_short/D")
	newTree_long.Branch("Y_toyMC_long",Y_toy_long,"Y_toy_long/D")
	newTree.Branch("Y_track_toyMC_short",Y_track_toy_short,"Y_track_toy_short/D")
	newTree_long.Branch("Y_track_toyMC_long",Y_track_toy_long,"Y_track_toy_long/D")
	newTree.Branch("drift_time_toyMC_short",drift_toy_short,"drift_toy_short/D")
	newTree_long.Branch("drift_time_toyMC_long",drift_toy_long,"drift_toy_long/D")

mamba_tree.GetEntry(0)

packID = mamba_tree.nevent
vmeID = packID - 1
vme_tree.GetEntry(vmeID)

# tracks

k_x = 0
k_y = 0

mambaTmp = 0

print('First event (i = {0} and vmeID = {1}) :'.format(mambaTmp, vmeID))

XCov = np.zeros((2, 2))
YCov = np.zeros((2, 2))

r = TRandom3()

for iEvt in xrange(n_mamba_events-mambaTmp):

	if abs(vme_missed_event) > 100000:
		break

	if (iEvt + mambaTmp + 1) % 100 == 0:
		sys.stdout.write("Processing event:\t" + str(iEvt + mambaTmp +1)+ "\r" )
		sys.stdout.flush()

	mamba_tree.GetEntry(iEvt + mambaTmp)
	
	vmeID = iEvt + mambaTmp
	vme_tree.GetEntry(vmeID)
	
	ntrk = mamba_tree.ntrk
	
	nt00 = vme_tree.nt00
	nt01 = vme_tree.nt01
	
	nt08 = vme_tree.nt08
	
	tt00 = vme_tree.tt00
	tt01 = vme_tree.tt01
	
	tt08 = vme_tree.tt08
	
	if ntrk == 1:
	
		#tracks
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

		if ((nt00==1) and (nt08==1)) :
			Y_track_short = Ytrak[4]
			Y_track_short_err = varY_short
			drift_short = conversion*(tt00[0]-tt08[0])
			if advVshapeFit:
				if enableToyMC:
					if Ytrak[4] < bordS[1] and Ytrak[4] > bordS[0]:
						if enableUniform:
							Ytrak_uni = r.Uniform(bordS[0],bordS[1])
							Ytrak_toy = r.Gaus(Ytrak_uni,sigmaYtrak)
							time_toy = r.Gaus(vshapeParametricMod(Ytrak_uni,parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak_uni,parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
						else:
							Ytrak_toy = r.Gaus(Ytrak[4],sigmaYtrak)
							time_toy = r.Gaus(vshapeParametricMod(Ytrak[4],parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak[4],parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
						Y_track_toy_short = Ytrak_toy
						drift_toy_short = time_toy
						if method == "Inter":
							if time_toy <= SinterTopMax or time_toy <= SinterBotMax:
								Y_toy_short = inverseVshapeParametricInter(time_toy,fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak_toy)
							else:
								choice = r.Rndm()
								if choice > 0.5:
									Y_toy_short = r.Uniform(bordS[0]-0.5,bordS[0])
								else:
									Y_toy_short = r.Uniform(bordS[1],bordS[1]+0.5)
						if method == "Pol6":
							Y_toy_short = inverseVshapeParametricMod(time_toy,parSinvTop,parSinvBot,SinvTopMin,SinvBotMin,Ytrak_toy)
				if method == "Pol6":
					Y_short = inverseVshapeParametricMod(conversion*(tt00[0]-tt08[0]),parSinvTop,parSinvBot,SinvTopMin,SinvBotMin,Ytrak[4])
				if method == "Inter":
					if conversion*(tt00[0]-tt08[0]) <= SinterTopMax or conversion*(tt00[0]-tt08[0]) <= SinterBotMax :
						Y_short = inverseVshapeParametricInter(conversion*(tt00[0]-tt08[0]),fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak[4])
					else:
						choice = r.Rndm()
						if choice > 0.5:
							Y_short = r.Uniform(bordS[0]-0.5,bordS[0])
						else:
							Y_short = r.Uniform(bordS[1],bordS[1]+0.5)
			else:
				Y_short = inverseVshapeParametric(conversion*(tt00[0]-tt08[0]),parS,Ytrak[4])
			newTree.Fill()

		if ((nt00>1) and (nt08==1)) :
			if advVshapeFit:
				bestPoints = chooseTheBest(nt00,tt00,tt08,vshapeParametricMod(Ytrak[4],parS),distSigma)
			else:
				bestPoints = chooseTheBest(nt00,tt00,tt08,vshapeParametric(Ytrak[4],parS),distSigma)
			if len(bestPoints) < 1: continue
			for best in bestPoints:
				Y_track_short = Ytrak[4]
				Y_track_short_err = varY_short
				drift_short = conversion*(tt00[best]-tt08[0])
				if advVshapeFit:
					if enableToyMC:
						if Ytrak[4] < bordS[1] and Ytrak[4] > bordS[0]:
							if enableUniform:
								Ytrak_uni = r.Uniform(bordS[0],bordS[1])
								Ytrak_toy = r.Gaus(Ytrak_uni,sigmaYtrak)
								time_toy = r.Gaus(vshapeParametricMod(Ytrak_uni,parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak_uni,parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
							else:
								Ytrak_toy = r.Gaus(Ytrak[4],sigmaYtrak)
								time_toy = r.Gaus(vshapeParametricMod(Ytrak[4],parS),((sigma2time-sigma1time)/(max(dtimeS)-min(dtimeS)))*vshapeParametricMod(Ytrak[4],parS)+(sigma1time*max(dtimeS)-min(dtimeS)*sigma2time)/(max(dtimeS)-min(dtimeS)))
							Y_track_toy_short = Ytrak_toy
							drift_toy_short = time_toy
							if method == "Inter":
								if time_toy <= SinterTopMax or time_toy <= SinterBotMax:
									Y_toy_short = inverseVshapeParametricInter(time_toy,fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak_toy)
								else:
									choice = r.Rndm()
									if choice > 0.5:
										Y_toy_short = r.Uniform(bordS[0]-0.5,bordS[0])
									else:
										Y_toy_short = r.Uniform(bordS[1],bordS[1]+0.5)
							if method == "Pol6":
								Y_toy_short = inverseVshapeParametricMod(time_toy,parSinvTop,parSinvBot,SinvTopMin,SinvBotMin,Ytrak_toy)
					if method == "Pol6":
						Y_short = inverseVshapeParametricMod(conversion*(tt00[best]-tt08[0]),parSinvTop,parSinvBot,SinvTopMin,SinvBotMin,Ytrak[4])
					if method == "Inter":
						if conversion*(tt00[best]-tt08[0]) <= SinterTopMax or conversion*(tt00[best]-tt08[0]) <= SinterBotMax :
							Y_short = inverseVshapeParametricInter(conversion*(tt00[best]-tt08[0]),fSinterTop,fSinterBot,SinterTopMin,SinterBotMin,SinterTopMax,SinterBotMax,Ytrak[4])
						else:
							choice = r.Rndm()
							if choice > 0.5:
								Y_short = r.Uniform(bordS[0]-0.5,bordS[0])
							else:
								Y_short = r.Uniform(bordS[1],bordS[1]+0.5)
				else:
					Y_short = inverseVshapeParametric(conversion*(tt00[best]-tt08[0]),parS,Ytrak[4])
				newTree.Fill()
				
		if ((nt01==1) and (nt08==1)) :
			Y_track_long = Ytrak[5]
			Y_track_long_err = varY_long
			drift_long = conversion*(tt01[0]-tt08[0])
			if advVshapeFit:
				if enableToyMC:
					if Ytrak[5] < bordL[1] and Ytrak[5] > bordL[0]:
						if enableUniform:
							Ytrak_uni = r.Uniform(bordL[0],bordL[1])
							Ytrak_toy = r.Gaus(Ytrak_uni,sigmaYtrak)
							time_toy = r.Gaus(vshapeParametricMod(Ytrak_uni,parL),((sigma2time-sigma1time)/(max(dtimeL)-min(dtimeL)))*vshapeParametricMod(Ytrak_uni,parL)+(sigma1time*max(dtimeL)-min(dtimeL)*sigma2time)/(max(dtimeL)-min(dtimeL)))
						else:
							Ytrak_toy = r.Gaus(Ytrak[5],sigmaYtrak)
							time_toy = r.Gaus(vshapeParametricMod(Ytrak[5],parL),((sigma2time-sigma1time)/(max(dtimeL)-min(dtimeL)))*vshapeParametricMod(Ytrak[5],parL)+(sigma1time*max(dtimeL)-min(dtimeL)*sigma2time)/(max(dtimeL)-min(dtimeL)))
						Y_track_toy_long = Ytrak_toy
						drift_toy_long = time_toy
						if method == "Inter":
							if time_toy <= LinterTopMax or time_toy <= LinterBotMax:
								Y_toy_long = inverseVshapeParametricInter(time_toy,fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak_toy)
							else:
								choice = r.Rndm()
								if choice > 0.5:
									Y_toy_long = r.Uniform(bordL[0]-0.5,bordL[0])
								else:
									Y_toy_long = r.Uniform(bordL[1],bordL[1]+0.5)
						if method == "Pol6":
							Y_toy_long = inverseVshapeParametricMod(time_toy,parLinvTop,parLinvBot,LinvTopMin,LinvBotMin,Ytrak_toy)
				if method == "Pol6":
					Y_long = inverseVshapeParametricMod(conversion*(tt01[0]-tt08[0]),parLinvTop,parLinvBot,LinvTopMin,LinvBotMin,Ytrak[5])
				if method == "Inter":
					if conversion*(tt01[0]-tt08[0]) <= LinterTopMax or conversion*(tt01[0]-tt08[0]) <= LinterBotMax :
						Y_long = inverseVshapeParametricInter(conversion*(tt01[0]-tt08[0]),fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak[5])
					else:
						choice = r.Rndm()
						if choice > 0.5:
							Y_long = r.Uniform(bordL[0]-0.5,bordL[0])
						else:
							Y_long = r.Uniform(bordL[1],bordL[1]+0.5)
			else:
				Y_long = inverseVshapeParametric(conversion*(tt01[0]-tt08[0]),parL,Ytrak[5])
			newTree_long.Fill()

		if ((nt01>1) and (nt08==1)) :
			if advVshapeFit:
				bestPoints = chooseTheBest(nt01,tt01,tt08,vshapeParametricMod(Ytrak[5],parL),distSigma)
			else:
				bestPoints = chooseTheBest(nt01,tt01,tt08,vshapeParametric(Ytrak[5],parL),distSigma)
			if len(bestPoints) < 1: continue
			for best in bestPoints:
				Y_track_long = Ytrak[5]
				Y_track_long_err = varY_long
				drift_long = conversion*(tt01[best]-tt08[0])
				if advVshapeFit:
					if enableToyMC:
						if Ytrak[5] < bordL[1] and Ytrak[5] > bordL[0]:
							if enableUniform:
								Ytrak_uni = r.Uniform(bordL[0],bordL[1])
								Ytrak_toy = r.Gaus(Ytrak_uni,sigmaYtrak)
								time_toy = r.Gaus(vshapeParametricMod(Ytrak_uni,parL),((sigma2time-sigma1time)/(max(dtimeL)-min(dtimeL)))*vshapeParametricMod(Ytrak_uni,parL)+(sigma1time*max(dtimeL)-min(dtimeL)*sigma2time)/(max(dtimeL)-min(dtimeL)))
							else:
								Ytrak_toy = r.Gaus(Ytrak[5],sigmaYtrak)
								time_toy = r.Gaus(vshapeParametricMod(Ytrak[5],parL),((sigma2time-sigma1time)/(max(dtimeL)-min(dtimeL)))*vshapeParametricMod(Ytrak[5],parL)+(sigma1time*max(dtimeL)-min(dtimeL)*sigma2time)/(max(dtimeL)-min(dtimeL)))
							Y_track_toy_long = Ytrak_toy
							drift_toy_long = time_toy
							if method == "Inter":
								if time_toy <= LinterTopMax or time_toy <= LinterBotMax:
									Y_toy_long = inverseVshapeParametricInter(time_toy,fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak_toy)
								else:
									choice = r.Rndm()
									if choice > 0.5:
										Y_toy_long = r.Uniform(bordL[0]-0.5,bordL[0])
									else:
										Y_toy_long = r.Uniform(bordL[1],bordL[1]+0.5)
							if method == "Pol6":
								Y_toy_long = inverseVshapeParametricMod(time_toy,parLinvTop,parLinvBot,LinvTopMin,LinvBotMin,Ytrak_toy)
					if method == "Pol6":
						Y_long = inverseVshapeParametricMod(conversion*(tt01[best]-tt08[0]),parLinvTop,parLinvBot,LinvTopMin,LinvBotMin,Ytrak[5])
					if method == "Inter":
						if conversion*(tt01[best]-tt08[0]) <= LinterTopMax or conversion*(tt01[best]-tt08[0]) <= LinterBotMax :
							Y_long = inverseVshapeParametricInter(conversion*(tt01[best]-tt08[0]),fLinterTop,fLinterBot,LinterTopMin,LinterBotMin,LinterTopMax,LinterBotMax,Ytrak[5])
						else:
							choice = r.Rndm()
							if choice > 0.5:
								Y_long = r.Uniform(bordL[0]-0.5,bordL[0])
							else:
								Y_long = r.Uniform(bordL[1],bordL[1]+0.5)
				else:
					Y_long = inverseVshapeParametric(conversion*(tt01[best]-tt08[0]),parL,Ytrak[5])
				newTree_long.Fill()

newTree.Write()
newTree_long.Write()

resol_file.Close()
