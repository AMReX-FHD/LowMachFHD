import sys

#################
# get arguments #
#################

# sys.argv[0] contains the name of this python script
# for sys.argv[1], provide thumbnails info file (required)
# for sys.argv[2], provide PNG filename (required)
# for sys.argv[3], which is optional, choose 
# either ONLY_FIRST_FRAME (to obtain the first frame for checking arrangement)
# or ALL_FRAMES (to obtain all frames for generating movie)

if (len(sys.argv)==3):
  filenameThumbnailsInfo = sys.argv[1]
  filenamePNG = sys.argv[2]
  onlyFirstFrame = False
elif (len(sys.argv)==4):
  filenameThumbnailsInfo = sys.argv[1]
  filenamePNG = sys.argv[2]
  if (sys.argv[3]=="ONLY_FIRST_FRAME"):
    onlyFirstFrame = True
  elif (sys.argv[3]=="ALL_FRAMES"):
    onlyFirstFrame = False
  else:
    print "Error: ", sys.argv
    print "filenamePNG ONLY_FIRST_FRAME|ALL_FRAMES"
    sys.exit()
else:
  print "Error: ", sys.argv
  print "filenamePNG ONLY_FIRST_FRAME|ALL_FRAMES"
  sys.exit()

if onlyFirstFrame:
  print "** filenamePNG="+filenamePNG+" (ONLY_FIRST_FRAME)"
else:
  print "** filenamePNG="+filenamePNG+" (ALL_FRAMES)"

#############################
# read thumbnails info file #
#############################

execfile(filenameThumbnailsInfo)

######################
# check input values #
######################

if len(sizeH)!=2*nC+1:
  print "Error: Check nC and sizeH"
  sys.exit()

if len(sizeV)!=2*nR+1:
  print "Error: Check nR and sizeV"
  sys.exit()

if len(listData)!=6*nC*nR:
  print "Error: Check listData"
  sys.exit()

if len(listIncludeTimeText)!=nC*nR:
  print "Error: Check listIncludeTimeText"
  sys.exit()

#############################
# definitions of data types #
#############################

from collections import namedtuple

Thumbnail = namedtuple("Thumbnail","xlo ylo xSize ySize plotInfo")

PlotInfo = namedtuple("PlotInfo",
  "databaseName domain variableName valMin valMax viewport \
  legendInfo timeTextInfo")

LegendInfo = namedtuple("LegendInfo",
  "visible xlo yhi xScale yScale nTick valFormat fontHeight")

TimeTextInfo = namedtuple("TimeTextInfo","set xlo ylo height")

LoHi2d = namedtuple("LoHi2d","xlo xhi ylo yhi") 

###########################
# construct listThumbnail #
###########################

listThumbnail = []

ylo = 0

for i in range(nR):

  bmargin = sizeV[2*i]-sizeV[2*i]/10 if i>0 else sizeV[0]
  tmargin = sizeV[2*i+2]/10 if i<nR-1 else sizeV[2*nR]
  ySize = bmargin+sizeV[2*i+1]+tmargin 

  ratio_yB = bmargin/float(ySize)
  ratio_yT = tmargin/float(ySize)
  ratio_yP = 1-ratio_yB-ratio_yT

  xlo = 0

  for j in range(nC):

    k = i*nC+j

    lmargin = sizeH[2*j]-sizeH[2*j]/10 if j>0 else sizeH[0]
    rmargin = sizeH[2*j+2]/10 if j<nC-1 else sizeH[2*nC]
    xSize = lmargin+sizeH[2*j+1]+rmargin 

    ratio_xL = lmargin/float(xSize)
    ratio_xR = rmargin/float(xSize)
    ratio_xP = 1-ratio_xL-ratio_xR

    databaseName = listData[6*k]
    variableName = listData[6*k+1]
    valMin = listData[6*k+2]
    valMax = listData[6*k+3]
    nTick = listData[6*k+4]
    visible = listData[6*k+5]

    domain = LoHi2d(domain_xlo,domain_xhi,domain_ylo,domain_yhi)
    viewport = LoHi2d(ratio_xL,1-ratio_xR,ratio_yB,1-ratio_yT)
    legendInfo = LegendInfo(visible,
      ratio_xL+param_cb1*ratio_xP,ratio_yB-param_cb2*ratio_yP,
      param_cb3*ratio_xP,param_cb4*ratio_yP,nTick,
      cb_val_format,param_cb5*ratio_yP)
    timeTextInfo = TimeTextInfo(listIncludeTimeText[k],
      ratio_xL+param_tt1*ratio_xP,ratio_yB-param_tt2*ratio_yP,
      param_tt3*ratio_yP) 

    plotInfo = PlotInfo(databaseName,domain,variableName,valMin,valMax,
      viewport,legendInfo,timeTextInfo)

    listThumbnail.append(Thumbnail(xlo,ylo,xSize,ySize,plotInfo))

    xlo += xSize 

  ylo += ySize

##############
# draw plots #
##############

from movie_modules import DrawPlot

for i in range(nR):
  for j in range(nC):

    k = i*nC+j

    if k>0:
      AddWindow()

    DrawPlot(k+1,listThumbnail[k].plotInfo)

################################
# save plots at each time step #
################################

from movie_modules import SetSaveWindowAttributes

sa = SetSaveWindowAttributes(nR,nC,listThumbnail)

nTimestep = TimeSliderGetNStates() if not onlyFirstFrame else 1
nPlot = nR*nC

for iTimestep in range(nTimestep):

  for iPlot in range(nPlot):
    SetActiveWindow(iPlot+1)
    TimeSliderSetState(iTimestep)

  tag = "%04d.png" % iTimestep
  sa.fileName = filenamePNG+tag
  visit.SetSaveWindowAttributes(sa)
  visit.SaveWindow()

#############
# terminate #
#############

sys.exit()
