import visit 

# Main function DrawPlot (defined at the bottom) uses the following functions:
# - SetPseudocolorAttributes
# - SetView2dAttributes
# - SetAnnotationAttribues

def SetPseudocolorAttributes(valMin,valMax):
# The minimum and maximum values of color box for the variable values
# are set as constant values valMin and valMax (i.e. not varying in time).

  pa = visit.PseudocolorAttributes()

  # options I like to set

  pa.minFlag = 1
  pa.min = valMin 
  pa.maxFlag = 1
  pa.max = valMax 

  # other options with default values

  pa.scaling = pa.Linear  # Linear, Log, Skew
  pa.skewFactor = 1
  pa.limitsMode = pa.OriginalData  # OriginalData, CurrentPlot
  pa.centering = pa.Natural  # Natural, Nodal, Zonal
  pa.colorTableName = "hot"
  pa.invertColorTable = 0
  pa.opacityType = pa.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
  pa.opacityVariable = ""
  pa.opacity = 1
  pa.opacityVarMin = 0
  pa.opacityVarMax = 1
  pa.opacityVarMinFlag = 0
  pa.opacityVarMaxFlag = 0
  pa.pointSize = 0.05
  pa.pointType = pa.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
  pa.pointSizeVarEnabled = 0
  pa.pointSizeVar = "default"
  pa.pointSizePixels = 2
  pa.lineType = pa.Line  # Line, Tube, Ribbon
  pa.lineStyle = pa.SOLID  # SOLID, DASH, DOT, DOTDASH
  pa.lineWidth = 0
  pa.tubeDisplayDensity = 10
  pa.tubeRadiusSizeType = pa.FractionOfBBox  # Absolute, FractionOfBBox
  pa.tubeRadiusAbsolute = 0.125
  pa.tubeRadiusBBox = 0.005
  pa.varyTubeRadius = 0
  pa.varyTubeRadiusVariable = ""
  pa.varyTubeRadiusFactor = 10
  pa.endPointType = pa.None  # None, Tails, Heads, Both
  pa.endPointStyle = pa.Spheres  # Spheres, Cones
  pa.endPointRadiusSizeType = pa.FractionOfBBox  # Absolute, FractionOfBBox
  pa.endPointRadiusAbsolute = 1
  pa.endPointRadiusBBox = 0.005
  pa.endPointRatio = 2
  pa.renderSurfaces = 1
  pa.renderWireframe = 0
  pa.renderPoints = 0
  pa.smoothingLevel = 0
  pa.legendFlag = 1
  pa.lightingFlag = 1

  visit.SetPlotOptions(pa)


def SetView2dAttributes(domain,viewport):
# set x- and y-ranges of domain from domain argument (e.g. (0,64,0,64))
# set relative position of plot from viewport (e.g. (0.06,0.94,0.11,0.99))

  va = visit.View2DAttributes()

  # options I like to set 

  va.windowCoords = (domain.xlo,domain.xhi,domain.ylo,domain.yhi)
  va.viewportCoords = (viewport.xlo,viewport.xhi,viewport.ylo,viewport.yhi)

  # other options with default values

  va.fullFrameActivationMode = va.Auto  # On, Off, Auto
  va.fullFrameAutoThreshold = 100
  va.xScale = va.LINEAR  # LINEAR, LOG
  va.yScale = va.LINEAR  # LINEAR, LOG
  va.windowValid = 1

  visit.SetView2D(va)


def SetAnnotationAttributes(domain,legendInfo):
# turn on legend (i.e. color bar)
# turn off user and database info
# turn on x- and y-axes with ticks but turn off all annotations

  aa = visit.AnnotationAttributes()

  # options I like to set

  if legendInfo.visible:
    aa.legendInfoFlag = 1
  else:
    aa.legendInfoFlag = 0
  aa.userInfoFlag = 0
  aa.databaseInfoFlag = 0

  # aa.axes2D.visible = 1
  aa.axes2D.visible = 0
  aa.axes2D.autoSetTicks = 0
  aa.axes2D.autoSetScaling = 0
  aa.axes2D.lineWidth = 0
  aa.axes2D.tickLocation = aa.axes2D.Inside  # Inside, Outside, Both
  aa.axes2D.tickAxes = aa.axes2D.All  # Off, Bottom, Left, BottomLeft, All

  aa.axes2D.xAxis.title.visible = 0
  aa.axes2D.xAxis.label.visible = 0
  aa.axes2D.xAxis.tickMarks.visible = 1
  aa.axes2D.xAxis.tickMarks.majorMinimum = domain.xlo 
  aa.axes2D.xAxis.tickMarks.majorMaximum = domain.xhi 
  aa.axes2D.xAxis.tickMarks.minorSpacing = (domain.xhi-domain.xlo)/16.
  aa.axes2D.xAxis.tickMarks.majorSpacing = (domain.xhi-domain.xlo)/4.
  aa.axes2D.xAxis.grid = 1

  aa.axes2D.yAxis.title.visible = 0
  aa.axes2D.yAxis.label.visible = 0
  aa.axes2D.yAxis.tickMarks.visible = 1
  aa.axes2D.yAxis.tickMarks.majorMinimum = domain.ylo
  aa.axes2D.yAxis.tickMarks.majorMaximum = domain.yhi 
  aa.axes2D.yAxis.tickMarks.minorSpacing = (domain.yhi-domain.ylo)/16. 
  aa.axes2D.yAxis.tickMarks.majorSpacing = (domain.yhi-domain.ylo)/4. 
  aa.axes2D.yAxis.grid = 1

  # other options with default values

  aa.backgroundColor = (255, 255, 255, 255)
  aa.foregroundColor = (0, 0, 0, 255)
  aa.gradientBackgroundStyle = aa.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
  aa.gradientColor1 = (0, 0, 255, 255)
  aa.gradientColor2 = (0, 0, 0, 255)
  aa.backgroundMode = aa.Solid  # Solid, Gradient, Image, ImageSphere

  aa.axes2D.xAxis.title.font.font = aa.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
  aa.axes2D.xAxis.title.font.scale = 1
  aa.axes2D.xAxis.title.font.useForegroundColor = 1
  aa.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
  aa.axes2D.xAxis.title.font.bold = 1
  aa.axes2D.xAxis.title.font.italic = 1
  aa.axes2D.xAxis.title.userTitle = 0
  aa.axes2D.xAxis.title.userUnits = 0
  aa.axes2D.xAxis.title.title = "X-Axis"
  aa.axes2D.xAxis.title.units = ""
  aa.axes2D.xAxis.label.font.font = aa.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
  aa.axes2D.xAxis.label.font.scale = 1
  aa.axes2D.xAxis.label.font.useForegroundColor = 1
  aa.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
  aa.axes2D.xAxis.label.font.bold = 1
  aa.axes2D.xAxis.label.font.italic = 1
  aa.axes2D.xAxis.label.scaling = 0

  aa.axes2D.yAxis.title.font.font = aa.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
  aa.axes2D.yAxis.title.font.scale = 1
  aa.axes2D.yAxis.title.font.useForegroundColor = 1
  aa.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
  aa.axes2D.yAxis.title.font.bold = 1
  aa.axes2D.yAxis.title.font.italic = 1
  aa.axes2D.yAxis.title.userTitle = 0
  aa.axes2D.yAxis.title.userUnits = 0
  aa.axes2D.yAxis.title.title = "Y-Axis"
  aa.axes2D.yAxis.title.units = ""
  aa.axes2D.yAxis.label.font.font = aa.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
  aa.axes2D.yAxis.label.font.scale = 1
  aa.axes2D.yAxis.label.font.useForegroundColor = 1
  aa.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
  aa.axes2D.yAxis.label.font.bold = 1
  aa.axes2D.yAxis.label.font.italic = 1
  aa.axes2D.yAxis.label.scaling = 0

  visit.SetAnnotationAttributes(aa)


def SetLegend(legendInfo):
# set color box (i.e. position, size, font, etc.)

  legend = visit.GetAnnotationObject(visit.GetPlotList().GetPlots(0).plotName)

  # options I like to set

  legend.active = 1
  legend.managePosition = 0
  legend.position = (legendInfo.xlo,legendInfo.yhi)
  legend.xScale = legendInfo.xScale 
  legend.yScale = legendInfo.yScale 
  legend.orientation = legend.HorizontalBottom # VerticalRight, VerticalLeft, HorizontalTop, HorizontalBottom
  legend.controlTicks = 1
  legend.numTicks = legendInfo.nTick 
  legend.minMaxInclusive = 1
  legend.drawLabels = legend.Values # None, Values, Labels, Both
  legend.numberFormat = legendInfo.valFormat 
  legend.drawTitle = 0
  legend.drawMinMax = 0
  legend.fontHeight = legendInfo.fontHeight 

  # other options

  legend.suppliedValues = (1,2,3,4,5)
  legend.suppliedLabels = ("A","B","C","D","E")


def SetTimeText(timeTextInfo):
# add text "Time=???" if set==True

  if timeTextInfo.set:

    text = visit.CreateAnnotationObject("Text2D")
    text.visible = 1
    text.active = 1
    text.position = (timeTextInfo.xlo,timeTextInfo.ylo)
    text.height = timeTextInfo.height 
    text.text = "Time=$time"


def DrawPlot(iWindow,plotInfo):
  # choose window
  visit.SetActiveWindow(iWindow)

  # open database
  visit.OpenDatabase(plotInfo.databaseName)

  # create pseudocolor plot
  visit.AddPlot("Pseudocolor",plotInfo.variableName)

  # set pseudocolor attributes
  SetPseudocolorAttributes(plotInfo.valMin,plotInfo.valMax)

  # set view2d attributtes
  SetView2dAttributes(plotInfo.domain,plotInfo.viewport)

  # set annotation attributes
  SetAnnotationAttributes(plotInfo.domain,plotInfo.legendInfo)

  # set legend
  SetLegend(plotInfo.legendInfo)

  # set time text
  SetTimeText(plotInfo.timeTextInfo)

  # draw plot
  visit.DrawPlots()


def SetSaveWindowAttributes(nR,nC,listThumbnail):

  sa = visit.SaveWindowAttributes()

  # options I want to set

  sa.width = listThumbnail[nR*nC-1].xlo+listThumbnail[nR*nC-1].xSize 
  sa.height = listThumbnail[nR*nC-1].ylo+listThumbnail[nR*nC-1].ySize 
  sa.format = sa.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
  sa.family = 0
  sa.saveTiled = 0
  sa.advancedMultiWindowSave = 1

  # other options with default values

  sa.outputToCurrentDirectory = 1
  sa.outputDirectory = "."
  sa.fileName = "visit" 
  sa.screenCapture = 0
  sa.quality = 80
  sa.progressive = 0
  sa.binary = 0
  sa.stereo = 0
  sa.compression = sa.PackBits  # None, PackBits, Jpeg, Deflate
  sa.forceMerge = 0
  sa.resConstraint = sa.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions

  # subwindow attributes

  for i in range(nR):
    for j in range(nC):

      k = i*nC+j+1
      lt = listThumbnail[k-1]
   
      if (k==1):
        win = sa.subWindowAtts.win1
      elif (k==2):
        win = sa.subWindowAtts.win2
      elif (k==3):
        win = sa.subWindowAtts.win3
      elif (k==4):
        win = sa.subWindowAtts.win4
      elif (k==5):
        win = sa.subWindowAtts.win5
      elif (k==6):
        win = sa.subWindowAtts.win6
      elif (k==7):
        win = sa.subWindowAtts.win7
      elif (k==8):
        win = sa.subWindowAtts.win8
      elif (k==9):
        win = sa.subWindowAtts.win9
      elif (k==10):
        win = sa.subWindowAtts.win10
      elif (k==11):
        win = sa.subWindowAtts.win11
      elif (k==12):
        win = sa.subWindowAtts.win12
      elif (k==13):
        win = sa.subWindowAtts.win13
      elif (k==14):
        win = sa.subWindowAtts.win14
      elif (k==15):
        win = sa.subWindowAtts.win15
      elif (k==16):
        win = sa.subWindowAtts.win16

      # options I like to set

      win.position = (lt.xlo,lt.ylo)
      win.size = (lt.xSize,lt.ySize)

      # other options with default values

      win.layer = 0
      win.transparency = 0
      win.omitWindow = 0

  visit.SetSaveWindowAttributes(sa)

  return sa
