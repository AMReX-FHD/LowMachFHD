import visit
import math
import sys

# input parameters

nx = 64
ny = 64
dx = 1.
dy = 1.
minSk = 0.9
maxSk = 1.1

datafile = "Hist3D.S_k.pair=1.vtk"

# auxiliary paramters

pi = math.pi
kxlo = -pi/dx*float(nx-1)/nx 
kxhi = pi/dx*float(nx+1)/nx 
kylo = -pi/dy*float(ny-1)/ny
kyhi = pi/dy*float(ny+1)/ny

# visit: open datafile and add pseudocolor plot of S(k) and add operator threeslice

visit.HideToolbars()

visit.OpenDatabase(datafile)
visit.AddPlot("Pseudocolor","S_k")

AddOperator("ThreeSlice")
tatts = ThreeSliceAttributes()
tatts.x = 0 
tatts.y = 0
tatts.z = 0
SetOperatorOptions(tatts)

# colorbar range

pa = visit.PseudocolorAttributes()
pa.minFlag = 1
pa.min = minSk 
pa.maxFlag = 1
pa.max = maxSk 
visit.SetPlotOptions(pa)

# x- and y-ranges of plot and position of plot

va = visit.View3DAttributes()
va.viewNormal = (0, 0, 1)
va.focus = (0, 0, 0)
va.viewUp = (0, 1, 0)
visit.SetView3D(va)

# legend(colorbox) and user/database info on/off | axes

aa = visit.AnnotationAttributes()
aa.legendInfoFlag = 1
aa.userInfoFlag = 0
aa.databaseInfoFlag = 0
aa.axes3D.visible = 0
aa.axes3D.autoSetTicks = 0
aa.axes3D.autoSetScaling = 0
aa.axes3D.lineWidth = 0
aa.axes3D.tickLocation = aa.axes2D.Inside  # Inside, Outside, Both
#aa.axes3D.tickAxes = aa.axes3D.All  # Off, Bottom, Left, BottomLeft, All
aa.axes3D.xAxis.title.visible = 0
aa.axes3D.xAxis.label.visible = 0
aa.axes3D.xAxis.tickMarks.visible = 1
aa.axes3D.xAxis.tickMarks.majorMinimum = kxlo
aa.axes3D.xAxis.tickMarks.majorMaximum = kxhi
aa.axes3D.xAxis.tickMarks.minorSpacing = (kxhi-kxlo)/16.
aa.axes3D.xAxis.tickMarks.majorSpacing = (kxhi-kxlo)/4.
aa.axes3D.xAxis.grid = 1
aa.axes3D.yAxis.title.visible = 0
aa.axes3D.yAxis.label.visible = 0
aa.axes3D.yAxis.tickMarks.visible = 1
aa.axes3D.yAxis.tickMarks.majorMinimum = kylo
aa.axes3D.yAxis.tickMarks.majorMaximum = kyhi
aa.axes3D.yAxis.tickMarks.minorSpacing = (kyhi-kylo)/16.
aa.axes3D.yAxis.tickMarks.majorSpacing = (kyhi-kylo)/4.
aa.axes3D.yAxis.grid = 1
visit.SetAnnotationAttributes(aa)

# colorbox setting

legend = visit.GetAnnotationObject(visit.GetPlotList().GetPlots(0).plotName)
legend.active = 1
legend.managePosition = 0
legend.position = (0.5,0.08)
legend.xScale = 1.2 
legend.yScale = 0.4
legend.orientation = legend.HorizontalBottom # VerticalRight, VerticalLeft, HorizontalTop, HorizontalBottom
legend.controlTicks = 1
legend.numTicks = 5
legend.minMaxInclusive = 1
legend.drawLabels = legend.Values # None, Values, Labels, Both
legend.numberFormat = "%.2f" 
legend.drawTitle = 0
legend.drawMinMax = 0
legend.fontHeight = 0.08

# draw plot

visit.DrawPlots()

# save figure

sa = visit.SaveWindowAttributes()
sa.width = 500 
sa.height = 500 
sa.format = sa.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
sa.family = 0
sa.saveTiled = 0
sa.fileName = "Sk"
sa.advancedMultiWindowSave = 0
visit.SetSaveWindowAttributes(sa)

visit.SaveWindow()

# terminite

sys.exit()
