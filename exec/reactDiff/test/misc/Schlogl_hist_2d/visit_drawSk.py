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

datafile = "Hist2D.S_k.pair=1.vtk"

# auxiliary paramters

pi = math.pi
kxlo = -pi/dx*float(nx-1)/nx 
kxhi = pi/dx*float(nx+1)/nx 
kylo = -pi/dy*float(ny-1)/ny
kyhi = pi/dy*float(ny+1)/ny

# visit: open datafile and add pseudocolor plot of S(k)

visit.OpenDatabase(datafile)
visit.AddPlot("Pseudocolor","S_k")

visit.HideToolbars()

# colorbar range

pa = visit.PseudocolorAttributes()
pa.minFlag = 1
pa.min = minSk 
pa.maxFlag = 1
pa.max = maxSk 
visit.SetPlotOptions(pa)

# x- and y-ranges of plot and position of plot

va = visit.View2DAttributes()
va.windowCoords = (kxlo,kxhi,kylo,kyhi)
va.viewportCoords = (0.06,0.94,0.11,0.99)
visit.SetView2D(va)

# legend(colorbox) and user/database info on/off | axes

aa = visit.AnnotationAttributes()
aa.legendInfoFlag = 1
aa.userInfoFlag = 0
aa.databaseInfoFlag = 0
aa.axes2D.visible = 1
aa.axes2D.autoSetTicks = 0
aa.axes2D.autoSetScaling = 0
aa.axes2D.lineWidth = 0
aa.axes2D.tickLocation = aa.axes2D.Inside  # Inside, Outside, Both
aa.axes2D.tickAxes = aa.axes2D.All  # Off, Bottom, Left, BottomLeft, All
aa.axes2D.xAxis.title.visible = 0
aa.axes2D.xAxis.label.visible = 0
aa.axes2D.xAxis.tickMarks.visible = 1
aa.axes2D.xAxis.tickMarks.majorMinimum = kxlo
aa.axes2D.xAxis.tickMarks.majorMaximum = kxhi
aa.axes2D.xAxis.tickMarks.minorSpacing = (kxhi-kxlo)/16.
aa.axes2D.xAxis.tickMarks.majorSpacing = (kxhi-kxlo)/4.
aa.axes2D.xAxis.grid = 1
aa.axes2D.yAxis.title.visible = 0
aa.axes2D.yAxis.label.visible = 0
aa.axes2D.yAxis.tickMarks.visible = 1
aa.axes2D.yAxis.tickMarks.majorMinimum = kylo
aa.axes2D.yAxis.tickMarks.majorMaximum = kyhi
aa.axes2D.yAxis.tickMarks.minorSpacing = (kyhi-kylo)/16.
aa.axes2D.yAxis.tickMarks.majorSpacing = (kyhi-kylo)/4.
aa.axes2D.yAxis.grid = 1
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
legend.numberFormat = "%.1f" 
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
