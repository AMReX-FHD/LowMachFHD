# you can execute this script as follows: visit -nowin -cli -s (name of this script)

import sys
import visit

# input parameters

databaseName = "movie.visit"
variableName = "n1"
valMin = 0.5
valMax = 2.5
colorTableName = "gray"
xlo = 0			# domain range
xhi = 256
ylo = 0
yhi = 256
width = 256		# image size
height = 256
filename = "test"	# image filename prefix

# visit script

visit.OpenDatabase(databaseName)
visit.AddPlot("Pseudocolor",variableName)

pa = visit.PseudocolorAttributes()
pa.minFlag = 1
pa.min = valMin
pa.maxFlag = 1
pa.max = valMax
pa.colorTableName = colorTableName 
visit.SetPlotOptions(pa)

aa = visit.AnnotationAttributes()
aa.axes2D.visible = 0
aa.userInfoFlag = 0
aa.databaseInfoFlag = 0
aa.legendInfoFlag = 0
visit.SetAnnotationAttributes(aa)

va = visit.View2DAttributes()
va.windowCoords = (xlo,xhi,ylo,yhi)
va.viewportCoords = (0,1,0,1)
visit.SetView2D(va)

visit.DrawPlots()

sa = visit.SaveWindowAttributes()
sa.width = width
sa.height = height 
sa.format = sa.PNG  
sa.family = 0

nTimestep = TimeSliderGetNStates()
for iTimestep in range(nTimestep):
  TimeSliderSetState(iTimestep)
  tag = "%04d.png" % iTimestep
  sa.fileName = filename+tag
  visit.SetSaveWindowAttributes(sa)
  visit.SaveWindow()

sys.exit(0)
