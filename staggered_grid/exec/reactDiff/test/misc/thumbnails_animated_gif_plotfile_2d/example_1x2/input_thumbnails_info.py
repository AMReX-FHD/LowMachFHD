##########################
# input: thumbnails info #
###############i##########

# numbers of rows and columns  (R=row,C=column)
nR = 1		
nC = 2

# size of plots and margins (H=horizontal,V=vertical)
# H=horizontal, 2*nC+1 components, from left to right
# (margin,plot,margin,plot, ... ,margin,plot,margin)
sizeH = (50,256,50,256,50)	
# V=vertical, 2*nR+1 components, from bottom to top
# (margin,plot,margin,plot, ... ,margin,plot,margin)
sizeV = (100,256,50)     

# plot number examples:
# 	(1x2) 	(2x1)	(2x2)	(2x3) 	(3x2)
# 	1 2   	2	3 4	4 5 6	5 6
#       	1       1 2	1 2 3	3 4
#					1 2

# databaseName, variableName, valMin, valMax, nTick, visible (color box)
listData = []
listData += ["../../Schlogl_spherical/RAD32.visit","n1",0.5,1.6,3,False]	# plot1
listData += ["../../Schlogl_spherical/RAD64.visit","n1",0.5,1.6,3,True]		# plot2

# domain range
domain_xlo = 0
domain_xhi = 256 
domain_ylo = 0
domain_yhi = 256

# including text "Time=$time" for each plot
listIncludeTimeText = [False,True]

# some parameters for color box legend
param_cb1 = -0.07	# (-) left (+) right w.r.t. left-side of plot
param_cb2 = 0.1 	# larger value for further below bottom of plot
param_cb3 = 2.2		# factor for the width of the legend 
param_cb4 = 1.2		# factor for the height of the legend		
param_cb5 = 0.08	# factor for the font height
cb_val_format = "%.2f"	# format for color box values

# some parameters for time text 
param_tt1 = 0.62	# increase value for moving text further to the right
param_tt2 = 0.2		# increase value for moving text further below
param_tt3 = 0.07	# factor for the font height

