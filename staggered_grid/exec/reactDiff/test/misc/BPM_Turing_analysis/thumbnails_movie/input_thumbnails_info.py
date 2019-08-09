##########################
# input: thumbnails info #
###############i##########

# numbers of rows and columns  (R=row,C=column)
nR = 5
nC = 3

# size of plots and margins (H=horizontal,V=vertical)
# H=horizontal, 2*nC+1 components, from left to right
# (margin,plot,margin,plot, ... ,margin,plot,margin)
sizeH = (20,128,20,128,20,128,20)
# V=vertical, 2*nR+1 components, from bottom to top
# (margin,plot,margin,plot, ... ,margin,plot,margin)
sizeV = (50,128,20,128,20,128,20,128,20,128,40)

# plot number examples:
# 	(1x2) 	(2x1)	(2x2)	(2x3) 	(3x2)
# 	1 2   	2	3 4	4 5 6	5 6
#       	1       1 2	1 2 3	3 4
#					1 2

# databaseName, variableName, valMin, valMax, nTick, visible (color box)
listData = []

listData += ["../../256x256_IMMID_dt0.5_factor4/movie.visit","n1",600,2000,3,False]		# plot1
listData += ["../../128x128_IMMID_dt1_factor2/movie.visit","n1",600,2000,3,False]		# plot2
listData += ["../../64x64_ImMid_dt1/movie.visit","n1",600,2000,3,True]				# plot3
listData += ["../../256x256_IMMID_dt0.25_factor4/movie.visit","n1",600,2000,3,False]		# plot4
listData += ["../../128x128_IMMID_dt0.5_factor2/movie.visit","n1",600,2000,3,False]		# plot5
listData += ["../../64x64_ImMid_dt0.5/movie.visit","n1",600,2000,3,False]			# plot6
listData += ["../../256x256_EXMID_dt0.025_factor4/movie.visit","n1",600,2000,3,False]		# plot7
listData += ["../../128x128_EXMID_dt0.1_factor2/movie.visit","n1",600,2000,3,False]		# plot8
listData += ["../../64x64_ExMid_dt0.25/movie.visit","n1",600,2000,3,False]			# plot9
listData += ["../../256x256_RDME_dt0.05_factor4/movie.visit","n1",600,2000,3,False]		# plot10
listData += ["../../128x128_RDME_dt0.2_factor2/movie.visit","n1",600,2000,3,False]		# plot11
listData += ["../../64x64_RDME_dt0.5/movie.visit","n1",600,2000,3,False]			# plot12
listData += ["../../256x256_Particle_dt0.05_factor4/movie.visit","Scalar1",600,2000,3,False]	# plot13
listData += ["../../128x128_Particle_dt0.2_factor2/movie.visit","Scalar1",600,2000,3,False]	# plot14
listData += ["../../64x64_Particle_dt0.5/movie.visit","Scalar1",600,2000,3,False]		# plot15

#listData += ["../../256x256_IMMID_dt0.5_factor4/movie.visit","n3",0,300,3,False]		# plot1
#listData += ["../../128x128_IMMID_dt1_factor2/movie.visit","n3",0,300,3,False]			# plot2
#listData += ["../../64x64_ImMid_dt1/movie.visit","n3",0,300,3,True]				# plot3
#listData += ["../../256x256_IMMID_dt0.25_factor4/movie.visit","n3",0,300,3,False]		# plot4
#listData += ["../../128x128_IMMID_dt0.5_factor2/movie.visit","n3",0,300,3,False]		# plot5
#listData += ["../../64x64_ImMid_dt0.5/movie.visit","n3",0,300,3,False]				# plot6
#listData += ["../../256x256_EXMID_dt0.025_factor4/movie.visit","n3",0,300,3,False]		# plot7
#listData += ["../../128x128_EXMID_dt0.1_factor2/movie.visit","n3",0,300,3,False]		# plot8
#listData += ["../../64x64_ExMid_dt0.25/movie.visit","n3",0,300,3,False]			# plot9
#listData += ["../../256x256_RDME_dt0.05_factor4/movie.visit","n3",0,300,3,False]		# plot10
#listData += ["../../128x128_RDME_dt0.2_factor2/movie.visit","n3",0,300,3,False]		# plot11
#listData += ["../../64x64_RDME_dt0.5/movie.visit","n3",0,300,3,False]				# plot12
#listData += ["../../256x256_Particle_dt0.05_factor4/movie.visit","Scalar3",0,300,3,False]	# plot13
#listData += ["../../128x128_Particle_dt0.2_factor2/movie.visit","Scalar3",0,300,3,False]	# plot14
#listData += ["../../64x64_Particle_dt0.5/movie.visit","Scalar3",0,300,3,False]			# plot15

#listData += ["../../256x256_IMMID_dt0.5_factor4/movie.visit","n3",-10,0,3,False]		# plot1
#listData += ["../../128x128_IMMID_dt1_factor2/movie.visit","n3",-10,0,3,False]			# plot2
#listData += ["../../64x64_ImMid_dt1/movie.visit","n3",-10,0,3,True]				# plot3
#listData += ["../../256x256_IMMID_dt0.25_factor4/movie.visit","n3",-10,0,3,False]		# plot4
#listData += ["../../128x128_IMMID_dt0.5_factor2/movie.visit","n3",-10,0,3,False]		# plot5
#listData += ["../../64x64_ImMid_dt0.5/movie.visit","n3",-10,0,3,False]				# plot6
#listData += ["../../256x256_EXMID_dt0.025_factor4/movie.visit","n3",-10,0,3,False]		# plot7
#listData += ["../../128x128_EXMID_dt0.1_factor2/movie.visit","n3",-10,0,3,False]		# plot8
#listData += ["../../64x64_ExMid_dt0.25/movie.visit","n3",-10,0,3,False]			# plot9
#listData += ["../../256x256_RDME_dt0.05_factor4/movie.visit","n3",-10,0,3,False]		# plot10
#listData += ["../../128x128_RDME_dt0.2_factor2/movie.visit","n3",-10,0,3,False]		# plot11
#listData += ["../../64x64_RDME_dt0.5/movie.visit","n3",-10,0,3,False]				# plot12
#listData += ["../../256x256_Particle_dt0.05_factor4/movie.visit","Scalar3",-10,0,3,False]	# plot13
#listData += ["../../128x128_Particle_dt0.2_factor2/movie.visit","Scalar3",-10,0,3,False]	# plot14
#listData += ["../../64x64_Particle_dt0.5/movie.visit","Scalar3",-10,0,3,False]			# plot15

# domain range
domain_xlo = 0
domain_xhi = 32
domain_ylo = 0
domain_yhi = 32

# including text "Time=$time" for each plot
listIncludeTimeText = [False,False,True,False,False,False,False,False,False,False,False,False,False,False,False]

# some parameters for color box legend
param_cb1 = -0.07	# (-) left (+) right w.r.t. left-side of plot
param_cb2 = 0.1 	# larger value for further below bottom of plot
param_cb3 = 2.2		# factor for the width of the legend 
param_cb4 = 1.2		# factor for the height of the legend		
param_cb5 = 0.08	# factor for the font height
cb_val_format = "%.0f"	# format for color box values

# some parameters for time text 
param_tt1 = 0.62	# increase value for moving text further to the right
param_tt2 = 0.2		# increase value for moving text further below
param_tt3 = 0.07	# factor for the font height

