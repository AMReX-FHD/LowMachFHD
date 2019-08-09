SRCVTK=$1
DESVTK=$2
EXEC=$3/coarsen_vtk_128_2
PART0=$3/64x64.vtk_part0

# check line numbers (start and end lines of each data)

cp $PART0 $DESVTK 
sed -n '43,1863p' $SRCVTK | $EXEC >> $DESVTK		# density
echo "FIELD FieldData 4" >> $DESVTK 
echo "Temperature 1 4096 double" >> $DESVTK 
sed -n '1866,3686p' $SRCVTK | $EXEC >> $DESVTK		# temperature
echo "Scalar1 1 4096 double" >> $DESVTK
sed -n '3688,5508p' $SRCVTK | $EXEC >> $DESVTK		# scalar1 = U
echo "Scalar2 1 4096 double" >> $DESVTK
sed -n '5510,7330p' $SRCVTK | $EXEC >> $DESVTK		# scalar2 = V
echo "Scalar3 1 4096 double" >> $DESVTK	
sed -n '7332,9152p' $SRCVTK | $EXEC >> $DESVTK		# scalar3 = W
echo "POINT_DATA 4225" >> $DESVTK 
echo "" >> $DESVTK
