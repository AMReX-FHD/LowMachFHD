SRCVTK=$1
DESVTK=$2
EXEC=$3/coarsen_vtk_256_4
PART0=$3/64x64.vtk_part0

# check line numbers (start and end lines of each data)

cp $PART0 $DESVTK 
sed -n '71,7352p' $SRCVTK | $EXEC >> $DESVTK		# density
echo "FIELD FieldData 4" >> $DESVTK 
echo "Temperature 1 4096 double" >> $DESVTK 
sed -n '7355,14636p' $SRCVTK | $EXEC >> $DESVTK		# temperature
echo "Scalar1 1 4096 double" >> $DESVTK
sed -n '14638,21919p' $SRCVTK | $EXEC >> $DESVTK	# scalar1 = U
echo "Scalar2 1 4096 double" >> $DESVTK
sed -n '21921,29202p' $SRCVTK | $EXEC >> $DESVTK	# scalar2 = V
echo "Scalar3 1 4096 double" >> $DESVTK	
sed -n '29204,36485p' $SRCVTK | $EXEC >> $DESVTK	# scalar3 = W
echo "POINT_DATA 4225" >> $DESVTK 
echo "" >> $DESVTK
