fname = "dimerization.S_k.pair=2.vtk"

outname0 = "Sk_vtk_out0"
outname1 = "Sk_vtk_out1"
outname2 = "Sk_vtk_out2"

out0 = open(outname0,'w')
out1 = open(outname1,'w')
out2 = open(outname2,'w')

status = 0
with open(fname) as f:
  for line in f:
    if (status==0):
      out0.write(line) 
      if line.startswith("LOOKUP_TABLE"):
        status = 1
    elif status==1:
      if line.startswith("POINT_DATA"):
        status = 2
        out2.write(line) 
      else:
        numbers = []
        for word in line.split():
          try:
            numbers.append(float(word))
          except ValueError:
            pass
        for i in range(len(numbers)):
          out1.write("%.12e\n" % numbers[i])
    else:
      out2.write(line)

out0.close()
out1.close()
out2.close()

print "** from %s: %s, %s, and %s generated" % (fname,outname0,outname1,outname2)
