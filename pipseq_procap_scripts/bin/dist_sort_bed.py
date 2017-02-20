import numpy as np

def dist_sort(file1,file1_ncol,file1_idcol,file2,file2_ncol,file2_idcol,output1_name,output2_name):

	data00=[]

	data1=open(file1,'r')
	data10=[]
	data11={}
	for records in data1:
		tmp=[x.strip() for x in records.split('\t')]
		data10.append(tmp[file1_idcol])
		data11[tmp[file1_idcol]]=tmp



	data2=open(file2,'r')
	data20=[]
	data21={}
	for records in data2:
		tmp=[x.strip() for x in records.split('\t')]
		data20.append(tmp[file2_idcol])
		data21[tmp[file2_idcol]]=tmp
		if tmp[file2_idcol] in data11:
			data00.append(tmp[file2_idcol])


	data30=[]
	data31=[]
	for records in data00:
		data30.append([ data11[records],data21[records]])
		if data11[records][5] == "+":
			data31.append( int(data11[records][file1_ncol])-int(data21[records][file2_ncol]) )
		else:
			data31.append( -int(data11[records][file1_ncol])+int(data21[records][file2_ncol]) )


	data31=np.array(data31)
	print(data31)
	data30=np.array(data30)
	data30 = data30[data31.astype(np.int).argsort()] 

	data1_sort=open(output1_name,'w')
	data2_sort=open(output2_name,'w')

	for records in data30:
		for i in range(0,len(records[0])-1):
			data1_sort.write(str(records[0][i])+'\t')
		data1_sort.write(str(records[0][len(records[0])-1])+'\n')
		for i in range(0,len(records[1])-1):
			data2_sort.write(str(records[1][i])+'\t')
		data2_sort.write(str(records[1][len(records[1])-1])+'\n')

	data1.close()
	data2.close()
	data1_sort.close()
	data2_sort.close()



import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ht:m:o:s:n:p:a:b:")
	except getopt.GetoptError:
		print 'python vlookup.py -t first_table_file -m first_table_ID -s second_table_file -n second_table_ID -o output_name'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python vlookup.py -t first_table_file -m first_table_ID -s second_table_file -n second_table_ID -o output_name'
			sys.exit()
		elif opt=="-t":
			file1=str(arg.strip())
		elif opt=="-m":
			file1_ncol=int(arg.strip())
		elif opt=="-o":
			file1_idcol=int(arg.strip())
		elif opt=="-s":
			file2=str(arg.strip())
		elif opt=="-n":
			file2_ncol=int(arg.strip())
		elif opt=="-p":
			file2_idcol=int(arg.strip())
		elif opt=="-a":
			output1_name=str(arg.strip())
		elif opt=="-b":
			output2_name=str(arg.strip())
	dist_sort(file1,file1_ncol,file1_idcol,file2,file2_ncol,file2_idcol,output1_name,output2_name)

if __name__=="__main__":
	main(sys.argv[1:])
