import numpy as np

def subtract_bg(signal_file,bg_file,output_file):
	data_out=open(output_file,'w')

	print('read signal data')
	data0=open(signal_file,'r')
	data_out.write(data0.readline())
	data00_name=[]
	data00_reads=[]
	for records in data0:
		tmp=records.split()
		data00_name.append(tmp[0])
		data00_reads.append(tmp[2:])
	data00_reads=np.array(data00_reads,dtype=float)

	print('read background data')
	data1=open(bg_file,'r')
	data1.readline()
	data10_reads=[]
	for records in data1:
		tmp=records.split()
		data10_reads.append(tmp[2:])
	data10_reads=np.array(data10_reads,dtype=float)

	print('subtract the background file')
	dif_reads=data00_reads-data10_reads
	dif_reads = dif_reads.clip(min=0)
	print('write bg_subtracted data')
	for name,reads in zip(data00_name,dif_reads):
		data_out.write(name+'\t'+name+'\t')
		for i in range(0,len(reads)-1):
			data_out.write(str(reads[i])+'\t')
		data_out.write(str(reads[-1])+'\n')

	data_out.close()
	data0.close()
	data1.close()

############################################################################
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hs:b:o:")
	except getopt.GetoptError:
		print 'python subtract_bg.py -s signal_file -b background_file -o output_file'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python subtract_bg.py -s signal_file -b background_file -o output_file'
			sys.exit()
		elif opt=="-s":
			signal_file=str(arg.strip())
		elif opt=="-b":
			bg_file=str(arg.strip())
		elif opt=="-o":
			output_file=str(arg.strip())

	subtract_bg(signal_file,bg_file,output_file)

if __name__=="__main__":
	main(sys.argv[1:])


