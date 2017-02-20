import numpy as np
def relocate_TSS(bed_file,readscount_file,output_file,searchTSS_win,output_win):
	data0=open(bed_file,'r')
	data00=[]

	for records in data0:
		data00.append(records.split())

	data01=data00
	data02=data00

	data1=open(readscount_file,'r')
	data1.readline()
	names=[]
	data10=[]

	for pk,records in zip(data01,data1):
		if pk[5]=='+':
			names.append(records.split()[0])
			data10.append(records.split()[2:])
		else:
			names.append(records.split()[0])
			data10.append(records.split()[2:][::-1])		


	data10=np.array(data10,dtype=float)
	print(data10.shape)

	center=(data10.shape[1]-1)/2
	start=center-searchTSS_win
	end=center+searchTSS_win
	data10_center_sig=data10[:,range(start,end+1)]

	print(data10_center_sig.shape)

	result1=open(output_file,'w')
	result2=open(output_file+'.'+str(output_win*2+1)+'.bed','w')
	result3=open(output_file+'.'+str(output_win*2+1)+'.procap_pos.bed','w')

	for pk,records in zip(data01,data10_center_sig):
		max_sig=np.amax(records)
		max_id=np.where(records==max_sig)[0]
		if len(max_id)==1:
			new_start=int(pk[1])+start+max_id[0]
		else:
			max_id_closerTSS=np.argmin((max_id-searchTSS_win)**2)
			new_start=int(pk[1])+start+max_id[max_id_closerTSS]

		result1.write(pk[0]+'\t'+str(new_start)+'\t'+str(new_start+1)+'\t'+pk[3]+'\t'+pk[4]+'\t'+pk[5]+'\n')
		result2.write(pk[0]+'\t'+str(new_start-output_win)+'\t'+str(new_start+1+output_win)+'\t'+pk[3]+'\t'+pk[4]+'\t'+pk[5]+'\n')
		if max_sig!=0:
			result3.write(pk[0]+'\t'+str(new_start-output_win)+'\t'+str(new_start+1+output_win)+'\t'+pk[3]+'\t'+pk[4]+'\t'+pk[5]+'\n')

	result1.close()
	result2.close()
	result3.close()
	data0.close()
	data1.close()


############################################################################
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hb:r:o:t:w:")
	except getopt.GetoptError:
		print 'python relocate_TSS.py -b bed_file -r readscount_file -o output_file -t searchTSS_win -w output_win'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python relocate_TSS.py -b bed_file -r readscount_file -o output_file -t searchTSS_win -w output_win'
			sys.exit()
		elif opt=="-b":
			bed_file=str(arg.strip())
		elif opt=="-r":
			readscount_file=str(arg.strip())
		elif opt=="-o":
			output_file=str(arg.strip())
		elif opt=="-t":
			searchTSS_win=int(arg.strip())
		elif opt=="-w":
			output_win=int(arg.strip())
	relocate_TSS(bed_file,readscount_file,output_file,searchTSS_win,output_win)

if __name__=="__main__":
	main(sys.argv[1:])



