def dist_sort(bed1,bed2,cdt1,cdt2):
	import numpy as np

	data01=open(bed1,'r')
	bed1start_strand=[]
	strand=[]
	for records in data01:
		bed1start_strand.append([int(records.split()[1]),records.split()[5],records.split()[3]])

	data02=open(bed2,'r')
	bed2start=[]
	for records in data02:
		bed2start.append([int(records.split()[1]),records.split()[3]])

	dist=[]
	label=[]
	for b1_s, b2 in zip(bed1start_strand, bed2start):
		if b1_s[1]=='+':
			dist.append([b1_s[0]-b2[0], b1_s[2], b2[1]])
		else:
			dist.append([-b1_s[0]+b2[0], b1_s[2], b2[1]])
	dist=np.array(dist)

	dist_sort=dist[np.argsort(dist[:, 0].astype(np.float))]
	dist_sort0=dist_sort
	dist_result=open(bed2+'.dist.txt','w')
	for records in dist_sort0:
		dist_result.write(str(records[0])+'\n')

	dist_result.close()



	dist1=dist_sort
	data11=open(cdt1,'r')
	cdt_1={}
	cdt_1_sort=open(cdt1+'.distsort.txt','w')
	#cdt_1_sort.write(data11.readline())
	#cdt_1_sort.write(data11.readline())
	for records in data11:
		cdt_1[records.split()[0]]=records

	for records in dist1:
		cdt_1_sort.write(cdt_1[records[1]])


	dist2=dist_sort
	data21=open(cdt2,'r')
	cdt_2={}
	cdt_2_sort=open(cdt2+'.distsort.txt','w')
	#cdt_2_sort.write(data21.readline())
	#cdt_2_sort.write(data21.readline())
	for records in data21:
		cdt_2[records.split()[0]]=records

	for records in dist2:
		cdt_2_sort.write(cdt_2[records[1]])

	data01.close()
	data02.close()
	data11.close()
	data21.close()
	cdt_1_sort.close()
	cdt_2_sort.close()


############################################################################
### python dist_sort.py -a Xu_2009_ORF-Ts_V64.procap.TSS.bed.1200 -b Xu_2009_ORF-Ts_V64.TSS -c Xu_2009_ORF-Ts_V64_procap_readc_sense.newTSS_K_G1 -d Xu_2009_ORF-Ts_V64_procap_readc_sense_K_G1
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:a:b:c:d:")
	except getopt.GetoptError:
		print 'dist_sort.py -a <1 TSS bed file> -b <2 TSS bed file> -c <1 TSS cdt file> -d <2 TSS cdt file>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'dist_sort.py -a <1 TSS bed file> -b <2 TSS bed file> -c <1 TSS cdt file> -d <2 TSS cdt file>'
			sys.exit()
		elif opt=="-a":
			bed1=str(arg.strip())
		elif opt=="-b":
			bed2=str(arg.strip())
		elif opt=="-c":
			cdt1=str(arg.strip())
		elif opt=="-d":
			cdt2=str(arg.strip())

	dist_sort(bed1,bed2,cdt1,cdt2)
if __name__=="__main__":
	main(sys.argv[1:])
