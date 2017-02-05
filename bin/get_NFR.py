import numpy as np
from scipy import signal

def normal_array(width, sigma, normalize=True):
	''' Returns an array of the normal distribution of the specified width '''
	import numpy as np
	import math
	sigma2 = float(sigma)**2
	def normal_func(x):
		return math.exp( -x * x / ( 2 * sigma2 ))
	# width is the half of the distribution
	values = map( normal_func, range(-width, width+1) )
	values = np.array( values, np.float )
	# normalization
	if normalize:
		values = 1.0/math.sqrt(2 * np.pi * sigma2) * values 
	return values

def get_nfr(cdt_file,win_size,sigma,minimum_nfr,maximum_nfr,minimum_dist):
	data0=open(cdt_file,'r')
	peak_length = len(data0.readline().split()[2:])

	data01=[]
	### smooth the tag
	win=win_size
	one_normaldist=normal_array(win, sigma)
	### get smoothed data at each position (weighted norm at each position)
	for records in data0:
		tag_tmp0=np.float32(records.split()[2:])
		tag_tmp01=( np.repeat([one_normaldist],tag_tmp0.shape[0],axis=0).transpose() * tag_tmp0 ).transpose()
		tag_tmp001=tag_tmp01[0]
		for s01 in range(1,tag_tmp01.shape[0]):
			### use norm unit to smooth np array
			tag_tmp001 = np.append(tag_tmp001,[0.0]) + np.append([0.0]*s01,tag_tmp01[s01])
		### remove padding
		tag_tmp001 = tag_tmp001[win:tag_tmp01.shape[0]+win]
		tag_f=signal.medfilt(tag_tmp001,kernel_size=15)

		data01.append(tag_f)
	data01=np.array(data01,dtype=float)

	gradient_pk=np.gradient(data01)
	gradient_pk=np.array(gradient_pk)

	print('write result!')

	np.savetxt(cdt_file+'.smooth.txt',data01,fmt='%10.5f',delimiter='\t')
	np.savetxt(cdt_file+'.grad1.txt',gradient_pk[1],fmt='%10.5f',delimiter='\t')

	data0.close()
	
	grad = gradient_pk[1]
	grad_zero = np.where(grad==0)
	grad_zero=np.array(grad_zero)

	###############################
	data_grad = gradient_pk[1]
	data_grad_0s = np.array(np.where(data_grad==0)) ### find local max/min pos

	### convert local max/min pos to dict
	data_grad_0s_dict = {}
	for rec in np.transpose(data_grad_0s):
		if rec[0] in data_grad_0s_dict:
			data_grad_0s_dict[rec[0]].append(rec[1])
		else:
			data_grad_0s_dict[rec[0]]=[rec[1]]

	### keep only the boundary of the local max/min pos
	data_grad_0s_boundary_dict={}
	for records in data_grad_0s_dict:
		data_grad_0s_boundary_dict[records] = []
		tmp=[]
		for rec in data_grad_0s_dict[records]:
			### put continuous 0s pos into tmp
			if tmp==[]:
				tmp.append(rec)
			elif rec-1 in tmp:
				tmp.append(rec)
			else:
				### keep the start and end OR empty
				if len(tmp) >2:
					data_grad_0s_boundary_dict[records].append([tmp[0],tmp[-1]])
				else:
					data_grad_0s_boundary_dict[records].append(tmp)
				tmp=[]
				tmp.append(rec)
		### keep the start and end OR empty LAST
		if len(tmp) >2:
			data_grad_0s_boundary_dict[records].append([tmp[0],tmp[-1]])
		else:
			data_grad_0s_boundary_dict[records].append(tmp)

	### distinguish local maximun (1) and local minimum (-1)
	data_grad_0s_LM_dict={}
	for records in data_grad_0s_boundary_dict:
		data_grad_0s_LM_dict[records]=[]
		for rec in data_grad_0s_boundary_dict[records]:
			if rec[0] >= 1 and rec[-1] <= data_grad.shape[1]-2: ### check all pos but peak ends
				if data_grad[records,rec[0]-1] > 0 and data_grad[records,rec[-1]+1] < 0:
					data_grad_0s_LM_dict[records].append(1)
				elif data_grad[records,rec[0]-1] < 0 and data_grad[records,rec[-1]+1] > 0:
					data_grad_0s_LM_dict[records].append(-1)
			elif rec[0] == 0: ### check peak start
				if data_grad[records,rec[0]+1] < 0:
					data_grad_0s_LM_dict[records].append(1)
				elif data_grad[records,rec[0]+1] > 0:
					data_grad_0s_LM_dict[records].append(-1)			
			elif rec[-1] == data_grad.shape[1]-1: ### check peak end
				if data_grad[records,rec[-1]-1] > 0:
					data_grad_0s_LM_dict[records].append(1)
				elif data_grad[records,rec[-1]-1] < 0:
					data_grad_0s_LM_dict[records].append(-1)

	### extract all local minimum matrix
	data_grad_local_LM_matrix_dict={}
	for records in data_grad_0s_LM_dict:
		data_grad_local_LM_matrix_dict[records]=[]
		for i in range(1,len(data_grad_0s_LM_dict[records])-1):
			### if minus peak label = -1; upstream of TSS; length of NFR > minimum_nfr; length of NFR > maximum_nfr
			if data_grad_0s_LM_dict[records][i] == -1 and (data_grad_0s_boundary_dict[records][i][0]+data_grad_0s_boundary_dict[records][i][-1])/2<=int(peak_length/2)-minimum_dist and data_grad_0s_boundary_dict[records][i+1][0]-data_grad_0s_boundary_dict[records][i-1][-1] >=minimum_nfr and data_grad_0s_boundary_dict[records][i+1][0]-data_grad_0s_boundary_dict[records][i-1][-1] <= maximum_nfr:
				data_grad_local_LM_matrix_dict[records].append( [ data_grad_0s_boundary_dict[records][i-1], data_grad_0s_boundary_dict[records][i], data_grad_0s_boundary_dict[records][i+1], data_grad_0s_boundary_dict[records][i+1][0]-data_grad_0s_boundary_dict[records][i-1][-1] ] )

	### extract minus1_nfr_plus1_dist matrix
	data_grad_local_minus1_nfr_plus1_dist_matrix_dict={}
	for records in data_grad_local_LM_matrix_dict:
		if len(data_grad_local_LM_matrix_dict[records]) > 0:
			data_grad_local_minus1_nfr_plus1_dist_matrix_dict[records]=data_grad_local_LM_matrix_dict[records][-1]
		else:
			data_grad_local_minus1_nfr_plus1_dist_matrix_dict[records]=[]

	### write results
	nfr_output=open(cdt_file+'.nucleosome.txt','w')
	no_sig = 0
	for i in range(0,len(data_grad_local_minus1_nfr_plus1_dist_matrix_dict)):
		if len(data_grad_local_minus1_nfr_plus1_dist_matrix_dict[i]) > 0:
			nfr_output.write(str((data_grad_local_minus1_nfr_plus1_dist_matrix_dict[i][0][0]+data_grad_local_minus1_nfr_plus1_dist_matrix_dict[i][0][-1])/2)+'\t') ### minus one nucleosome
			nfr_output.write(str((data_grad_local_minus1_nfr_plus1_dist_matrix_dict[i][1][0]+data_grad_local_minus1_nfr_plus1_dist_matrix_dict[i][1][-1])/2)+'\t') ### nfr
			nfr_output.write(str((data_grad_local_minus1_nfr_plus1_dist_matrix_dict[i][2][0]+data_grad_local_minus1_nfr_plus1_dist_matrix_dict[i][2][-1])/2)+'\t') ### plus one nucleosome
			nfr_output.write(str(data_grad_local_minus1_nfr_plus1_dist_matrix_dict[i][3])+'\n')
		else:
			nfr_output.write('na'+'\t'+'na'+'\t'+'na'+'\t'+'na'+'\n')
			no_sig = no_sig+1
	nfr_output.close()
	print(no_sig)

############################################################################
### python get_NFR.py -f SGD_features_ORF_Verified_modified_h3_merged_MN_seq_sort_read1_combined.cdt -w 100 -s 25 -l 150 -u 500 -d 100
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:f:w:s:l:u:d:")
	except getopt.GetoptError:
		print 'python get_NFR.py -f <cdt file> -w <smooth window size> -s <Gaussian smooth sigma> -l <NFR minimum length> -u <NFR maximum length> -d <minimum dist to TSS>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python get_NFR.py -f <cdt file> -w <smooth window size> -s <Gaussian smooth sigma> -l <NFR minimum length> -u <NFR maximum length> -d <minimum dist to TSS>'
			sys.exit()
		elif opt=="-f":
			cdt_file=str(arg.strip())
		elif opt=="-w":
			win_size=int(arg.strip())
		elif opt=="-s":
			sigma=int(arg.strip())
		elif opt=="-l":
			minimum_nfr=int(arg.strip())
		elif opt=="-u":
			maximum_nfr=int(arg.strip())
		elif opt=="-d":
			minimum_dist=int(arg.strip())

	get_nfr(cdt_file,win_size,sigma,minimum_nfr,maximum_nfr,minimum_dist)
if __name__=="__main__":
	main(sys.argv[1:])




