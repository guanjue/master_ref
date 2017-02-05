import numpy as np
from rdp import rdp
import matplotlib.pyplot as plt
from scipy import signal
from sklearn import mixture
import matplotlib.mlab as mlb
from scipy import stats

#data0s=open('Xu_2009_ORF-Ts_V64_SRR1951311_readc_sense_K_A1.cdt','r')
def GMM(input_file_sense,input_file_anti,scale_n):
	input_file=input_file_sense
	data0s=open(input_file_sense+'.cdt','r')
	data0s.readline()
	data0s.readline()
	data0s1=[]
	data0s1_n=[]
	for records in data0s:
		data0s1_n.append(records.split()[0])
		data0s1.append(records.split()[3:])
	data0s1=np.array(data0s1,dtype=float)

	data1s=open(input_file_anti+'.cdt','r')
	data1s.readline()
	data1s.readline()
	data1s1=[]
	data1s1_n=[]
	for records in data1s:
		data1s1_n.append(records.split()[0])
		data1s1.append(records.split()[3:])
	data1s1=np.array(data1s1,dtype=float)

	data0s1=data0s1+data1s1

	data0s1_sum=np.sum(data0s1,axis=0) ### get the sum of each position
	data0s1_sum=data0s1_sum#/np.sum(data0s1_sum)*1000000 ### normalize the the same total number of reads count in the pks

	print(data0s1_sum)
	print(data0s1_sum.shape)

	data0s1_sum_smooth=signal.medfilt(data0s1_sum,kernel_size=11) ### smooth by the median number within a window
	uplimit=np.percentile(data0s1_sum_smooth[np.nonzero(data0s1_sum_smooth)], 100)
	data0s1_sum_smooth=np.clip(data0s1_sum_smooth,0, uplimit )

	plt.plot(range(-600,601),data0s1_sum_smooth,'b-')
	plt.xlim(-600,600)
	plt.grid(True)
	plt.savefig(input_file+'.composite.curve.pdf')
	plt.close()

	#scale_n=data0s1.shape[0]/20

	data0s1_sum_bin=np.repeat(np.sum(np.split(data0s1_sum_smooth[0:data0s1_sum.shape[0]-1,],240),axis=1),5) ### bin the reads count
	data0s1_sum_bin=data0s1_sum_bin[0:data0s1_sum.shape[0]-1,]
	#scale_n=np.sum(data0s1_sum_bin)/1000
	print(scale_n)
	data0s1_sum_bin=data0s1_sum_bin/scale_n*5

	#data0s1_sum_bin=signal.medfilt(data0s1_sum_bin,kernel_size=11)
	plt.plot(range(-600,600),data0s1_sum_bin,'b-')
	plt.xlim(-600,600)
	plt.grid(True)
	plt.savefig(input_file+'.composite.bin.curve.pdf')
	plt.close()

	### For GMM model
	### Convert data to GMM input format
	data0s1_sum_smooth=data0s1_sum_smooth/scale_n*2
	yourdata=[]
	for i in range(0,len(data0s1_sum_bin)):
		yourdata=np.append(yourdata,np.array([i+1]*int(round(data0s1_sum_smooth[i]))) )
	yourdata=(np.transpose([yourdata]))

	print(yourdata.shape)
	### fit GMM
	gmm = mixture.GaussianMixture(n_components=10, covariance_type='full', max_iter=500, tol=0.0000)
	gmm.fit(yourdata)
	clf = gmm

	means=[]
	sds=[]
	cs=[]
	ws=[]

	x=np.array(range(0,len(data0s1_sum_bin)))

	histdist = plt.hist(yourdata, 240, normed=True, histtype='stepfilled')

	y_real_sec=data0s1_sum_bin
	y_predict=np.repeat(0.0,len(data0s1_sum_bin))
	y_predict_derivative=np.repeat(0.0,len(data0s1_sum_bin))

	for i in range(0,len(clf.means_)):
		sd_tmp=np.sqrt(clf.covariances_[i])[0]# / clf.weights_[i] #np.sqrt()
		if (1==1) :#sd_tmp < sds[0]:
			means.append(clf.means_[i])
			#means=[clf.means_[i]]
			ws.append(clf.weights_[i])
			cs.append(clf.covariances_[i])
			sds.append(sd_tmp)
			#sds=[sd_tmp]
		plotgauss1 = lambda x: plt.plot(x,(clf.weights_[i])*mlb.normpdf(x,clf.means_[i],np.sqrt(clf.covariances_[i]))[0], linewidth=1.0,color='0.75',linestyle='--')
		#print(plotgauss1)
		plotgauss1(histdist[1])
		y_real_sec=y_real_sec-np.sum(data0s1_sum_bin)*(clf.weights_[i])*mlb.normpdf(x,clf.means_[i],np.sqrt(clf.covariances_[i]))[0]
		y_predict=y_predict+(clf.weights_[i])*mlb.normpdf(x,clf.means_[i],np.sqrt(clf.covariances_[i]))[0]#*np.sum(data0s1_sum)
		#print(np.exp(-(np.array(x)-clf.means_[i])/2/clf.covariances_[i]))
		y_predict_derivative=y_predict_derivative+ ( clf.weights_[i] * 1/np.sqrt(2*3.14*clf.covariances_[i]) *  np.exp(-(np.array(x)-clf.means_[i])**2/2/clf.covariances_[i]) * (-2/2/clf.covariances_[i] * (x-clf.means_[i])) )

	y_predict_derivative=(y_predict_derivative)[0][200:-200]

	derivative_max=np.max(y_predict_derivative)
	derivative_max_id=np.argmax(y_predict_derivative)
	derivative_min=np.min(y_predict_derivative)
	derivative_min_id=np.argmin(y_predict_derivative)

	plt.plot(x, y_predict, 'r--', linewidth=1.5)
	weight_p=(clf.weights_) / np.sum(clf.weights_)
	mean_all=np.sum( (weight_p) * np.transpose(clf.means_) )
	WSD=np.sqrt(np.sum(abs(yourdata - stats.mode(yourdata)[0][0][0])**2)/(len(yourdata)-1))
	#variance = (weight_p * len(yourdata) / (weight_p * len(yourdata)-1)) * np.transpose(clf.covariances_)

	std_within_components=np.sqrt(np.sum( weight_p * ( np.transpose(clf.covariances_) ) ) )
	#std_all=np.sqrt(  np.sum( weight_p * ( (np.transpose(clf.means_))**2 + variance ) ) - mean_all**2  )
	#std_outside=np.sqrt(  np.sum( weight_p * ( (np.transpose(clf.means_))**2  ) ) - mean_all**2  )
	#std_all_modify=np.sqrt(np.sum( weight_p * ( variance ) ) )
	print('Highest Pos: ')
	max_indice_f=np.argmax(data0s1_sum_bin,axis=0)

	print(max_indice_f)
	print(WSD)
	print(std_within_components)
	#std_all_newref=np.sqrt(np.sum( weight_p * ( np.transpose(clf.means_)**2 + np.transpose(clf.covariances_) ) ) - max_indice_f**2)
	#plt.ylim(0,0.015)
	plt.savefig(input_file+'.fit_curve.pdf')
	plt.close()
	'''
	plt.plot(x[100:-100], y_predict_derivative, 'k-', linewidth=2)
	plt.ylim(-0.00008,0.00008)
	plt.savefig(input_file+'.fit_curve.derivative.pdf')
	plt.close()
	'''
	return WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components

data_type=['other_transcripts']#'CUTS','ORF-Ts','SUTS','other_transcripts']#,'All']

for dt in data_type:
	print(dt)
	WSD0=[]
	max_indice_f0=[]
	derivative_max0=[]
	derivative_max_id0=[]
	derivative_min0=[]
	derivative_min_id0=[]
	std_within_components0=[]

	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_sua7_readc_sense_K_G1','Xu_2009_'+dt+'_V64_sua7_readc_anti_K_G1',500)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)
	
	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_sua7_readc_sense.newTSS_K_G1','Xu_2009_'+dt+'_V64_sua7_readc_anti.newTSS_K_G1',500)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)

	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_proseq_readc_anti_K_G1','Xu_2009_'+dt+'_V64_proseq_readc_anti_K_G1',200)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)
	
	
	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_proseq_readc_anti.newTSS_K_G1','Xu_2009_'+dt+'_V64_proseq_readc_anti.newTSS_K_G1',200)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)


	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_sua7_pipseq_readc_sense_K_G1','Xu_2009_'+dt+'_V64_sua7_pipseq_readc_anti_K_G1',500)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)
	
	
	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_sua7_pipseq_readc_sense.newTSS_K_G1','Xu_2009_'+dt+'_V64_sua7_pipseq_readc_anti.newTSS_K_G1',500)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)

	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_h3_mnseq_readc_combined_K_G1','Xu_2009_'+dt+'_V64_h3_mnseq_readc_combined_K_G1',10)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)
	
	
	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_h3_mnseq_readc_combined.newTSS_K_G1','Xu_2009_'+dt+'_V64_h3_mnseq_readc_combined.newTSS_K_G1',10)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)

	'''
	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_procap_readc_sense_K_G1','Xu_2009_'+dt+'_V64_procap_readc_anti_K_G1',500)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)
	
	
	WSD, max_indice_f, derivative_max, derivative_max_id, derivative_min, derivative_min_id, std_within_components=GMM('Xu_2009_'+dt+'_V64_procap_readc_sense.newTSS_K_G1','Xu_2009_'+dt+'_V64_procap_readc_sense.newTSS_K_G1',500)
	WSD0.append(WSD)
	max_indice_f0.append(max_indice_f)
	derivative_max0.append(derivative_max)
	derivative_max_id0.append(derivative_max_id)
	derivative_min0.append(derivative_min)
	derivative_min_id0.append(derivative_min_id)
	std_within_components0.append(std_within_components)
	'''

	TF_name=['Xu_2009','Booth_2016']
	TF_name=['Xu_2009.Sua7','Booth_2016.Sua7','Xu_2009.Pro-seq','Booth_2016.Pro-seq','Xu_2009.Pip-seq','Booth_2016.Pip-seq','Xu_2009.MN-seq','Booth_2016.MN-seq']

	col=['bp','rp','mp','gp','yp','kp','wp','cp','bp','rp','mp','gp','yp','kp','wp','cp']
	fig=plt.figure()#facecolor="white")
	ax = fig.add_subplot(111)
	max_pos_i=0
	tf_label0=[]
	for x,y in zip(np.array(max_indice_f0)-500, WSD0):
		tf_label, =ax.plot(x,y,col[max_pos_i])#,label=TF_name[max_pos_i])
		tf_label0.append(tf_label)
		max_pos_i=max_pos_i+1

	plt.legend(tf_label0, TF_name,loc=3,fontsize=10)
	#print(len(np.sum(data03_f_n_WSD,axis=1)))
	#plt.plot(np.sum(data03_f_n_WSD,axis=0), 'bp', -np.sum(data03_r_n_WSD,axis=0), 'rp')
	plt.xlim(-600,601)
	#plt.ylim(0, 300)
	plt.grid(True)
	plt.axhline(0, color='black', lw=2)
	plt.savefig(dt+'.max_pos.WSD_all.pdf')
	plt.close()




	fig=plt.figure()#facecolor="white")
	ax = fig.add_subplot(111)
	max_pos_i=0
	tf_label0=[]
	for x,y in zip(np.array(max_indice_f0)-500, std_within_components0):
		tf_label, =ax.plot(x,y,col[max_pos_i])#,label=TF_name[max_pos_i])
		tf_label0.append(tf_label)
		max_pos_i=max_pos_i+1

	plt.legend(tf_label0, TF_name,loc=3,fontsize=10)
	#print(len(np.sum(data03_f_n_WSD,axis=1)))
	#plt.plot(np.sum(data03_f_n_WSD,axis=0), 'bp', -np.sum(data03_r_n_WSD,axis=0), 'rp')
	plt.xlim(-600,601)
	#plt.ylim(0, 300)
	plt.grid(True)
	plt.axhline(0, color='black', lw=2)
	plt.savefig(dt+'.max_pos.std_within_components0.pdf')
	plt.close()


	fig=plt.figure()#facecolor="white")
	ax = fig.add_subplot(111)
	max_pos_i=0
	tf_label0=[]
	for x,y in zip(np.array(derivative_max_id0)-400, derivative_max0):
		tf_label, =ax.plot(x,y,col[max_pos_i])#,label=TF_name[max_pos_i])
		tf_label0.append(tf_label)
		max_pos_i=max_pos_i+1
	max_pos_i=0
	for x,y in zip(np.array(derivative_min_id0)-400, derivative_min0):
		tf_label, =ax.plot(x,y,col[max_pos_i])#,label=TF_name[max_pos_i])
		#tf_label0.append(tf_label)
		max_pos_i=max_pos_i+1

	plt.legend(tf_label0, TF_name,loc=3,fontsize=10)
	#print(len(np.sum(data03_f_n_WSD,axis=1)))
	#plt.plot(np.sum(data03_f_n_WSD,axis=0), 'bp', -np.sum(data03_r_n_WSD,axis=0), 'rp')
	plt.xlim(-600,601)
	#plt.ylim(-0.00008, 0.00008)
	plt.grid(True)
	plt.axhline(0, color='black', lw=2)
	plt.savefig(dt+'.GMM_derivative.pdf')
	plt.close()


'''
print('RDP smoothing epsilon!!!!!!')
rdp_spsilon=np.sum(data0s1_sum)/len(data0s1_sum)
print(rdp_spsilon)
data0s1_sum_rdp=rdp(np.transpose(np.array([np.array(range(0,len(data0s1_sum))), data0s1_sum])),epsilon=rdp_spsilon)
print(data0s1_sum_rdp.shape)
print(data0s1_sum_rdp)

print('RDP smoothing!!!!!!!!')

rdp_slope=[]
for i in range(0,data0s1_sum_rdp.shape[0]-1):
	rdp_slope.append( [(data0s1_sum_rdp[i+1,0]+data0s1_sum_rdp[i,0])/2 , (data0s1_sum_rdp[i+1,1]-data0s1_sum_rdp[i,1]) / (data0s1_sum_rdp[i+1,0]-data0s1_sum_rdp[i,0])] )
	#rdp_slope.append( [y_real_rdp[i,0] , (y_real_rdp[i+1,1]-y_real_rdp[i,1]) / (y_real_rdp[i+1,0]-y_real_rdp[i,0])] )

rdp_slope=np.array(rdp_slope)
#print(rdp_slope)

#print(np.gradient(y_real_rdp))
#plt.plot(np.transpose(data0s1_sum_rdp)[0]-600,np.transpose(data0s1_sum_rdp)[1],'b-')
plt.plot(np.transpose(data0s1_sum)[0]-600,np.transpose(data0s1_sum)[1],'b-')

plt.xlim(-700,700)
plt.grid(True)
plt.savefig('test'+'.rdp.curve.pdf')
plt.close()
'''

