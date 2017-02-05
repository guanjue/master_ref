import numpy as np

data0=open('Xu_2009_SUTs_V64.TSS.bed','r')
data00=[]

for records in data0:
	data00.append(records.split())

data01=data00
data02=data00

data1=open('Xu_2009_SUTs_V64_procap_readc_sense.tabular','r')
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
start=center-250
end=center+250
data10_center_sig=data10[:,range(start,end+1)]

print(data10_center_sig.shape)

result1=open('new_TSS.bed','w')
result2=open('new_TSS.1200.bed','w')

for pk,records in zip(data01,data10_center_sig):
	max_sig=np.amax(records)
	max_id=np.where(records==max_sig)[0]
	if len(max_id)==1:
		new_start=int(pk[1])+start+max_id[0]
	else:
		max_id_closerTSS=np.argmin((max_id-250)**2)
		new_start=int(pk[1])+start+max_id[max_id_closerTSS]

	result1.write(pk[0]+'\t'+str(new_start)+'\t'+str(new_start+1)+'\t'+pk[3]+'\t'+pk[4]+'\t'+pk[5]+'\n')
	result2.write(pk[0]+'\t'+str(new_start-600)+'\t'+str(new_start+1+600)+'\t'+pk[3]+'\t'+pk[4]+'\t'+pk[5]+'\n')

result1.close()
result2.close()
data0.close()
data1.close()




