import numpy as np

data_out=open('ensGeneTSS1kwin_SRR3031844_1_readc_bgsub.tabular','w')

data0=open('SRR3031844_1/ensGeneTSS1kwin_SRR3031844_1_readc_sense.tabular','r')
data_out.write(data0.readline())
data00_name=[]
data00_reads=[]
for records in data0:
	tmp=records.split()
	data00_name.append(tmp[0])
	data00_reads.append(tmp[1:])
data00_reads=np.array(data00_reads)

data1=open('SRR3031845_1/ensGeneTSS1kwin_SRR3031845_1_readc_sense.tabular','r')
data1.readline()
data10_reads=[]
for records in data1:
	tmp=records.split()
	data10_reads.append(tmp[1:])
data10_reads=np.array(data10_reads)

dif_reads=data00_reads-data10_reads
dif_reads=dif_reads.clip(min=0)

for name,reads in zip(data00_name,dif_reads):
	data_out.write(name+'\t'+name+'\t')
	for read in reads:
		data_out.write(read+'\t')
	data_out.write('\n')

data_out.close()
data0.close()
data1.close()




# reads processed: 10779040
# reads with at least one reported alignment: 1900389 (17.63%)
# reads that failed to align: 305903 (2.84%)
# reads with alignments suppressed due to -m: 8572748 (79.53%)
Reported 1900389 alignments to 1 output stream(s)

# reads processed: 10779040
# reads with at least one reported alignment: 1894434 (17.58%)
# reads that failed to align: 294435 (2.73%)
# reads with alignments suppressed due to -m: 8590171 (79.69%)
Reported 1894434 alignments to 1 output stream(s)