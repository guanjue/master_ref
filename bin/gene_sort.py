def gene_sort(target_file,target_id_colnum,order_file,order_id_colnum,output_name):
	data0 = open(target_file,'r')
	data01 = {}
	data02 = []

	for records in data0:
		tmp=records.split()
		data01[tmp[target_id_colnum]] = tmp
		data02.append(tmp)

	data0.close()


	data1 = open(order_file,'r')
	data11 = []

	for records in data1:
		data11.append(records.split()[order_id_colnum])

	data1.close()

	result = open(output_name,'w')
	result1 = open('nomatch_'+output_name,'w')
	for records in data11:
		if records in data01:
			for rec in data01[records]:
				result.write(rec+'\t')
			result.write('\n')
	
	for records in data02:
		if not(records[target_id_colnum] in data11):
			for rec in records:
				result1.write(rec+'\t')
			result1.write('\n')			
	
	result.close()
	result1.close()


############################################################################
### python sort_bed.py -t Xu_2009_ORF-Ts_V64_sort_pk_plus1.bed -a 3 -r nfr_based_genename_order_nfrbased.txt -b 0 -o Xu_2009_ORF-Ts_V64_sort_pk_plus1_nfrbased_sort.bed
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:t:a:r:b:o:")
	except getopt.GetoptError:
		print 'python sort_bed.py -t <target_file> -a <target_id_colnum> -r <order_file> -b <order_id_colnum> -o <output_name>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python sort_bed.py -t <target_file> -a <target_id_colnum> -r <order_file> -b <order_id_colnum> -o <output_name>'
			sys.exit()
		elif opt=="-t":
			target_file=str(arg.strip())
		elif opt=="-a":
			target_id_colnum=int(arg.strip())
		elif opt=="-r":
			order_file=str(arg.strip())
		elif opt=="-b":
			order_id_colnum=int(arg.strip())
		elif opt=="-o":
			output_name=str(arg.strip())


	gene_sort(target_file,target_id_colnum,order_file,order_id_colnum,output_name)
if __name__=="__main__":
	main(sys.argv[1:])
