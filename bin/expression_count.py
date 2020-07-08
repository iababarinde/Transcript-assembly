def extract_all2(indir,column,outfile):
	import os
	my_list=os.listdir(indir)
	my_list.sort()
	gene_dict=dict()
	first_sample=my_list[0]
	print 'First sample:',first_sample
	print 'Number of sample:',len(my_list)
	name=''
	for subfile in os.listdir(os.path.join(indir,first_sample)):
		if subfile[-4:]=='.gtf':
			name=subfile
			break
	yui=open(os.path.join(indir,first_sample,name),'r')
	total=0
	nonzero=0
	for line in yui:
		split=line.split()
		if split[2]!='transcript':
			continue
		n=6
		while n<len(split):
			if split[n]=='transcript_id':
				tid=split[n+1][1:-2]
			elif split[n]==column:
				value=str(float(split[n+1][1:-2]))
				break
			n+=1
		gene_dict[tid]={first_sample:value}
		total+=1
		if float(value)!=0:
			nonzero+=1
	print len(gene_dict),'genes found in the first file:',total,nonzero
	for one_file in my_list[1:]:
		print one_file,
		name=''
		for subfile in os.listdir(os.path.join(indir,one_file)):
			if subfile[-4:]=='.gtf':
				name=subfile
				break
		current=open(os.path.join(indir,one_file,name),'r')
		total=0
		nonzero=0
		for line in current:
			split=line.split()
			if split[2]!='transcript':
				continue
			n=6
			while n<len(split):
				if split[n]=='transcript_id':
					tid=split[n+1][1:-2]
				elif split[n]==column:
					value=str(float(split[n+1][1:-2]))
					break
				n+=1
			gene_dict[tid][one_file]=value
			total+=1
			if float(value)!=0:
				nonzero+=1
		print '	'+str(total)+'	'+str(nonzero)
	new=open(outfile,'a')
	new.write('gene_id')
	for name in my_list:
		new.write('	'+name.replace('.bam','').replace('.gtf',''))
	new.write('\n')
	for tid in gene_dict:
		new.write(tid)
		for name in my_list:
			new.write('	'+gene_dict[tid][name])
		new.write('\n')
	new.close()
	print 'Completed!'





from optparse import OptionParser
import os,sys
parser = OptionParser()
parser.add_option('-i','--indir', dest='indir',help='Absolute path to the input gtf directory with values per sample [required]',type='str',default='')
parser.add_option('-c','--column', dest='column',help='The column of the gtf file to be extracted [Options: FPKM(default),cov,TPM]',type='str',default='FPKM')
parser.add_option('-o','--outfile', dest='outfile',help='Absolute path to the output file',type='str',default='')
(options,args)=parser.parse_args()
if (options.indir=='') or (options.outfile==''):
	print '\nRequired filed(s) not supplied\n# The program extracts the expression matrix of the ingtf from indir directory.\n The column to extract can be selected\n'
	parser.print_help()
	sys.exit(1)
extract_all2(options.indir,options.column,options.outfile)