def compare_GTF(ingtf,refgtf,outfile,match_file):
	# Compare two gtf files
	yui=open(refgtf,'r')
	ref_dict=dict()
	current=[]
	trans=[]
	strand=''
	name=''
	tcount=0
	ecount=0
	last_chro=''
	for line in yui:
		line=line.replace('chr','')
		split=line.split()
		if (len(split)<4) or (line[0]=='#'):
			continue
		if split[0] not in ref_dict:
			ref_dict[split[0]]=[]
		if split[2]=='transcript':
			if name!='':
				exon_list.sort()
				ref_dict[last_chro].append([trans,strand,exon_list,name,length,splice])
			tcount+=1
			trans=[int(split[3]),int(split[4])]
			trans.sort()
			strand=split[6]
			for i in range(8,len(split),1):
				if split[i]=='transcript_id':
					name=split[i+1][1:-2]
					break
			exon_list=[]
			length=0
			splice=0
		elif split[2]=='exon':
			ecount+=1
			coord=[int(split[3]),int(split[4])]
			splice+=2
			coord.sort()
			length+=(coord[1]-coord[0])+1
			exon_list.append(coord)
		last_chro=split[0]
	yui.close()
	print 'Number of transcripts in reference:',tcount
	print 'Number of exons in reference:',ecount
	if name!='':
		ref_dict[last_chro].append([trans,strand,exon_list,name,length,splice])
	for chromo in ref_dict:
		ref_dict[chromo].sort()
	yui=open(ingtf,'r')
	new=open(outfile,'a')
	newm=open(match_file,'a')
	name=''
	last_chro=''
	tcount=0
	ecount=0
	written=0
	unique=0
	match=0
	exonlap=0
	noverl=0
	tlen=0
	trans=[]
	for line in yui:
		split=line.split()
		if (len(split)<4) or (line[0]=='#'):
			continue
		found=0
		if split[2]=='transcript':
			tcount+=1
			if name!='':
				if last_chro not in ref_dict:
					continue
				exon_list.sort()
				max_exon_splice=[]
				for one in ref_dict[last_chro]:
					if one[0][0]>trans[1]:
						break
					elif one[0][1]<trans[0]:
						continue
					elif strand!=one[1]:
						continue
					else:
						found=1
						nmatch=0
						splicem=0
						ref_exon_list=one[2]
						overl_ec=0
						sm=[]
						for inexon in exon_list:
							syes,eyes=0,0
							for refexon in ref_exon_list:
								if (refexon[0]-10)<=inexon[0]<=(refexon[0]+10): # splices are the same if found within 10bp
									splicem+=1
									syes=1
								if (refexon[1]-10)<=inexon[1]<=(refexon[1]+10):
									splicem+=1
									eyes=1
								if inexon[1]<refexon[0]:
									break
								elif inexon[0]>refexon[1]:
									continue
								else:
									if exon_list.index(inexon)==ref_exon_list.index(refexon):
										overl_ec+=1
									start=max([inexon[0],refexon[0]])
									end=min([inexon[1],refexon[1]])
									nmatch+=(end-start+1)	
							sm.extend([syes,eyes])					
						if overl_ec==0:
							continue
						new.write(last_chro+'	'+str(one[0][0])+'	'+str(one[0][1])+'	'+str(trans[0])+'	'+str(trans[1])+'	'+strand+'	'+name+'	'+one[3]+'	'+str(float(overl_ec)*2*100/(len(ref_exon_list)+len(exon_list)))+'	'+str(float(nmatch)*100/one[4])+'	'+str(float(splicem)*100/one[5])+'\n')
						written+=1
						max_exon_splice.append([(float(nmatch)*100/max([tlen,one[4]]))+(float(splicem)*100/(2*max([len(ref_exon_list),len(exon_list)]))),(float(nmatch)*100/max([tlen,one[4]])),(float(splicem)*100/(2*max([len(ref_exon_list),len(exon_list)]))),(float(overl_ec)*100/max([len(ref_exon_list),len(exon_list)])),one[3],sm,exon_list,ref_exon_list])
						# exon_match_perc, splice_match_perc and perc_overl_exon calculated here
				if len(max_exon_splice)!=0:
					max_exon_splice.sort(reverse=True)
					if len(max_exon_splice[0][-1])!=len(max_exon_splice[0][-2]) or max_exon_splice[0][3]<100:
						decision='different'
					elif max_exon_splice[0][1]>=75:
						decision='same'
					else:
						decision='different'
					if max_exon_splice[0][1]>0:
						exonlap+=1
					newm.write(name+'	'+max_exon_splice[0][4]+'	'+str(max_exon_splice[0][3])+'	'+str(max_exon_splice[0][1])+'	'+str(max_exon_splice[0][2])+'	'+decision+'\n')
					match+=1
				else:
					noverl+=1
				if found==1:
					unique+=1
			trans=[int(split[3]),int(split[4])]
			trans.sort()
			strand=split[6]
			for i in range(8,len(split),1):
				if split[i]=='transcript_id':
					name=split[i+1][1:-2]
					break
			exon_list=[]
			tlen=0
		elif split[2]=='exon':
			ecount+=1
			coord=[int(split[3]),int(split[4])]
			coord.sort()
			tlen+=coord[1]-coord[0]+1
			exon_list.append(coord)
		last_chro=split[0]
	if (name!='') and (last_chro in ref_dict):
		exon_list.sort()
		for one in ref_dict[last_chro]:
			if one[0][0]>trans[1]:
				break
			elif one[0][1]<trans[0]:
				continue
			elif strand!=one[1]:
				continue
			else:
				found=1
				nmatch=0
				splicem=0
				ref_exon_list=one[2]
				overl_ec=0
				sm=[]
				for inexon in exon_list:
					syes,eyes=0,0
					for refexon in ref_exon_list:
						if (refexon[0]-10)<=inexon[0]<=(refexon[0]+10): # splices are the same if found within 10bp
							splicem+=1
							syes=1
						if (refexon[1]-10)<=inexon[1]<=(refexon[1]+10):
							splicem+=1
							eyes=1
						if inexon[1]<refexon[0]:
							break
						elif inexon[0]>refexon[1]:
							continue
						else:
							overl_ec+=1
							start=max([inexon[0],refexon[0]])
							end=min([inexon[1],refexon[1]])
							nmatch+=(end-start+1)
					sm.extend([syes,eyes])
				if overl_ec==0:
					continue
				new.write(last_chro+'	'+str(one[0][0])+'	'+str(one[0][1])+'	'+str(trans[0])+'	'+str(trans[1])+'	'+strand+'	'+name+'	'+one[3]+'	'+str(float(overl_ec)*2*100/(len(ref_exon_list)+len(exon_list)))+'	'+str(float(nmatch)*100/one[4])+'	'+str(float(splicem)*100/one[5])+'\n')
				written+=1
				max_exon_splice.append([(float(nmatch)*100/max([tlen,one[4]]))+(float(splicem)*100/(2*max([len(ref_exon_list),len(exon_list)]))),(float(nmatch)*100/max([tlen,one[4]])),(float(splicem)*100/(2*max([len(ref_exon_list),len(exon_list)]))),(float(overl_ec)*100/max([len(ref_exon_list),len(exon_list)])),one[3],sm,exon_list,ref_exon_list])
				# exon_match_perc, splice_match_perc and perc_overl_exon calculated here
		if len(max_exon_splice)!=0:
			max_exon_splice.sort(reverse=True)
			if len(max_exon_splice[0][-1])!=len(max_exon_splice[0][-2]) or max_exon_splice[0][3]<100:
				decision='different'
			elif max_exon_splice[0][1]>=75:
				decision='same'
			else:
				decision='different'
			if max_exon_splice[0][1]>0:
				exonlap+=1
			newm.write(name+'	'+max_exon_splice[0][4]+'	'+str(max_exon_splice[0][3])+'	'+str(max_exon_splice[0][1])+'	'+str(max_exon_splice[0][2])+'	'+decision+'\n')
			match+=1
		else:
			noverl+=1
		if found==1:
			unique+=1
	yui.close()
	new.close()
	newm.close()
	print 'Number of trancripts:',tcount
	print 'Number of exons:',ecount
	print 'Transcripts with no overlap:',noverl
	print 'Number of pairs written:',written
	print 'Number of matched trancripts:',unique
	print 'Number of matched overlap:', match
	print 'Number of exon-overlap transcripts:',exonlap
	print 'Number of unmatched trancripts:',tcount-unique
	print 'Number of nonoverlapping transcripts:',noverl


def count_threshold(infile,threshold):
	yui=open(infile,'r')
	exon_dict=dict()
	ec_dict=dict()
	splice_dict=dict()
	all_dict=dict()
	total=0
	for line in yui:
		total+=1
		split=line.split()
		ec=float(split[2])
		exon=float(split[3])
		splice=float(split[4])
		if ec>=threshold:
			ec_dict[split[1]]=0
		if exon>=threshold:
			exon_dict[split[1]]=0
		if splice>=threshold:
			splice_dict[split[1]]=0
			if exon>=threshold and ec>=threshold:
				all_dict[split[1]]=0
	yui.close()
	print 'Total:',total
	print 'Exon count threshold pass:',len(ec_dict),float(len(ec_dict))/total
	print 'Exon length threshold pass:',len(exon_dict),float(len(exon_dict))/total
	print 'Splice threshold pass:',len(splice_dict),float(len(splice_dict))/total
	print 'Exon count, length and splice threshold pass:',len(all_dict),float(len(all_dict))/total



from optparse import OptionParser
import os,sys
parser = OptionParser()
parser.add_option('-i','--ingtf', dest='ingtf',help='Absolute path to the input gtf file [required]',type='str',default='')
parser.add_option('-r','--refgtf', dest='refgtf',help='Absolute path to the reference gtf file [required]',type='str',default='')
parser.add_option('-o','--outfile', dest='outfile',help='Absolute path to the output file',type='str',default='')
parser.add_option('-m','--match_file', dest='match_file',help='Absolute path to the match file',type='str',default='')
parser.add_option('-t','--threshold', dest='threshold',help='Threshold completeness for the matched file',type='float',default=100.0)
(options,args)=parser.parse_args()
if (options.ingtf=='') or (options.refgtf=='') or (options.outfile==''):
	print '\nRequired filed(s) not supplied\n# The program compares the input gtf file to the reference gtf file\n'
	parser.print_help()
	sys.exit(1)
if options.match_file=='':
	options.match_file=options.outfile+'.match'

compare_GTF(options.ingtf,options.refgtf,options.outfile,options.match_file)

print 'Getting the summaries...'

print 'Threshold considered:',options.threshold

count_threshold(options.match_file,options.threshold)
