def compare_GTF(ingtf,refgtf,outgtf,outfile,match_file):
	#Specifically for GENCODE reference gtf, use GENCODE gtf to update ingtf
	yui=open(refgtf,'r')
	ref_dict=dict()
	tid_dict=dict()
	current=[]
	trans=[]
	strand=''
	name=''
	tcount=0
	ecount=0
	last_chro=''
	found_exon=0
	for line in yui:
		split=line.split()
		if (len(split)<4) or (line[0]=='#'):
			continue
		split[0]=split[0].replace('chr','')
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
			(gn,tn,egid,etbt)=('','','','')
			for i in range(8,len(split),1):
				if split[i]=='gene_id':
					egid=split[i+1][1:-2]
				if split[i]=='transcript_id':
					name=split[i+1][1:-2]
				if split[i]=='gene_name':
					gn=split[i+1][1:-2]
				if split[i]=='transcript_name':
					tn=split[i+1][1:-2]
				if split[i]=='transcript_type':
					etbt=split[i+1][1:-2]
			tid_dict[name]=(gn,tn,egid,etbt)
			exon_list=[]
			length=0
			splice=0
		elif split[2]=='exon':
			eeid=''
			for i in range(8,len(split),1):
				if split[i]=='exon_id':
					eeid=split[i+1][1:-2]
					found_exon+=1
					break
			ecount+=1
			coord=[int(split[3]),int(split[4])]
			splice+=2
			coord.sort()
			coord.append(eeid)
			length+=(coord[1]-coord[0])+1
			exon_list.append(coord)
		last_chro=split[0]
	yui.close()
	print 'Last exon ID:',eeid
	print 'Number of transcripts in reference:',tcount
	print 'Number of exon ids:',found_exon
	print 'Number of exons in reference:',ecount
	if name!='':
		ref_dict[last_chro].append([trans,strand,exon_list,name,length,splice])
	for chromo in ref_dict:
		ref_dict[chromo].sort()
	yui=open(ingtf,'r')
	new=open(outfile,'a')
	newm=open(match_file,'a')
	newm.write('tid	GENCODE_id	perc_exon_count	perc_overl_exon_cov	perc_overl_splice_count	decision\n')
	new.write('chro	start	end	GENCODE_start	GENCODE_end	strand	tid	GENCODE_id	perc_exon_count	perc_overl_exon_cov	perc_overl_splice_count	decision\n')
	name=''
	last_chro=''
	tcount=0
	ecount=0
	written=0
	unique=0
	match=0
	tmatch=0
	bmatch=0
	tmatch_dict=dict()
	ematch_dict=dict()
	ecount_dict=dict()
	for line in yui:
		split=line.split()
		if len(split[0])>2 or split[0]=='MT':
			continue
		if (len(split)<4) or (line[0]=='#'):
			continue
		found=0
		if split[2]=='transcript':
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
						pos=0
						oec=0
						for inexon in exon_list:
							overlexon=0
							pos+=1
							exonid=name+'e'+str(pos)
							for refexon in ref_exon_list:
								if inexon[1]<refexon[0]:
									break
								elif inexon[0]>refexon[1]:
									continue
								else:
									overlexon=1
									start=max([inexon[0],refexon[0]])
									end=min([inexon[1],refexon[1]])
									nmatch+=(end-start+1)
									ematch=(end-start+1)
									if (exonid not in ematch_dict) or (ematch_dict[exonid][0]<ematch):
										ematch_dict[exonid]=[ematch,refexon[2],0]
										if (refexon[0]-1)<=inexon[0]<=(refexon[0]+1):
											ematch_dict[exonid][2]+=1
										if (refexon[1]-1)<=inexon[1]<=(refexon[1]+1):
											ematch_dict[exonid][2]+=1
							if overlexon==1:
								oec+=1
						sm=[]
						pos=0
						for inexon in exon_list:
							pos+=1
							syes,eyes=0,0
							for refexon in ref_exon_list:
								if (refexon[0]-1)<=inexon[0]<=(refexon[0]+1):
									splicem+=1
									syes=1
								if (refexon[1]-1)<=inexon[1]<=(refexon[1]+1):
									splicem+=1
									eyes=1
							sm.extend([syes,eyes])
						if len(exon_list)!=len(ref_exon_list):
							decide='different'
						elif len(exon_list)==1:
							comb_exon=[[exon_list[0][0],'e'],[exon_list[0][1],'e'],[ref_exon_list[0][0],'r'],[ref_exon_list[0][1],'r']]
							len_list=[abs(comb_exon[1][0]-comb_exon[0][0]),abs(comb_exon[3][0]-comb_exon[2][0])]
							len_list.sort()
							comb_exon.sort()
							if comb_exon[1][1]==comb_exon[2][1]:
								decide='same'
							elif float(len_list[0])/float(len_list[1])>=0.75:
								decide='same'
							else:
								decide='different'
						elif 0 in sm[1:-1]:
							decide='different'
						else:
							decide='same'
						if decide=='same':
							tmatch+=1
						new.write(last_chro+'	'+str(one[0][0])+'	'+str(one[0][1])+'	'+str(trans[0])+'	'+str(trans[1])+'	'+strand+'	'+name+'	'+one[3]+'	'+str(float(oec)*100/len(exon_list))+'	'+str(float(nmatch)*100/one[4])+'	'+str(float(splicem)*100/one[5])+'	'+decide+'\n')
						written+=1
						max_exon_splice.append([float(oec)*100/len(exon_list),(float(nmatch)*100/one[4]),(float(splicem)*100/one[5]),one[3],sm,exon_list,ref_exon_list])
				if len(max_exon_splice)!=0:
					max_exon_splice.sort(reverse=True)
					ecount_dict[name]=[max_exon_splice[0][0],len(exon_list)]
					if len(max_exon_splice[0][-1])!=len(max_exon_splice[0][-2]):
						decision='different'
					elif len(max_exon_splice[0][-1])==1:
						comb_exon=[[max_exon_splice[0][-2][0][0],'e'],[max_exon_splice[0][-2][0][1],'e'],[max_exon_splice[0][-1][0][0],'r'],[max_exon_splice[0][-1][0][1],'r']]
						comb_exon.sort()
						len_list=[abs(max_exon_splice[0][-2][0][1]-max_exon_splice[0][-2][0][0]),abs(max_exon_splice[0][-1][0][1]-max_exon_splice[0][-1][0][0])]
						len_list.sort()
						if comb_exon[1][1]==comb_exon[2][1]:
							decision='same'
						elif max_exon_splice[0][1]>=75:
							decision='same'
						else:
							decision='different'
					elif 0 in max_exon_splice[0][-3][1:-1]:
						decision='different'
					else:
						if max_exon_splice[0][1]>=75:
							decision='same'
						else:
							decision='different'
					if decision=='same':
						bmatch+=1
					newm.write(name+'	'+max_exon_splice[0][3]+'	'+str(max_exon_splice[0][0])+'	'+str(max_exon_splice[0][1])+'	'+str(max_exon_splice[0][2])+'	'+decision+'\n')
					tmatch_dict[name]=(max_exon_splice[0][3],str(max_exon_splice[0][1]),str(max_exon_splice[0][2]),decision)
					match+=1
				else:
					ecount_dict[name]=[0,len(exon_list)]
				if found==1:
					unique+=1
			tcount+=1
			trans=[int(split[3]),int(split[4])]
			trans.sort()
			strand=split[6]
			for i in range(8,len(split),1):
				if split[i]=='transcript_id':
					name=split[i+1][1:-2]
					break
			exon_list=[]
		elif split[2]=='exon':
			ecount+=1
			coord=[int(split[3]),int(split[4])]
			coord.sort()
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
				pos=0
				oec=0
				for inexon in exon_list:
					pos+=1
					exonid=name+'e'+str(pos)
					overlexon=0
					for refexon in ref_exon_list:
						if inexon[1]<refexon[0]:
							break
						elif inexon[0]>refexon[1]:
							continue
						else:
							overlexon=1
							start=max([inexon[0],refexon[0]])
							end=min([inexon[1],refexon[1]])
							nmatch+=(end-start+1)
							ematch=(end-start+1)
							if (exonid not in ematch_dict) or (ematch_dict[exonid][0]<ematch):
								ematch_dict[exonid]=[ematch,refexon[2],0]
								if (refexon[0]-1)<=inexon[0]<=(refexon[0]+1):
									ematch_dict[exonid][2]+=1
								if (refexon[1]-1)<=inexon[1]<=(refexon[1]+1):
									ematch_dict[exonid][2]+=1
					if overlexon==1:
						oec+=1
				sm=[]
				for inexon in exon_list:
					syes,eyes=0,0
					for refexon in ref_exon_list:
						if (refexon[0]-1)<=inexon[0]<=(refexon[0]+1):
							splicem+=1
							syes=1
						if (refexon[1]-1)<=inexon[1]<=(refexon[1]+1):
							splicem+=1
							eyes=1
					sm.extend([syes,eyes])
				if len(exon_list)!=len(ref_exon_list):
					decide='different'
				elif len(exon_list)==1:
					comb_exon=[[exon_list[0][0],'e'],[exon_list[0][1],'e'],[ref_exon_list[0][0],'r'],[ref_exon_list[0][1],'r']]
					len_list=[abs(comb_exon[1][0]-comb_exon[0][0]),abs(comb_exon[3][0]-comb_exon[2][0])]
					len_list.sort()
					comb_exon.sort()
					if comb_exon[1][1]==comb_exon[2][1]:
						decide='same'
					elif max_exon_splice[0][1]>=75:
						decide='same'
					else:
						decide='different'
				elif 0 in sm[1:-1]:
					decide='different'
				else:
					if max_exon_splice[0][1]>=75:
						decide='same'
					else:
						decide='different'
				if decide=='same':
					tmatch+=1
				new.write(last_chro+'	'+str(one[0][0])+'	'+str(one[0][1])+'	'+str(trans[0])+'	'+str(trans[1])+'	'+strand+'	'+name+'	'+one[3]+'	'+str(float(oec)*100/len(exon_list))+'	'+str(float(nmatch)*100/one[4])+'	'+str(float(splicem)*100/one[5])+'	'+decide+'\n')
				written+=1
				max_exon_splice.append([float(oec)*100/len(exon_list),(float(nmatch)*100/one[4]),(float(splicem)*100/one[5]),one[3],sm,exon_list,ref_exon_list])
		if len(max_exon_splice)!=0:
			max_exon_splice.sort(reverse=True)
			ecount_dict[name]=[max_exon_splice[0][0],len(exon_list)]
			if len(max_exon_splice[0][-1])!=len(max_exon_splice[0][-2]):
				decision='different'
			elif len(max_exon_splice[0][-1])==1:
				comb_exon=[[max_exon_splice[0][-2][0][0],'e'],[max_exon_splice[0][-2][0][1],'e'],[max_exon_splice[0][-1][0][0],'r'],[max_exon_splice[0][-1][0][1],'r']]
				comb_exon.sort()
				len_list=[abs(max_exon_splice[0][-2][0][1]-max_exon_splice[0][-2][0][0]),abs(max_exon_splice[0][-1][0][1]-max_exon_splice[0][-1][0][0])]
				len_list.sort()
				if comb_exon[1][1]==comb_exon[2][1]:
					decision='same'
				elif float(len_list[0])/float(len_list[1])>=0.75:
					decision='same'
				else:
					decision='different'
			elif 0 in max_exon_splice[0][-3][1:-1]:
				decision='different'
			else:
				decision='same'
			if decision=='same':
				bmatch+=1
			newm.write(name+'	'+max_exon_splice[0][3]+'	'+str(max_exon_splice[0][0])+'	'+str(max_exon_splice[0][1])+'	'+str(max_exon_splice[0][2])+'	'+decision+'\n')
			tmatch_dict[name]=(max_exon_splice[0][3],str(max_exon_splice[0][1]),str(max_exon_splice[0][2]),decision)
			match+=1
		else:
			ecount_dict[name]=[0,len(exon_list)]
		if found==1:
			unique+=1
	yui.close()
	new.close()
	newm.close()
	print 'Number of trancripts:',tcount
	print 'Number of exons:',ecount
	print 'Number of pairs written:',written
	print 'Number of matched trancripts:',unique
	print 'Number of matched overlap:', match
	print 'Number of unmatched trancripts:',tcount-unique
	print 'Number of total same transcripts:',tmatch
	print 'Number of best same transcripts:',bmatch
	print 'Number of collected exon data:',len(ematch_dict)
	for one in ematch_dict:
		print 'Example exon:',one, ematch_dict[one]
		break
	yui=open(ingtf,'r')
	new=open(outgtf,'a')
	tt,te,mt,me=0,0,0,0
	tma,tno,tva=0,0,0
	for line in yui:
		split=line.split()
		if len(split[0])>2 or split[0]=='MT':
			continue
		if split[2]=='transcript':
			tt+=1
			for i in range(8,len(split),1):
				if split[i]=='transcript_id':
					tid=split[i+1][1:-2]
					break
			try:
				details=tid_dict[tmatch_dict[tid][0]]
				mt+=1
				if ecount_dict[tid][0]==0:
					cla='novel'
					tno+=1
				elif ecount_dict[tid][0]==100 and tmatch_dict[tid][3]=='same':
					cla='matching'
					tma+=1
				else:
					cla='variant'
					tva+=1
				new.write(line[:-1]+' transcript_class "'+cla+'";'+' exon_count "'+str(ecount_dict[tid][1])+'";'+' overlapping_exon_count "'+str(ecount_dict[tid][0])+'";'+' percent_exon_match "'+str(tmatch_dict[tid][1])+'";'+' percent_splice_match "'+str(tmatch_dict[tid][2])+'";'+' decision "'+str(tmatch_dict[tid][3])+'";'+' gene_name "'+details[0]+'";'+' transcript_name "'+details[1]+'";'+' GENCODE_gene_id "'+details[2]+'";'+' GENCODE_transcript_id "'+tmatch_dict[tid][0]+'";'+' GENCODE_transcript_type "'+details[3]+'";\n')
			except:
				cla='novel'
				tno+=1
				new.write(line[:-1]+' transcript_class "'+cla+'";'+' exon_count "'+str(ecount_dict[tid][1])+'";'+' overlapping_exon_count "'+str(ecount_dict[tid][0])+'";'+' percent_exon_match "0";'+' percent_splice_match "0";'+' decision "different";'+'";\n')
		elif split[2]=='exon':
			te+=1
			for i in range(8,len(split),1):
				if split[i]=='transcript_id':
					tid=split[i+1][1:-2]
				if split[i]=='exon_number':
					ename=tid+'e'+split[i+1][1:-2]
			try:
				details=tid_dict[tmatch_dict[tid][0]]
				try:
					perc_ematch=float(ematch_dict[ename][0])*100/(abs(int(split[3])-int(split[4]))+1)
					me+=1		
					new.write(line[:-1]+' transcript_class "'+cla+'";'+' transcript_percent_exon_match "'+str(tmatch_dict[tid][1])+'";'+' transcript_percent_splice_match "'+str(tmatch_dict[tid][2])+'";'+' decision "'+str(tmatch_dict[tid][3])+'";'+' gene_name "'+details[0]+'";'+' transcript_name "'+details[1]+'";'+' GENCODE_gene_id "'+details[2]+'";'+' GENCODE_transcript_id "'+tmatch_dict[tid][0]+'";'+' GENCODE_exon_id "'+ematch_dict[ename][1]+'";'+' GENCODE_exon_overlap_perc "'+str(perc_ematch)+'";'+' GENCODE_splice_overlap_perc "'+str(float(ematch_dict[ename][2])*100/2)+'";'+' GENCODE_transcript_type "'+details[3]+'";\n')
				except:
					new.write(line[:-1]+' transcript_class "'+cla+'";'+' transcript_percent_exon_match "'+str(tmatch_dict[tid][1])+'";'+' transcript_percent_splice_match "'+str(tmatch_dict[tid][2])+'";'+' decision "'+str(tmatch_dict[tid][3])+'";'+' gene_name "'+details[0]+'";'+' transcript_name "'+details[1]+'";'+' GENCODE_gene_id "'+details[2]+'";'+' GENCODE_transcript_id "'+tmatch_dict[tid][0]+'";'+' GENCODE_exon_overlap_perc "0";'+' GENCODE_splice_overlap_perc "0";'+' GENCODE_transcript_type "'+details[3]+'";\n')
			except:
				new.write(line[:-1]+' transcript_class "'+cla+'";'+' transcript_percent_exon_match "0";'+' transcript_percent_splice_match "0";'+' decision "different";'+' GENCODE_exon_overlap_perc "0";'+' GENCODE_splice_overlap_perc "0";\n')
	yui.close()
	new.close()
	print 'Number of transcript:',tt
	print '	Matching transcript:',mt
	print 'Number of exons:',te
	print '	Matching exons:',me
	print 'Example ename:',ename
	print 'Number of matching transcripts:',tma
	print 'Number of variant transcripts:',tva
	print 'Number of novel transcripts:',tno



def count_threshold(infile,threshold):
	yui=open(infile,'r')
	exon_dict=dict()
	splice_dict=dict()
	count_dict=dict()
	all_dict=dict()
	total=0
	for line in yui:
		split=line.split()
		try:
			ecount=float(split[2])
			exon=float(split[3])
			splice=float(split[4])
		except:
			continue
		total+=1
		if exon>=threshold:
			exon_dict[split[1]]=0
		if ecount>=threshold:
			count_dict[split[1]]=0
		if splice>=threshold:
			splice_dict[split[1]]=0
			if exon>=threshold and ecount>=threshold:
				all_dict[split[1]]=0
	yui.close()
	print 'Total:',total
	print 'Exon count threshold pass:',len(exon_dict),float(len(count_dict))/total
	print 'Exon coverage threshold pass:',len(exon_dict),float(len(exon_dict))/total
	print 'Splice threshold pass:',len(splice_dict),float(len(splice_dict))/total
	print 'Exon count, exon coverage and splice threshold pass:',len(all_dict),float(len(all_dict))/total



from optparse import OptionParser
import os,sys
parser = OptionParser()
parser.add_option('-i','--ingtf', dest='ingtf',help='Absolute path to the input gtf file [required]',type='str',default='')
parser.add_option('-r','--refgtf', dest='refgtf',help='Absolute path to the reference, this case GENCODE, gtf file [required]',type='str',default='')
parser.add_option('-g','--outgtf', dest='outgtf',help='Absolute path to the output updated gtf file [required]',type='str',default='')
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

compare_GTF(options.ingtf,options.refgtf,options.outgtf,options.outfile,options.match_file)

print 'Getting the summaries...'

print 'Threshold considered:',options.threshold

count_threshold(options.match_file,options.threshold)
