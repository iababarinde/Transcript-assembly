def update(sr_gtf,lr_gtf,match,outgtf):
	yui=open(match,'r')
	sr_dict=dict()
	overlap=dict()
	lr_dict=dict()
	sr_value=dict()
	lr_value=dict()
	for line in yui:
		split=line.split()
		if split[-1]=='same':
			if (split[1] not in lr_value) or (float(split[-2])+float(split[-3])>lr_value[split[1]][1]):
				lr_value[split[1]]=[split[0],float(split[-2])+float(split[-3])]
			if (split[0] not in sr_value) or (float(split[-2])+float(split[-3])>sr_value[split[0]][1]):
				sr_value[split[0]]=[split[1],float(split[-2])+float(split[-3])]
		else:
			overlap[split[0]]=split[1].split('.')[1]
	yui.close()
	norep=0
	for trans in lr_dict:
		if (lr_dict[trans][0] in sr_dict) and sr_dict[lr_dict[trans][0]][0]==trans:
			norep+=1
			lr_dict[trans]=lr_dict[trans][0]
			sr_dict[lr_dict[trans][0]]=trans
		else:
			overlap[lr_dict[trans][0]]=trans.split('.')[1]
	print len(sr_dict),'matched short read transcripts'
	print len(lr_dict),'matched long read transcripts'
	print len(overlap),'overlapping different transcript'
	print norep,'unique matched transcripts'
	chro_dict=dict()
	yui=open(lr_gtf,'r')
	count=0
	total=0
	chro_list=[]
	lr_exon=dict()
	lr_trans=dict()
	for line in yui:
		split=line.split()
		if split[2]=='transcript':
			total+=1
		for i in range(7,len(split),1):
			if split[i]=='transcript_id':
				tid=split[i+1][1:-2]
				if tid in lr_dict:
					if split[2]=='transcript':
						coord=[int(split[3]),int(split[4])]
						coord.sort()
						lr_trans[tid]=coord
						lr_exon[tid]=[]
					elif split[2]=='exon':
						coord=[int(split[3]),int(split[4])]
						coord.sort()
						lr_exon[tid].append(coord)
					break
				else:
					if split[0] not in chro_dict:
						chro_dict[split[0]]=[]
						chro_list.append(split[0])
					if tid in overlap:
						for i in range(7,len(split),1):
							if split[i]=='gene_id':
								gid=split[i+1]
								break
						line=line.replace(gid,'"HPSCSR'+overlap[tid]+'";')
					if split[2]=='transcript':
						count+=1
						chro_dict[split[0]].append([int(split[3]),int(split[4]),line])
					elif split[2]=='exon':
						chro_dict[split[0]][-1].append(line)
				break
	yui.close()
	print 'Total Number of LR:',total
	print 'Number in the dictionary:',len(chro_dict)
	print count,'additional transcripts to be written'
	yui=open(sr_gtf,'r')
	count=0
	for line in yui:
		split=line.split()
		if split[0] not in chro_dict:
			chro_dict[split[0]]=[]
			chro_list.append(split[0])
		if split[2]=='transcript':
			count+=1
			chro_dict[split[0]].append([int(split[3]),int(split[4]),line])
		elif split[2]=='exon':
			chro_dict[split[0]][-1].append(line)
	yui.close()
	print count,'SR transcripts found'
	new=open(outgtf,'a')
	written=0
	sr=0
	lr=0
	srlr=0
	for chromo in chro_list:
		current=chro_dict[chromo]
		current.sort()
		for transcript in current:
			written+=1
			for line in transcript[2:]:
				split=line.split()
				newl=''
				for i in range(8,len(split),1):
					if split[i]=='transcript_id':
						tid=split[i+1][1:-2]
						break
				if split[2]=='transcript':
					if tid in sr_dict:
						tcoord=[int(split[3]),int(split[4])]
						tcoord.sort()
						try:
							lr_id=sr_dict[tid]
							trans=lr_trans[lr_id]
							if tcoord[0]<tcoord[1]<trans[0] or trans[0]<trans[1]<tcoord[0]:
								pass
							else:
								overlen,tstart,tend=(min([tcoord[1],trans[1]])-max([tcoord[0],trans[0]])+1),min([tcoord[0],trans[0]]),max([tcoord[1],trans[1]])
								use_overlen,use_start,use_end=overlen,tstart,tend
								split[3],split[4]=str(use_start),str(use_end)
						except:
							pass
					for one in split[:8]:
						newl+=one+'	'
					for one in split[8:-6]:
						newl+=one+' '
				elif split[2]=='exon':
					ecoord=[int(split[3]),int(split[4])]
					ecoord.sort()
					use_overlen,use_start,use_end=0,0,0
					try:
						lr_id=sr_dict[tid]
						for exon in lr_exon[lr_id]:
							if ecoord[0]>exon[1]:
								continue
							elif ecoord[1]<exon[0]:
								break
							else:
								overlen,estart,eend=(min([ecoord[1],exon[1]])-max([ecoord[0],exon[0]])+1),min([ecoord[0],exon[0]]),max([ecoord[1],exon[1]])
								use_overlen,use_start,use_end=overlen,estart,eend
						split[3],split[4]=str(use_start),str(use_end)
					except:
						pass
					for one in split[:8]:
						newl+=one+'	'
					for one in split[8:-2]:
						newl+=one+' '
				if tid in sr_dict:
					newl+='evidence "SR+LR"; LR_ID "'+sr_dict[tid]+'";'
					if split[2]=='transcript':
						srlr+=1
				elif 'HPSCSR' in tid:
					newl+='evidence "SR";'
					if split[2]=='transcript':
						sr+=1
				elif 'HPSCLR' in tid:
					newl+='evidence "LR";'
					if split[2]=='transcript':
						lr+=1
				new.write(newl+'\n')
	new.close()
	print written,'transcripts written'
	print 'SR only:',sr
	print 'LR only:',lr
	print 'Both SR and LR:',srlr
	print tid,'last transcript'


from optparse import OptionParser
import os,sys
parser = OptionParser()
parser.add_option('-s','--sr_gtf', dest='sr_gtf',help='Absolute path to the short read gtf file [required]',type='str',default='')
parser.add_option('-l','--lr_gtf', dest='lr_gtf',help='Absolute path to the long read gtf file [required]',type='str',default='')
parser.add_option('-m','--match_file', dest='match_file',help='Absolute path to the match file obtained from comparegtf',type='str',default='')
parser.add_option('-o','--outgtf', dest='outgtf',help='Absolute path to the output updated gtf file',type='str',default='')
(options,args)=parser.parse_args()

if (options.sr_gtf=='') or (options.lr_gtf=='') or (options.match_file=='') or (options.outgtf==''):
	print '\nRequired filed(s) not supplied\n# The program builds an updated gtf file from short and long read gtf files\n'
	parser.print_help()
	sys.exit(1)

update(options.sr_gtf,options.lr_gtf,options.match_file,options.outgtf)