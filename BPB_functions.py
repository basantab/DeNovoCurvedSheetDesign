"""
This library contains the utility functions for generating the curved sheet folds
described in "Principles for designing proteins with cavities formed by curved beta-sheets"
Enrique Marcos - emarcos82@gmail.com
Benjamin Basanta - basantab@uw.edu
"""

import itertools
import re
import sys
import os
import copy
from Blueprint import Blueprint
import numpy as np

def Bulged(strand):
        """
        Return absolute position of bulge in input strand if present. Otherwise, return None.
        """
        flag=False ; bulgepos=None
        for i in range(1,len(strand.bp_data)-1):
                prev_res = strand.bp_data[i-1]
                res = strand.bp_data[i]
                next_res = strand.bp_data[i+1]
                if 'EA' == res[2] and 'EB' == prev_res[2] and 'EB' == next_res[2]:
                    flag=True
                    bulgepos = res[0]
        return bulgepos

def MakePlainBlueprint(ss,comb,bluefile):
        """
        Generates the first blueprint, which serves as a guide for placing bulges
        in the following steps
        """
        c = comb
        out_file = open(bluefile,'w')
        struct = {'H':'HA', 'E':'EB', 'L':'LA'}
        total_length=sum(c)
        k=0
        curr_ss = ss[k]
        for i in range(1,total_length+1):
                if i>sum(c[:k+1]):
                        k+=1
                        curr_ss = ss[k]
                out_file.write('0  V  %s  R\n' %(curr_ss))
        out_file.close()

def MakeRefBlueprint(ss,comb,bulges,**kwargs):
        """
        Generates reference blueprint at the beginning of the process
        """
        refbluefile = kwargs.get('refblue')
        c = comb
        out_file = open(refbluefile,'w')
        struct = {'H':'HA', 'E':'EB', 'L':'LA'}
        total_length=sum(c)
        k=0
        curr_ss = ss[k]
        for i in range(1,total_length+1):
                if i>sum(c[:k+1]):
                        k+=1
                        curr_ss = ss[k]
                out_file.write('0  V  %s  R\n' %(curr_ss))
        out_file.close()

        # Put bulges
        blue = Blueprint(refbluefile)
       	bluelist = []
       	for seg in blue.segments:
       	        if seg.id in bulges.keys():
       	                for j, res in enumerate(seg.bp_data):
       	                        if j == bulges[seg.id]-1:
       	                                res[2] = 'EA'
       	blue.dump_blueprint(refbluefile)

def Shift(**kwargs):
        """
        Calculates register shifts between strands based on their pairing and length
        """
        refbluefile = kwargs.get('refblue')
        seg1 = kwargs.get('seg1')
        seg2 = kwargs.get('seg2')
        refblue = Blueprint(refbluefile)
        nseg1 = len( refblue.segment_dict[seg1].bp_data )
        nseg2 = len( refblue.segment_dict[seg2].bp_data )
        shift = nseg1 - nseg2
        return shift

def MakeFirstBlueprint(**kwargs):
        """
        Generates the first blueprint in the protein assembly process. This blueprint
        contains only the two central strands in the sheet
        """
        tail = [[0, 'V', 'L', 'R']]
        refbluefile = kwargs.get('refblue')
        segments = kwargs.get('segments')
        specific_abego = kwargs.get('specific_abego')
        newbluefile = kwargs.get('newblue')
        ss_pairing = kwargs.get('ss_pairing')
        hs_pairing = kwargs.get('hs_pairing')
        hh_pairing = kwargs.get('hh_pairing')
        seg_adapt = kwargs.get('adapt')

        refblue = Blueprint(refbluefile)

        seg1 = seg_adapt.keys()[0] # strand to adapt
        seg2 = seg_adapt[seg1] # referemce strand
        nseg1 = len( refblue.segment_dict[seg1].bp_data )
        nseg2 = len( refblue.segment_dict[seg2].bp_data )
        shift = nseg1 - nseg2

        bp_data_new = []
        for seg in segments:
                if shift == 0:
                        for res in refblue.segment_dict[seg].bp_data:
                                bp_data_new.append(res)
                elif shift > 0:
                        if seg==seg1:
                                for k in range(0,nseg1-shift):
                                        bp_data_new.append(refblue.segment_dict[seg1].bp_data[k])
                        else:
                                for res in refblue.segment_dict[seg].bp_data:
                                        bp_data_new.append(res)
                elif shift < 0:
                        if seg==seg2:
                                for k in range(shift,nseg2):
                                        bp_data_new.append(refblue.segment_dict[seg2].bp_data[k])
                                        bp_data_new.append(refblue.segment_dict[seg1].bp_data[k])
                        else:
                                for res in refblue.segment_dict[seg].bp_data:
                                        bp_data_new.append(res)

	refblue.bp_data = tail + bp_data_new + tail

	# write strand pairing
	ss_line = "SSPAIR "
	# Get first strand to set "start"
	segments.sort()
	npairs=0
	for i,key in enumerate(ss_pairing.keys()):
            for st in ss_pairing[key]:
		s1=key
		s2 = st[:st.find('.')]
		orient = st[-1]
		# compute effective length of strands for calculating register shift
		#------
		if Bulged(refblue.segment_dict[s1]):
			n1 = len( refblue.segment_dict[s1].bp_data ) - 1
		else:
			n1 = len( refblue.segment_dict[s1].bp_data )
                if Bulged(refblue.segment_dict[s2]):
                        n2 = len( refblue.segment_dict[s2].bp_data ) - 1
                else:
                        n2 = len( refblue.segment_dict[s2].bp_data )
		#------
		if orient == 'A':
			shift = n1-n2

		if npairs==0:
			ss_line += '%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)
		else:
			ss_line += ';%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)
		npairs+=1

	header = [ss_line]
	if hs_pairing != None:
		header.append( hs_pairing )

	if hh_pairing !=None:
		header.append( hh_pairing )

	refblue.dump_blueprint(newbluefile,header_lines=header)
	blue = Blueprint('%s' %(newbluefile))
	if specific_abego: # only abegos for specific positions
 		for seg in specific_abego:
			for position in specific_abego[seg]:
				pos,letter =  position
				blue.segment_dict[seg].bp_data[pos][2]+=letter
				if letter == 'E' or letter == 'G':
					blue.segment_dict[seg].bp_data[pos][1] = 'G'

	blue.dump_blueprint(newbluefile,header_lines=header)
	# Abego 'B' for building strand
	os.system("sed  's/ E / EB/g;s/ H / HA/g' %s > %s.b" %(newbluefile,newbluefile))

def AddSegmentToBlueprint(**kwargs):
        """
        Generates subsequent blueprints, adding secondary structure elements
        """
	tail = [[0, 'V', 'L', 'R']]
	refbluefile = kwargs.get('refblue')
	segments = kwargs.get('segments')
	blue0file= kwargs.get('blue0')
	newbluefile = kwargs.get('newblue')
	append = kwargs.get('append')
	insert_between_first = kwargs.get('insert_between_first')
	insert_between_last = kwargs.get('insert_between_last')
	ss_pairing = kwargs.get('ss_pairing')
	ss_pairing_shift = kwargs.get('ss_pairing_shift')
	hs_pairing = kwargs.get('hs_pairing')
	hh_pairing = kwargs.get('hh_pairing')
	seg_abego = kwargs.get('seg_abego')
	specific_abego = kwargs.get('specific_abego')
	insert = kwargs.get('insert')
	only_remodel = kwargs.get('only_remodel',False)

	blue0 = Blueprint(blue0file)
	blue0.reindex_blueprint(start=1)
	blue0.freeze_all()

	refblue = Blueprint(refbluefile)
	bp_data_new = []
	for seg in segments:
		for res in refblue.segment_dict[seg].bp_data:
			bp_data_new.append(res)

        if append==True:
                blue0.bp_data[-2][3]='R'
                blue0.bp_data[-1][3]='R' ; blue0.bp_data[-1][2]='L'
                bp_data_new = blue0.bp_data + bp_data_new[1:] + tail

        elif append==False:
                blue0.bp_data[0][3]='R' #; blue0.bp_data[0][2]='L'
                blue0.bp_data[1][3]='R'
		if segments[-1][0] != 'L':
			blue0.bp_data[0][2]= bp_data_new[-1][2]
                bp_data_new =  tail + bp_data_new[:-1] + blue0.bp_data

	elif insert_between_first:
		# replace original residues by insert ones. so we can use original bp as blue0
		index1 = insert_between_first - 1
		index2 = insert_between_last - 1
		if index1>1:  # we need at least one point fixed (not R)
			blue0.bp_data[index1][3]='R' #; blue0.bp_data[index1][2]='LX'

		bp_data_new = blue0.bp_data[:index1+1] + bp_data_new # + blue0.bp_data[index2:]
		shift = index2-index1-1
		for k in range(index2,len(blue0.bp_data)):
			blue0.bp_data[k][0]-=shift
		blue0.bp_data[index2][3]='R'
		bp_data_new = bp_data_new + blue0.bp_data[index2:]

        else:
                bp_data_new = blue0.bp_data

	blue0_top = blue0.topology() # for abego conversion later
	blue0_cp = copy.deepcopy(blue0)

        # if we dont add residues then we dont need to replace bp_data. it is already update with R
        if only_remodel==False or isinstance(only_remodel,list):
                blue0.bp_data = bp_data_new

        # write strand pairing
	if ss_pairing != None:
         ss_line = "SSPAIR "
  	 keys = ss_pairing.keys()
	 keys.sort()
	 npairs=0
         for i,key in enumerate(keys):
	    for st in ss_pairing[key]:
		s1=key
	        s2 = st[:st.find('.')]
                orient = st[-1]

                # compute effective length of strands for calculating register shift
                #------
                if Bulged(refblue.segment_dict[s1]):
                        n1 = len( refblue.segment_dict[s1].bp_data ) - 1
                else:
                        n1 = len( refblue.segment_dict[s1].bp_data )
                if Bulged(refblue.segment_dict[s2]):
                        n2 = len( refblue.segment_dict[s2].bp_data ) - 1
                else:
                        n2 = len( refblue.segment_dict[s2].bp_data )
                #------

		shift=99
		if npairs==0:
                        ss_line += '%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)
                else:
                        ss_line += ';%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)
		npairs+=1

	header=[]
	if ss_pairing_shift != None:
	        header.append( ss_pairing_shift )
	else:
		header.append( ss_line )
        if hs_pairing != None:
		#for i in hs_pairing:
		#	header.append( i )
                header.append( hs_pairing )

        if hh_pairing !=None:
                header.append( hh_pairing )

	blue0.dump_blueprint(newbluefile, header_lines=header)


        # insertion
        if insert:
                seg = insert.keys()[0]
                shift = insert[seg]
                blue = Blueprint('%s' %(newbluefile))
                blue_aux = copy.deepcopy(blue)
                blue_aux.reindex_blueprint()
                insert_index = blue_aux.segment_dict[seg].bp_data[-1][0] # we want to insert just before the 1st residue of the next seg, so that we insert at the very end

                for k in range(shift):
                        blue.bp_data.insert(insert_index, [0, 'V', 'E', 'R'])

                blue.dump_blueprint(newbluefile,header_lines=header)

	# abego loop
	# Adapt global abego motif to current stage
	#---------------------------------
	if seg_abego != None:
		new_abego={};conversor={}
		r=re.compile('[HEL]')
		top= blue0_cp.topology()
		curr_ss = r.findall(top)
		counter=0
		loop_counter=0
		for seg in segments:
			ss = seg[0]
			if append:
				if ss=='L':
					newindex = curr_ss.count(ss)+loop_counter
					loop_counter+=1
				else:
					newindex = curr_ss.count(ss) + 1
					curr_ss.append(ss)
			elif append==False:
				curr_ss.insert(counter,ss)
				newindex = curr_ss[:counter+1].count(ss)
				if ss=='L': newindex+=1 # because Nter is loop L1 (and is not added through segments)
				counter+=1
			elif insert_between_first:
				r=re.compile('[HEL]\d') ; seg_list = r.findall(top)
				flag=False
				for k,sg in enumerate(seg_list): # find segment whhere insert_first is located
					for res in blue0_cp.segment_dict[sg].bp_data:
						if res[0]==insert_between_first+1:
							flag=True
							break
					if flag:
						break

				first_seg=k
				newindex = curr_ss[:first_seg+1+counter].count(ss)
				counter+=1

			conversor[seg]='%s%s' %(ss,newindex)

		for seg in seg_abego:
			new_abego[conversor[seg]] = seg_abego[seg]

		print conversor
        	blue = Blueprint('%s' %(newbluefile))
        	if new_abego != None:
        	        for seg in new_abego.keys():
        	        	abego=new_abego[seg]
        	        	for i,res in enumerate(blue.segment_dict[seg].bp_data):
        	                	res[2]+=abego[i]
					# Added on 20160612:
					if abego[i]=='G' or abego[i]=='E':
						res[1]='G'

        	blue.dump_blueprint(newbluefile,header_lines=header)

	if specific_abego: # only abegos for specific positions
		blue = Blueprint('%s' %(newbluefile))
		for seg in specific_abego:
			pos,letter =  specific_abego[seg]
			blue.segment_dict[seg].bp_data[pos][2]+=letter
			blue.segment_dict[seg].bp_data[pos][1] = 'A'

		blue.dump_blueprint(newbluefile,header_lines=header)

	# only remodel.
	blue = Blueprint('%s' %(newbluefile))
        if isinstance(only_remodel,list):
                for seg in only_remodel:
                        blue.remodel_segment(id=seg)

		blue.dump_blueprint(newbluefile,header_lines=header)

	# Abego 'B' for building strand
	os.system("sed  's/ E / EB/g;s/ H / HA/g' %s > %s.b" %(newbluefile,newbluefile))

def write_dummy_pdb(filename):
    """
    Creates initial dummy PDB file that is used by the blueprint build mover
    """
    dummy_pdb = 'ATOM      1  N   GLY A   1       0.346   1.210   0.744  1.00  0.00 \n\
ATOM      2  CA  GLY A   1       1.687   1.135   0.174  1.00  0.00 \n\
ATOM      3  C   GLY A   1       2.383   2.488   0.222  1.00  0.00 \n\
ATOM      4  O   GLY A   1       2.996   2.918  -0.752  1.00  0.00 \n\
ATOM      5 1H   GLY A   1      -0.448   0.981   0.185  1.00  0.00 \n\
ATOM      6 2H   GLY A   1       0.109   0.647   1.536  1.00  0.00 \n\
ATOM      7 3H   GLY A   1       0.005   2.079   1.101  1.00  0.00 \n\
ATOM      8 1HA  GLY A   1       2.277   0.413   0.741  1.00  0.00 \n\
ATOM      9 2HA  GLY A   1       1.615   0.810  -0.864  1.00  0.00 \n\
ATOM     10  N   GLY A  2       2.284   3.158   1.368  1.00  0.00 \n\
ATOM     12  CA  GLY A  2       3.918   4.450   2.676  1.00  0.00 \n\
ATOM     13  C   GLY A  2       3.551   4.379   3.850  1.00  0.00 \n\
ATOM     14  O   GLY A  2       1.859   5.564   1.752  1.00  0.00 \n\
ATOM     15 1H   GLY A  2       1.061   5.930   0.512  1.00  0.00 \n\
ATOM     16 2H   GLY A  2      -0.012   6.930   0.747  1.00  0.00 \n\
ATOM     17 3H   GLY A  2      -0.863   7.182  -0.405  1.00  0.00 \n\
ATOM     18 1HA  GLY A  2      -0.564   8.037  -1.404  1.00  0.00 \n\
ATOM     19 2HA  GLY A  2       0.540   8.749  -1.379  1.00  0.00 \n'
    out = open(filename, 'w')
    out.write(dummy_pdb)

def Bulges(blue):
        """
        Returns all bulges in strainds
        """
        bulge_list=[]
        for i,res in enumerate(blue.bp_data):
                if res[2]=='EA':
                        bulge_list.append(res[0])
        return bulge_list

def XMLReplaceXXXYYY(**kwargs):
        """
        Replaces strings in XML lines.
        """
        xml_lines = kwargs.get('xml_lines')
        identifier = kwargs.get('identifier')
        xxx = kwargs.get('xxx')
        yyy = kwargs.get('yyy')
	zzz = kwargs.get('zzz')

        for i,line in enumerate(xml_lines):
                if identifier in line:
                        if xxx:
                                line = line.replace('xxx','%s' %(xxx))
				xml_lines[i] = line

                        if yyy:
                                line = line.replace('yyy','%s' %(yyy))
				xml_lines[i] = line

			if zzz:
				line = line.replace('zzz','%s' %(zzz))


def AmbiguousConstraints(list1,list2):
        """
        Returns the string declaring ambigous constraints between
        all atoms N of residues in list1 and O atoms of residues in list2
        """
	st='AmbiguousConstraint\n'
	for res1 in list1:
		for res2 in list2:
			st += "AtomPair N %i O %i BOUNDED 3.5 4.5 0.5\n" %(res1,res2)
	st+='END_AMBIGUOUS\n'
	return st

def MotifTopology(topol,motif_filename):
        """
        Deprecated function
        """
	elements = re.compile("[HEL]")
	ss = elements.findall(topol)
	relengths = re.compile("(\d+)-(\d+)")
	relengths2 = re.compile("(\d+),(\d+)") # for specific lengths, not ranges
	lengths = relengths.findall(topol)
	lengths2 = relengths2.findall(topol)

	## for interpreting specific lengths
	index=-1
	specific=[]
	for st in topol:
		if st=='[':
			index+=1
		if st==',':
			specific.append(index)

	comb=[] ; j=0 ; k=0
	for i in range(len(ss)):
		if i in specific:
			frag = lengths2[j]
			comb.append( [int(l) for l in frag] )
			j+=1
		else:
			fragment = lengths[k]
			comb.append(range(int(fragment[0]),int(fragment[1])+1))
			k+=1

	# index ss elements for reading motif
	nl=0; nh=0; ne=0
	s_index=[]
	for s in ss:
		if s=='L':
			nl+=1
			s_index.append(s+'%s' %(nl))
		if s=='H':
			nh+=1
			s_index.append(s+'%s' %(nh))
		if s=='E':
			ne+=1
			s_index.append(s+'%s' %(ne))

	#----------------
	# read motifile
	#----------------
	motifile = open(motif_filename)
	header = motifile.readline()
	motif = motifile.readline()
	dic_motif={}
	for a,b in zip(header.split(),motif.split()):
		dic_motif[a] = b

	dic_abego={}
	for k,s in enumerate(s_index):
		if s in dic_motif.keys():
			motif_value = dic_motif[s]
			if motif_value.isdigit():
				comb[k]=[int(motif_value)]
			else:
				comb[k]=[len(motif_value)]
				dic_abego[s] = motif_value


	combinations = list(itertools.product(*comb))
	print 'Number of combinations: %s' %(len(combinations))

	return ss,combinations, dic_abego

def GetCombinations(topol):
        """
        Generates a list with all combinations of lengths of secondary structure elements
        """
	elements = re.compile("[HEL]")
	ss = elements.findall(topol)
	relengths = re.compile("(\d+)-(\d+)")
	relengths2 = re.compile("(\d+),(\d+)") # for specific lengths, not ranges
	lengths = relengths.findall(topol)
	lengths2 = relengths2.findall(topol)
	print lengths2
	## for interpreting specific lengths
	index=-1
	specific=[]
	for st in topol:
		if st=='[':
			index+=1
		if st==',':
			specific.append(index)

	comb=[] ; j=0 ; k=0
	for i in range(len(ss)):
		if i in specific:
			frag = lengths2[j]
			comb.append( [int(l) for l in frag] )
			j+=1
		else:
			fragment = lengths[k]
			comb.append(range(int(fragment[0]),int(fragment[1])+1))
			k+=1

	combinations = list(itertools.product(*comb))
	print 'Number of combinations: %s' %(len(combinations))
	return ss,combinations

def HBondConstraints(donor,acceptor):
    """
    Returns string for backbone H-bond constraints between donor residue and
    acceptor residue
    """
    hb_ang = np.deg2rad(160.)
    hb_ang_tol=np.deg2rad(20.0)
    #hb_ang_tol=np.deg2rad(30.0)
    st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor)
    #st = "AtomPair N %i O %i BOUNDED 2.7 3.3 0.2 \n" %(donor,acceptor)
    st+= "Angle N %i H %i O %i HARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
    st+= "Angle H %i O %i C %i HARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
    return st

def CircularHBondConstraints(donor,acceptor):
    """
    Returns string for backbone H-bond constraints between donor residue and
    acceptor residue, with angles angles around 360 and 0 being counted as a
    continuum
    """
    hb_ang = np.deg2rad(180.)
    #hb_ang_tol=np.deg2rad(20.0) -> original
    hb_ang_tol=np.deg2rad(30.0)
    st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor)
    #st = "AtomPair N %i O %i BOUNDED 2.7 3.3 0.2 \n" %(donor,acceptor)
    st+= "Angle N %i H %i O %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
    st+= "Angle H %i O %i C %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
    return st

def PairConstraints(a,b,value,tol,tag):
    """
    Returns properly formated string for declaring flat-bottomed distance constraints
    between CA carbons of resiedues a and b, at distance "value", and "tol" flat bottom
    width.
    """
    st = "AtomPair CA %i CA %i BOUNDED %3.1f %3.1f %3.1f 0.5 %s\n" %(a,b,value-tol,value+tol,tol/2,tag)
    return st

def HarmonicPairConstraints(a,b,value,sd):
    """
    Returns properly formated string for declaring harmonic distance constraints
    between CA carbons of resiedues a and b, at distance "value", and "sd" standard deviation.
    """
    st = "AtomPair CA %i CA %i HARMONIC %3.1f %3.1f\n" %(a,b,value,sd)
    return st


def AngleConstraints(a,b,c,value,tol,tag):
    """
    Returns properly formated string for declaring flat-bottomed angle constraints
    between CA carbons of resiedues a, b and c, at angle "value", and "tol" flat bottom
    width.
    """
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Angle CA %i CA %i CA %i BOUNDED %3.1f %3.1f 0.5 %s\n" %(a,b,c,ang-ang_tol,ang+ang_tol,tag)
    return st

def HarmonicAngleConstraints(a,b,c,value,tol):
    """
    Returns properly formated string for declaring harmonic angle constraints
    between CA carbons of resiedues a, b and c, at angle "value", and "tol" standard
    deviation. An angle of 360 is treated as being next to 0 (CIRCULARHARMONIC).
    """
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Angle CA %i CA %i CA %i CIRCULARHARMONIC %3.1f %3.1f\n" %(a,b,c,ang,ang_tol)
    return st


def CstTypeAngleConstraints(a,b,c,value,tol,cst_type):
    """
    Returns properly formated string for declaring 'harmonic' or 'bounded'
    angle constraints between CA carbons of resiedues a, b and c, at angle "value",
    and "tol" standard deviation (harmonic) or flat bottom width (bounded).
    """
    st=''
    if cst_type=='harmonic':
	    ang = np.deg2rad(value)
	    ang_tol = np.deg2rad(tol)
	    st = "Angle CA %i CA %i CA %i HARMONIC %3.1f %3.1f\n" %(a,b,c,ang,ang_tol)
    elif cst_type=='bounded':
            ang = np.deg2rad(value)
            ang_tol = np.deg2rad(tol)
	    sd = ang_tol/2
	    tag="ang_%i.%i.%i" %(a,b,c)
	    st = "Angle CA %i CA %i CA %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,b,c,ang-ang_tol,ang+ang_tol,sd,tag)
    return st

def DihedralConstraints(a,b,value,tol,tag):
    """
    Returns properly formated string for declaring flat-bottomed (bounded)
    dihedral angle constraints between CBa-CAa-CAb-CBb carbons of resiedues
    a and b at angle "value", and "tol" flat bottom width.
    """
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    sd = ang_tol/2
    st = "Dihedral CB %i CA %i CA %i CB %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,a,b,b,ang-ang_tol,ang+ang_tol,sd,tag)
    return st

def CaDihedralConstraints(a,b,c,d,value,tol):
    """
    Returns properly formated string for declaring flat-bottomed (bounded)
    dihedral angle constraints between CA carbons of resiedues
    a, b, c and d at angle "value", and "tol" flat bottom width.
    """
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    sd = ang_tol/2
    tag="ang_%i.%i.%i.%i" %(a,b,c,d)
    st = "Dihedral CA %i CA %i CA %i CA %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,b,c,d,ang-ang_tol,ang+ang_tol,sd,tag)
    return st

def CircularHarmonicCaDihedralConstraints(a,b,c,d,value,tol):
    """
    Returns properly formated string for declaring harmonic
    dihedral angle constraints between CA carbons of resiedues
    a, b, c and d at angle "value", and "tol" standard deviation.
    Angles of 0 and 360 are treated as being next to each other.
    """
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CA %i CA %i CA %i CA %i CIRCULARHARMONIC %3.1f %3.1f\n" %(a,b,c,d,ang,ang_tol)
    return st

def HarmonicDihedralConstraints(a,b,value,tol):
    """
    Returns properly formated string for declaring harmonic
    dihedral angle constraints between CBa-CAa-CAb-CBb carbons of resiedues
    a and b at angle "value", and "tol" standard deviation.
    """
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CB %i CA %i CA %i CB %i HARMONIC %3.1f %3.1f\n" %(a,a,b,b,ang,ang_tol)
    return st

def CircularHarmonicDihedralConstraints(a,b,value,tol):
    """
    Returns properly formated string for declaring harmonic
    dihedral angle constraints between CBa-CAa-CAb-CBb carbons of resiedues
    a and b at angle "value", and "tol" standard deviation.
    Angles of 0 and 360 are treated as being next to each other.
    """
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CB %i CA %i CA %i CB %i CIRCULARHARMONIC %3.1f %3.1f\n" %(a,a,b,b,ang,ang_tol)
    return st

def PerfectHelixCst(bluefile,helixn):
    """
    Returns properly formated string for declaring harmonic dihedral angle
    constraints for PHI and PSI angles of all residues of helix number 'helixn',
    forcing them in the most abundant helical conformation (straight helix).
    It leaves out the PHI of the first helix residue and the PSI os the last one.
    """
    blue = Blueprint(bluefile)
    blue.reindex_blueprint(start=1)
    hnum = 'H'+str(helixn)
    H = blue.segment_dict[hnum]
    HN = int(H.bp_data[0][0])
    HC = int(H.bp_data[-1][0])
    cst = ''
    #cst += "Dihedral N %i CA %i C %i N %i HARMONIC -0.78 0.09 \n" %(HN,HN,HN,HN+1) #psi
    cst += "Dihedral N %i CA %i C %i N %i HARMONIC -0.78 0.87 \n" %(HN,HN,HN,HN+1) #psi
    for pos in range(HN+1,HC-1):
		    #cst += "Dihedral C %i N %i CA %i C %i HARMONIC -1.05 0.09 \n" %(pos-1,pos,pos,pos) #phi
		    #cst += "Dihedral N %i CA %i C %i N %i HARMONIC -0.78 0.09 \n" %(pos,pos,pos,pos+1) #psi
            cst += "Dihedral C %i N %i CA %i C %i HARMONIC -1.05 0.87 \n" %(pos-1,pos,pos,pos) #phi
            cst += "Dihedral N %i CA %i C %i N %i HARMONIC -0.78 0.87 \n" %(pos,pos,pos,pos+1) #psi
            #cst += "Dihedral C %i N %i CA %i C %i HARMONIC -1.05 0.09 \n" %(HC-1,HC,HC,HC) #phi
    cst += "Dihedral C %i N %i CA %i C %i HARMONIC -1.05 0.87 \n" %(HC-1,HC,HC,HC) #phi
    return cst



def AmbiguousCst(cst_lst):
        """
        Returns properly formated string for declaring ambiguous constraints,
        it takes a list of Rosetta-formated constraints string.
        """
        header = 'AmbiguousConstraint\n'
        for cst in cst_lst:
                header += cst
        header += 'END_AMBIGUOUS\n'
        return header


def RegularStrandCurvature(**kwargs):
        """
        Returns angle and dihedral constraints strings that dictate the curvature of strands
        without bulges. Given a strand, place curvature and twist constraints every 2*level positions.
        Angle constraints, for curvature, are placed on CA atoms. Dihedral constraints are placed between
        CA-CB vectors. Values for twist and curvature are set by the user.
        """
        segment = kwargs.get("strand")
        level = kwargs.get("level") # 1 or 2
        bend = kwargs.get("global_bend",None)
        bend_tol = kwargs.get("global_bend_tol",10.0)
        twist = kwargs.get("global_twist",None )
        twist_tol = kwargs.get("global_twist_tol",5.0 )

        bend_area = kwargs.get("bend_area_value",None )
        bend_area_value = kwargs.get("bend_area_value",None )
        bend_area_value_tol = kwargs.get("bend_area_value_tol",5.0)  # n-ter, middle, c-ter

        twist_area = kwargs.get("twist_area")  # n-ter, middle, c-ter
        twist_area_value = kwargs.get("twist_area_value",None )
        twist_area_value_tol = kwargs.get("twist_area_value_tol",5.0 )

	bend_positions = kwargs.get("bend_positions",None )
	bend_positions_value = kwargs.get("bend_positions_value",None )
	bend_positions_value_tol = kwargs.get("bend_positions_value_tol",None )

	constraint_type = kwargs.get("constraint_type","harmonic" )

	blue = kwargs.get("blueprint")
	blue.reindex_blueprint(start=1)
	seg = blue.segment_dict[segment]
	cst_st=''
	#################
	n = len(seg.bp_data)

	############
	# BENDING
	############

	#----------------------
	# Set all triads for bend calculation
	#----------------------
	bend_triads = []
	step=level*2
	for i in range(0,n-step):
		 if i+step*2 < n:
			pos1 = i
			pos2 = i + step*1
			pos3 = i + step*2
		        bend_triads.append([pos1,pos2,pos3])

        #----------------------
        # Define bend areas
        #----------------------
	# By positions. This is used especially for positions paired to bulges, where we expect higher bend than the rest of triads
	bend_positions_triads=[]
	if bend_positions:
		for relpos in bend_positions:
			if relpos < 0: # when giving rel position from C-terminal (for E3)
				pos = n+relpos
			else:
				pos=relpos
			for triad in bend_triads:
				if pos==triad[1]: # central position of triad
					bend_positions_triads.append( triad )

	# By areas
	bend_area_triads=[]
        if bend_area =='n-term':
                bend_area_triads = [ bend_triads[0] ]
        elif bend_area =='c-term':
                bend_area_triads = [ bend_triads[-1] ]
        elif bend_area =='center':
                index =  len(bend_triads)/2 - 1
                if len(bend_triads) % 2 == 0:
                        bend_area_triads = [ pair for k in bend_triads[index:index+2] ]
                else:
                        bend_area_triads = [ bend_triads[index] ]

        ############
        # TWIST
        ############

      	#----------------------
       	# Set all pairs for twist calculation
	#----------------------
	twist_pairs = []
	step=level*2
	for i in range(0,n-step):
		if i+step < n:
			pos1 = i
			pos2 = i + step*1
			twist_pairs.append([pos1,pos2])

	#----------------------
	# Define twist areas
	#----------------------
	twist_area_pairs=[]
	if twist_area =='n-term':
		twist_area_pairs = [ twist_pairs[0] ]
	elif twist_area =='c-term':
		twist_area_pairs = [ twist_pairs[-1] ]
        elif twist_area =='center':
		index =  len(twist_pairs)/2 - 1
		if len(twist_pairs) % 2 == 0:
			twist_area_pairs = [ pair for pair in twist_pairs[index:index+2] ]
		else:
			twist_area_pairs = [ twist_pairs[index] ]

	#########################
	# INTRODUCE CONSTRAINTS
	#########################
    	# After indentifying combination of positions for bending and twist... Put constraints
	# Calculate Bends
	for triad in bend_triads:
		a,b,c = triad
		pos1 = seg.bp_data[b][0]
		pos2 = seg.bp_data[a][0]
		pos3 = seg.bp_data[c][0]
		if bend_positions and triad in bend_positions_triads:
			print triad
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend_positions_value,bend_positions_value_tol,constraint_type) ; cst_st += st
		elif bend_area_value and triad in bend_area_value_triads and triad not in bend_positions_triads:
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend_area_value,bend_area_value_tol,constraint_type) ; cst_st += st
		elif bend:
			print triad
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend,bend_tol,constraint_type) ; cst_st += st

	# Calculate Twists
	for pair in twist_pairs:
		a,b = pair
		pos1 = seg.bp_data[a][0]
		pos2 = seg.bp_data[b][0]
		if twist_area_value and pair in twist_area_pairs:
			st = DihedralConstraints(pos1,pos2,twist_area_value,twist_area_value_tol,'dih%i.%i' %(pos1,pos2)) ; cst_st += st
		elif twist:
			st = DihedralConstraints(pos1,pos2,twist,twist_tol,'dih%i.%i' %(pos1,pos2)) ; cst_st += st

	return cst_st

def AbegoPhiPsiConstraints(pos,blue):
    """
    Returns angle constraints for the Phi and Psi angles of position "pos",
    limiting it to the ABEGO bin designated in the Blueprint "blue"
    """
    abego = blue.bp_data[pos-1][2][1]
    if abego=="A":
        phi_min = -180.0
        phi_max = 0.0
        psi_min = -75.0
        psi_max = 50.0
    if abego=="B":
        phi_min = -180.0
        phi_max = 0.0
        psi_min = 50.0
        psi_max = 285.5
    if abego=="E":
        phi_min = 0.0
        phi_max = 180.5
        psi_min = 100.0
        psi_max = 260.5
    if abego=="G":
        phi_min = 0.0
        phi_max = 180.5
        psi_min = -100.0
        psi_max = 100.0
    if abego=="S":
        phi_min = -180.0
        phi_max = 0.0
        psi_min = 100.0
        psi_max = 195.0
    if abego=="P":
        phi_min = -180.0
        phi_max = 0.0
        psi_min = 100.0
        psi_max = 195.0
    if abego=="D":
        phi_min = -180.0
        phi_max = 0.0
        psi_min = 195.0
        psi_max = 285.5
    if abego=="Z":
        phi_min = -180.0
        phi_max = -100.0
        psi_min = 50.0
        psi_max = 100.0
    if abego=="Y":
        phi_min = -100.0
        phi_max = 0.0
        psi_min = 50.0
        psi_max = 100.0
    if abego=="M":
        phi_min = -180.0
        phi_max = -90.0
        psi_min = -75.0
        psi_max = -50.0
    if abego=="N":
        phi_min = -90.0
        phi_max = 0.0
        psi_min = -75.0
        psi_max = -50.0

    st=''
    phi_min = np.deg2rad(phi_min)
    phi_max = np.deg2rad(phi_max)
    psi_min = np.deg2rad(psi_min)
    psi_max = np.deg2rad(psi_max)
    if pos >= 2:
	    st+="Dihedral C %i N %i CA %i C %i BOUNDED %3.1f %3.1f 0.5 phi_%s\n" %(pos-1,pos,pos,pos,phi_min,phi_max,pos)
    if pos < len(blue.bp_data):
	    st+="Dihedral N %i CA %i C %i N %i BOUNDED %3.1f %3.1f 0.5 psi_%s\n" %(pos,pos,pos,pos+1,psi_min,psi_max,pos)
    return st

def AbegoConstraintSegment(segments,blue):
        """
        Returns angle constraints for the Phi and Psi angles of all positions in the
        input segment, constraining them to the designated ABEGO bin in the Bleuprint "blue".
        """
	blue.reindex_blueprint(start=1)
	cst=''
	for segment in segments:
		seg = blue.segment_dict[segment]
		for res in seg.bp_data:
			pos = res[0]
			if len(blue.bp_data[pos-1][2])==2: # abego specified.
				st = AbegoPhiPsiConstraints(pos,blue)
				cst+=st
	return cst

def HairpinPairingResidues(blue,segment1,segment2):
    """
    Returns a list of paired positions in strands designated by "segment1" and
    "segment2"
    """
    # segment1 is the first strand in the hairpin according to sequence
    s1 = blue.segment_dict[segment1]
    s2 = blue.segment_dict[segment2]
    map1={} ; map2={}
    pairs=[]
    if Bulged(s1):
        b1pos = Bulged(s1)
        shift = int(s1.bp_data[-1][0]) - b1pos
        for i in range(len(s2.bp_data)):
            pos2 = int(s2.bp_data[i][0]) # Antiparallel
            if i < shift:
                pos1 = int(s1.bp_data[-1-i][0])
            elif i+1 < len(s1.bp_data):
                pos1 = int(s1.bp_data[-1-i-1][0])
            else:
                break
            pairs.append([pos1,pos2])
    # The bulgepos+1 is the one paired to the second strand

    if Bulged(s2):
        b2pos = Bulged(s2)
        shift = b2pos - int(s2.bp_data[0][0])
        for i in range(len(s1.bp_data)):
            pos1 = int(s1.bp_data[-1-i][0])
            if i < shift:
                pos2 = int(s2.bp_data[i][0])
            elif i+1 < len(s2.bp_data):
                pos2 = int(s2.bp_data[i+1][0])
            else:
                break
            pairs.append([pos1,pos2])

    if Bulged(s1)==None and Bulged(s2)==None:
         for i in range(len(s2.bp_data)):
            pos2 = int(s2.bp_data[i][0]) # Antiparallel
            if i < len(s1.bp_data):
                pos1 = int(s1.bp_data[-1-i][0])
            else:
                break
            pairs.append([pos1,pos2])


    return pairs

def HbondsBulgedStrand(**kwargs):
        """
        Returns a list of Rosetta constraint strings describing hydrogen Hbonsds between
        two paired strands when one of them is bulged.
        """
	strand1 = kwargs.get('strand1')
	strand2 = kwargs.get('strand2')
	blue = kwargs.get('blueprint')

	pair1 = HairpinPairingResidues(blue,strand1,strand2)

	seg1 = blue.segment_dict[strand1]
	seg2 = blue.segment_dict[strand2]

	if Bulged(seg1):
		seg = seg1
	elif Bulged(seg2):
		seg = seg2
	else:
		print 'Warning: None of the strands is bulged'

	b1pos = Bulged(seg)

	hblist=[]
	# to the left of the bulge
	i=0
	pos = b1pos
	while pos > seg.bp_data[0][0]:
	    pos = (b1pos-2)-2*i
	    i+=1
	    hblist.append(pos)

	pos = b1pos
	# to the right of the bulge
	i=0
	while pos < seg.bp_data[-1][0]:
	    pos = (b1pos+1)+2*i
	    i+=1
	    hblist.append(pos)

	hbpairs=[]
	for pair in pair1:
	    hbpair = list( set(hblist).intersection(set(pair)) )
	    if len(hbpair)==1:
	        hbpairs.append(pair)

	cst_st = ''
	for hbpair in hbpairs:
		pos1,pos2 = hbpair
		st = CircularHBondConstraints(pos1,pos2) ; cst_st += st
		st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
		# Add hbond for bulge position (which is not C-alpha paired)
		if (b1pos+1) in hbpair: # position b1pos+1 is in the pairing list (not the bulgepos, so we add the hbond here)
			# identify paired position
			paired_pos = list( set(hbpair).difference(set([b1pos+1])) )[0]
			st = CircularHBondConstraints(b1pos,paired_pos) ; cst_st += st


	return cst_st

def AllSheetSegmentPairs(blue):
        """
        Returns a dictionary where each position in a pair maps to the strand segment
        it's paird to, and that segment pairs to the position in itself where the initial
        position pairs to.
        """
        cst_st=''
        pair1 = HairpinPairingResidues(blue,'E1','E2')
        pair2 = HairpinPairingResidues(blue,'E2','E3')
        pair3 = HairpinPairingResidues(blue,'E3','E4')

        pairs=[]
        pairs.extend(pair1)
        pairs.extend(pair2)
        pairs.extend(pair3)

	dic_pairs={}
	for pair in pairs:
		a,b = pair
		seg_a = blue.residue_segment(a)
		seg_b = blue.residue_segment(b)
		dic_pairs.setdefault(a,{})
		dic_pairs.setdefault(b,{})
		dic_pairs[a][seg_b]=b
		dic_pairs[b][seg_a]=a

	return dic_pairs


def HbondsRegularHairpin(**kwargs):
        """
        Returns a list of Rosetta constraint strings describing hydrogen Hbonsds between
        two paired strands when none of them is bulged.
        """
        strand1 = kwargs.get('strand1')
        strand2 = kwargs.get('strand2')
        blue = kwargs.get('blueprint')

        pair1 = HairpinPairingResidues(blue,strand1,strand2)

        seg1 = blue.segment_dict[strand1]
        seg2 = blue.segment_dict[strand2]


        # to the left of the bulge
        inipos = seg1.bp_data[-1][0] # start counting from hairpin loop
        pos=inipos
        hblist=[pos]
        # all residues of seg1 are hbonded to seg2
        while pos >= seg1.bp_data[0][0]+2:
            pos -= 2
            hblist.append(pos)

        hbpairs=[]
        for pair in pair1:
            hbpair = list( set(hblist).intersection(set(pair)) )
            if len(hbpair)==1:
                hbpairs.append(pair)

        cst_st = ''
        for hbpair in hbpairs:
                pos1,pos2 = hbpair
                st = CircularHBondConstraints(pos1,pos2) ; cst_st += st
                st = CircularHBondConstraints(pos2,pos1) ; cst_st += st


        return cst_st

if __name__ == '__main__':
	print "Functions for building curved beta sheets"
