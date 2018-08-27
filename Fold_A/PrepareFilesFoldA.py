""
#Script for writing blueprint files from a string defining the protein topology
#H[n-m]L[n-m]
#Length: xxx
""
# USE: python build_blueprints.v2.py -xml template_bb+design.xml -blueresfile

import sys
sys.path.append("../")

import itertools
from argparse import ArgumentParser
import re
import os
import copy
from Blueprint import Blueprint
import numpy as np
from BPB_functions import *
import subprocess

#==============================
# INPUT PARAMETERS
#==============================
parser = ArgumentParser()
parser.add_argument('-xml', type=str, help="xml template")
parser.add_argument('-pdb', type=str, help="input pdb")
parser.add_argument('-wts', type=str, help="wts file for design")
parser.add_argument('-blueresfile', action='store_true')
parser.add_argument('-resfile', type=str)
args = parser.parse_args()

template_xml = args.xml
topol = "L[1-1]H[20-20]L[3-3]E[10-10]L[2-2]E[10-10]L[2-2]E[12-12]L[2-2]E[9-9]L[1-1]"

ss,combinations = GetCombinations(topol)

print topol

bulges =  {'E1':5,'E4':5} # Specific definition. Bulges are located at positions defined by the user


# Make directory and bluprints for each combination along with bulge positions
for comb in combinations:
	## Build directories and add bulge combinations
	#---------------------------------------------------------
	# Make the directory name
	filename = 'r0-all_rev-' ; strand=0 ; resnum=0
	for i,s in enumerate(ss):
		filename+='%s%i' %(ss[i],comb[i])

	# First build a plain blueprint file to after create a bp object to make easier the bulge combinations

	MakePlainBlueprint(ss,comb,'bp')
	blue = Blueprint('bp')

	# Bulges
	keys = bulges.keys() # bulged strand names
	keys.sort()
	bpos_dic = {} # all positions considered for each bulged strand
	for key in keys:
		for seg in blue.segments:
			if seg.id == key:
				st_len = len(seg.bp_data)
				if bulges[key] == 'n-center':
					if st_len %2 == 0: # even strand length. Bulge must be at odd position
						bposs = range(1,st_len,2)[1:-1]
					else:
						bposs = range(2,st_len,2)[1:-1]
				elif bulges[key] == 'c-center': # it depends on the orientation of the first residue. In this case...
					bposs = range(1,st_len,2)[1:-1]
				else:
					bposs=[bulges[key]]

				bpos_dic[key] = bposs

	# Take all bulge combinations
	bcomb=[]
	for key in keys:
	        bcomb.append( bpos_dic[key] )

	bcombinations = list(itertools.product(*bcomb))

	for j,bulcomb in enumerate(bcombinations):
			picked_bulges={}
			pathname=filename
			for k,key in enumerate(keys):
				pathname += '-b%s.%i' %(key,bulcomb[k]) # this is the new filename
				picked_bulges[key]=bulcomb[k]

			# Make directory of the topology
			if not os.path.exists(pathname):
				os.mkdir(pathname)
			os.chdir(pathname)

			## Build blueprints
			MakeRefBlueprint(ss,comb,picked_bulges,refblue = 'bp')
			# split blueprint according to different phases in backbone generation
			# pairing according to current blueprint
			###########
			# ATTENTION!!! to make more efficient the biasing with the SSPAIR
			# headers of blueprints it is better only write the pairing that is
			# actually built with this blueprint, then 100% of biasing is focused on the region being constructed.
			# Pairings must be indexed according to the current blueprint not the global one in contrast to segments
			#########
			# step1
			MakeFirstBlueprint(refblue = 'bp', segments = ['E2','L4','E3'], \
					   newblue = 'bp1', ss_pairing={'E1':['E2.A']},adapt={'E3':'E2'})

			# step2
			AddSegmentToBlueprint(refblue = 'bp', segments = ['E1','L3'], blue0 = 'bp1', newblue='bp2', \
						      append=False, ss_pairing_shift="SSPAIR 1-2.A.1", ss_pairing={'E1':['E2.A'],'E2':['E3.A']})
			# step3
			shift = Shift(refblue = 'bp', seg1 = 'E3', seg2 = 'E2')
			AddSegmentToBlueprint(refblue = 'bp', segments = ['L5','E4'], \
					      blue0 = 'bp2', newblue='bp3', append=True, \
					      ss_pairing={'E1':['E2.A'],'E2':['E3.A'],'E3':['E4.A']},insert={'E3':shift})

			# step4
			AddSegmentToBlueprint(refblue = 'bp', segments = ['H1','L2'], \
					      blue0 = 'bp3', newblue='bp4', append=False, ss_pairing={'E1':['E2.A'],'E2':['E3.A'],'E3':['E4.A']}, \
					      seg_abego={'L2':'GBA'}) # Added B to L6

			#---------------------------------------------------------
			# Blueprints are already written. Now XMLs, cst, etc...
			# Write constraints for good HB pairing of built strands

			# Write initial input PDB:
			if args.pdb:
				os.system('cp ../%s input.pdb' %(args.pdb))
			else:
				write_dummy_pdb('input.pdb')

			if args.wts:
				os.system('cp ../%s .' %(args.wts))
			# take the sheet pdb reset to one
			# RESFILE
			if args.resfile:
				os.system( 'python ../manual_resfile.py bp4.b %s' %( args.resfile ) )

			# XML
		    	os.system('cp ../%s foo.xml' %(template_xml))
			xml_lines = open('../%s' %(template_xml),'r').readlines()
   			# names of blueprints must be consistent between here and the template.xml
			# Move above dir for another topology

			################################################
			## Define bulged hairping shifts to find hbond partner
			blue = Blueprint('bp4') ; blue.reindex_blueprint(start=1) # Get the global blueprint.
			# First bulged hairpin: E3-E4. I define it as a "N-bulged hairpin"
			# For a N-bulged hairpin
                        bulged_strand = 'E1'
                        seg = blue.segment_dict[bulged_strand]
                        bulge1_shift = ( len(seg.bp_data) - bulges[bulged_strand] ) - 1

                        # Second bulged hairpin: E5-E6. I define it as a "C-bulged hairpin"
                        # For a C-bulged hairpin
                        bulged_strand = 'E4'
                        seg = blue.segment_dict[bulged_strand]
                        bulge2_shift = bulges[bulged_strand]

			#----------------------
			# Write flags file
			#----------------------
			flags_out = open('flags','w')
			################################################

			# step1
			blue = Blueprint('bp1') ; blue.reindex_blueprint(start=1)
			fileout = open('cst1','w') # E2//E3
			fileout_min = open('cst1_min','w') # E2//E3
                        #st11 = RegularStrandCurvature(level=3,blueprint=blue,strand='E1',\
                        #global_twist=40.0, global_twist_tol=30.0 )
                        #st22 = RegularStrandCurvature(level=3,blueprint=blue,strand='E2',\
                        #global_twist=50.0, global_twist_tol=35.0 )


                        st1 = RegularStrandCurvature(level=2,blueprint=blue,strand='E1',\
                        global_bend=55.0,global_bend_tol=10.0,\
			global_twist=35.0,global_twist_tol=10.0,\
			constraint_type='bounded' )

                        st2 = RegularStrandCurvature(level=2,blueprint=blue,strand='E2',\
			global_twist=35.0,global_twist_tol=10.0,\
                        global_bend=55.0,global_bend_tol=10.0,constraint_type='bounded' )

			# Vars for Stage #1
			#SS_s1 = "-parser:script_vars s1_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp1"])) ; flags_out.write(SS_s1)
			SS_s1 = "-parser:script_vars s1_SS=%s\n" %(''.join([ pos[2][0] for pos in Blueprint('bp1').bp_data])) ; flags_out.write(SS_s1)
			bp1_name = "-parser:script_vars bp1=../bp1.b\n" ; flags_out.write(bp1_name)
			bp1b_name = "-parser:script_vars bp1.b=../bp1.b\n" ; flags_out.write(bp1b_name)
			bp1 = Blueprint('bp1'); r_res = [ i+1 for i,pos in enumerate(bp1.bp_data) if pos[-1] == 'R' ]
			R_in_bp1= "-parser:script_vars R_in_bp1=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp1)
			cst1_name = "-parser:script_vars cst1=../cst1\n" ; flags_out.write(cst1_name)
			cst1_min_name = "-parser:script_vars cst1_min=../cst1_min\n" ; flags_out.write(cst1_min_name)
			flags_out.write('-restore_talaris_behavior\n')
			sthb = HbondsRegularHairpin(strand1="E1",strand2="E2",blueprint=blue) ; fileout.write(sthb)

			fileout.write(st1) ; fileout.write(st2) # ; fileout.write(st11) ; fileout.write(st22)
			fileout_min.write(st1) ; fileout_min.write(st2)

			fileout.close()

			#############################################

                        # step2
                        fileout = open('cst2','w') # E1

			# bulge constraints
			blue = Blueprint('bp2.b') ; blue.reindex_blueprint(start=1)
			s1 = blue.segment_dict['E1']
			s2 = blue.segment_dict['E2']

			st = HbondsBulgedStrand(strand1="E1",strand2="E2",blueprint=blue) ; fileout.write(st)
                        # Specifiy movemap for minimization in XML
                        pos1 = s1.bp_data[0]
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='MoveMap2 Flexible',xxx=1, yyy=s2.bp_data[0][0]-1)
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='MoveMap2 Rigid',xxx=s2.bp_data[0][0], yyy=blue.bp_data[-1][0])
                        fileout.close()

			# Vars for Stage #2
			#SS_s2 = "-parser:script_vars s2_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp2"])) ; flags_out.write(SS_s2)
			SS_s2 = "-parser:script_vars s2_SS=%s\n" %(''.join([ pos[2][0] for pos in Blueprint('bp2').bp_data])) ; flags_out.write(SS_s2)
			bp2_name = "-parser:script_vars bp2=../bp2.b\n" ; flags_out.write(bp2_name)
			bp2b_name = "-parser:script_vars bp2.b=../bp2.b\n" ; flags_out.write(bp2b_name)
                        bp2 = Blueprint('bp2'); r_res = [ i+1 for i,pos in enumerate(bp2.bp_data) if pos[-1] == 'R' ]
			R_in_bp2= "-parser:script_vars R_in_bp2=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp2)
			cst2_name = "-parser:script_vars cst2=../cst2\n" ; flags_out.write(cst2_name)

                        #############################################

                        # step3 # in this test i dont reead cst3 in xml.  try without these new constraints only for E4
                        fileout = open('cst3','w') # E4

                        # bulge constraints
                        blue = Blueprint('bp3.b') ; blue.reindex_blueprint(start=1)
                        s1 = blue.segment_dict['E1']
                        s2 = blue.segment_dict['E2']
                        s3 = blue.segment_dict['E3']
                        s4 = blue.segment_dict['E4']

			# hbonds
			st = HbondsBulgedStrand(strand1="E3",strand2="E4",blueprint=blue) ; fileout.write(st)
			st = HbondsRegularHairpin(strand1="E2",strand2="E3",blueprint=blue) ; fileout.write(st)

			# Inter-strand twist: E2-E3
			dic_pairs = AllSheetSegmentPairs(blue)
			p1 = s3.bp_data[-1][0]
			p2 = p1-4
			p3 = dic_pairs[p2]['E2']
			p4 = p3-2
			twist = -35.0 ; twist_tol=10.0
			#st = CircularHarmonicCaDihedralConstraints(p1,p2,p3,p4,twist,twist_tol) ; fileout.write(st)
			#st = CaDihedralConstraints(p1,p2,p3,p4,twist,twist_tol) ; fileout.write(st)
			################ DIHEDRAL CST FOR BETTER TWIST IN E5 and E6 ###########
			## 20140415 changed residues and angle for E3 to better reproduce native twist
			cE3 = s3.bp_data[-1][0]
			bulge_list = Bulges(blue)
                        bulge2 = bulge_list[1]
			#st = HarmonicDihedralConstraints(cE3-4,cE3,82,10) ; fileout.write(st)
			#st = HarmonicDihedralConstraints(bulge2,bulge2+4,125,10) ; fileout.write(st)
			######################################################################

                        # Specifiy movemap for minimization in XML
                        pos1 = s1.bp_data[0]
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='MoveMap3 Rigid',xxx=blue.bp_data[0][0], yyy=s2.bp_data[-2][0])
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='MoveMap3 Flexible',xxx=s2.bp_data[-1][0], yyy=blue.bp_data[-1][0])

                        fileout.close()

			# Vars for Stage #3
			#SS_s3 = "-parser:script_vars s3_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp3"])) ; flags_out.write(SS_s3)
			SS_s3 = "-parser:script_vars s3_SS=%s\n" %(''.join([ pos[2][0] for pos in Blueprint('bp3').bp_data])) ; flags_out.write(SS_s3)
			bp3_name = "-parser:script_vars bp3=../bp3.b\n" ; flags_out.write(bp3_name)
			bp3b_name = "-parser:script_vars bp3.b=../bp3.b\n" ; flags_out.write(bp3b_name)
			bp3 = Blueprint('bp3'); r_res = [ i+1 for i,pos in enumerate(bp3.bp_data) if pos[-1] == 'R' ]
			R_in_bp3= "-parser:script_vars R_in_bp3=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp3)
			cst3_name = "-parser:script_vars cst3=../cst3\n" ; flags_out.write(cst3_name)

			##### STEP 4 #####
                        blue = Blueprint('bp4')
			blue.reindex_blueprint(start=1)
                        bulge_list = Bulges(blue)
			bulge2 = bulge_list[1]
			e1 = blue.segment_dict['E1']
			e1_hb = e1.bp_data[0][0]
			seg1 = blue.segment_dict['H1']
			fileout = open('cst4','w')
			fileoutb = open('cst4b','w')
			fileoutc = open('cst4c','w')
			# 2 pairs of hydrogen bonding between E1 and E6

			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist4a',xxx=e1_hb, yyy=bulge2)

			## CSTs for H3-E3 loop:

			# Perfect Helices
			st = PerfectHelixCst('bp4',1); fileoutb.write(st); fileout.write(st);

			#Enforce ABEGO at positions 10 and 9:
			l2 = blue.segment_dict['L2']
			l2N = int(l2.bp_data[0][0])
			h1 = blue.segment_dict['H1']
			h1n = h1.bp_data[0][0] + 1
			l3 = blue.segment_dict['L3'].bp_data[0][0]
			h1c = h1.bp_data[-1][0]
			e4c = blue.segment_dict['E4'].bp_data[-1][0]
			st = HarmonicPairConstraints(h1n,l3-1,10,8) ; fileout.write(st); fileoutb.write(st);
			st = HarmonicPairConstraints(h1n,l3+2,10,8) ; fileout.write(st); fileoutb.write(st);
			st = HarmonicPairConstraints(h1c,e4c,8,8) ; fileout.write(st); fileoutb.write(st);
			st = HarmonicPairConstraints(h1c,e4c-2,8,8) ; fileout.write(st); fileoutb.write(st);
			fileout.close()
			fileoutb.close()
			fileoutc.close()

			# Vars for Stage #4
			SS_s4 = "-parser:script_vars s4_SS=%s\n" %(''.join([ pos[2][0] for pos in Blueprint('bp4').bp_data])) ; flags_out.write(SS_s4)
			#bp4_name = "-parser:script_vars bp4=../bp4\n" ; flags_out.write(bp4_name)
			bp4b_name = "-parser:script_vars bp4.b=../bp4.b\n" ; flags_out.write(bp4b_name)
			bp4 = Blueprint('bp4'); r_res = [ i+1 for i,pos in enumerate(bp4.bp_data) if pos[-1] == 'R' ]
                        R_in_bp4= "-parser:script_vars R_in_bp4=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp4)
			cst4_name = "-parser:script_vars cst4=../cst4\n" ; flags_out.write(cst4_name)
			cst4b_name = "-parser:script_vars cst4b=../cst4b\n" ; flags_out.write(cst4b_name)
			cst4c_name = "-parser:script_vars cst4c=../cst4c\n" ; flags_out.write(cst4c_name)
			flags_out.write("-picking_old_max_score 1\n")
			flags_out.write("-rama_map ../../../Rama_XPG_3level.txt\n")
			flags_out.close()
			#-----------------------
			# Write modified xml
			#-----------------------
			xml_out = open('input.xml','w')
			for line in xml_lines:
				xml_out.write(line)
			xml_out.close()

			os.chdir('../')
