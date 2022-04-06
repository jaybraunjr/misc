 #!/usr/bin/env python3
import argparse
import re
import os
import subprocess

def skipline(buff,j):
	while(j<len(buff) and buff[j] != '\n'):
		j+=1
	return j+1

def gmxGrompp(mdpfile,grofile,topfile,tprfile):
	myexec = 'gmx'
	args = [myexec, 'grompp', 
	        #'-quiet', 
	        '-f', mdpfile,
	        '-c', grofile,
	        '-p', topfile,
	        '-o', tprfile,
	        '-maxwarn', '1'
	        ]
	p = subprocess.Popen(args, 
		                 stdin=subprocess.PIPE,
		                 #stdout=subprocess.PIPE,
		                 #stderr=subprocess.PIPE
		                 )
	p.wait()
	os.remove('mdout.mdp')

def gmxMdrun(trrfile,tprfile,logfile):
	myexec = 'gmx'
	args = [myexec, 'mdrun', 
	        #'-quiet', 
	        '-o', trrfile,
	        '-s', tprfile,
	        '-g', logfile
	        ]
	p = subprocess.Popen(args, 
		                 stdin=subprocess.PIPE,
		                 #stdout=subprocess.PIPE,
		                 #stderr=subprocess.PIPE
		                 )
	p.wait()

def gmxEnergy(enerfile,xvgfile):
	myexec = 'gmx'
	args = [myexec, 'energy', 
	        #'-quiet', 
	        '-f', enerfile,
	        '-o', xvgfile
	        ]
	p = subprocess.Popen(args, 
		                 stdin=subprocess.PIPE,
		                 #stdout=subprocess.PIPE,
		                 #stderr=subprocess.PIPE
		                 )
	p.stdin.write(b'Potential\n')
	p.communicate()[0]
	p.stdin.close()
	p.wait()

def runDihedral(dihedral,args):
	tprfile = 'dihedral.tpr'
	trrfile = 'dihedral.trr'
	logfile = 'dihedral.log'
	topfile = args.t
	grofile = args.c
	mdpfile = args.m
	enerfile = args.e
	xvgfile = str(dihedral)+'.xvg'
		# Edit top file, dihedral restraint
	with open(topfile, 'r') as file:
		buff = file.read()
	j=buff.index('[dihedral_restraints]')
	j=skipline(buff,j)
	k=skipline(buff,j)
	fields=re.split(' +',buff[j:k])
	if fields[0]=='':
		fields=fields[1:]
	fields[5]=str(dihedral)
	line='   '.join(fields)
	line='  '+line
	buff=buff[:j]+line+buff[k:]
	with open(topfile, 'w') as file:
		file.write(buff)

	gmxGrompp(mdpfile,grofile,topfile,tprfile)
	os.remove('confout.gro')
	gmxMdrun(trrfile,tprfile,logfile)
	gmxEnergy(enerfile,xvgfile)

	os.remove(tprfile)
	os.remove(trrfile)
	os.remove(enerfile)
	os.remove(logfile)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', help='GROMACS (.gro) coordinate file', required=True)
	parser.add_argument('-m', help='GROMACS (.mdp) configuration file', required=True)
	parser.add_argument('-t', help='GROMACS (.top) topology   file', required=True)
	parser.add_argument('-u', help='GROMACS (.tpr) tpr   file', required=True)
	parser.add_argument('-e', help='GROMACS (.edr) energy   file', required=True)
	args = parser.parse_args()

	dihedralStart=0
	dihedralInc=10
	numInc=36
	dihedral=dihedralStart
	count=0
	while count<numInc:
		runDihedral(dihedral,args)
		args.c='confout.gro'
		dihedral+=dihedralInc
		count+=1



main()