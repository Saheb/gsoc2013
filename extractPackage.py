#!/usr/bin/env python
# extractPackage
#
# Part of OGDF
# (c) 2007 - 2012
# Author: Markus Chimani, markus.chimani@cs.uni-dortmund.de
#########################################################

import os, sys, shutil
from datetime import date

if len(sys.argv) != 3:
	print '===================================================='
	print 'extractPackage'
	print '--------------'
	print '   Use to generate a distributable package'
	print '   This tool is part of OGDF'
	print ' '
	print '   Usage: python extractPackage.py <TAG> <TARGETDIR>'
	print '     <TAG> = ( lic_gpl | lic_com )'
	print '     <TARGETDIR> = the new base directory'
	print '===================================================='
	sys.exit(1)

# parameters
checkTag = sys.argv[1]
targetDir = sys.argv[2]

svnCommand = 'svn propget ' + checkTag + ' '

def Copy( fullname ):
	runArgs[len(runArgs)-1] = fullname
	os.spawnv( os.P_WAIT, runCommand, runArgs)

def DirWalk( curdir, exportTo ):
	names = os.listdir( curdir )
	names.sort()
	createdDir = False

	for name in names:
		if name.startswith('.') or name.startswith('_') or (name=='Debug') or (name=='Release') or (name=='html'):
			continue
		fullname = os.path.normpath(os.path.join(curdir, name))
		if os.path.isdir(fullname) and not os.path.islink(fullname):
			objs = DirWalk( fullname, exportTo+name+'/' )
		else:
			f = os.popen(svnCommand+fullname, 'r')
			ch = f.read(1)
			if ch:
				if not createdDir:
					try:
						os.makedirs(exportTo)
					except:
						pass
					createdDir = True
				shutil.copy(fullname,exportTo)
			f.close()

# Call recursive function
DirWalk( '.', targetDir+'/')

# Fix names and permissions of files
os.rename(targetDir+'/makeMakefile.config.default',targetDir+'/makeMakefile.config')
os.rename(targetDir+'/makeVCProj.config.default',targetDir+'/makeVCProj.config')
os.rename(targetDir+'/makeVCXProj.config.default',targetDir+'/makeVCXProj.config')
os.chmod(targetDir+'/makeMakefile.sh',484)

# Generate version file
version_h = open(targetDir + '/include/ogdf/internal/version.h','w')
version_h.write('#ifndef OGDF_VERSION_H\n')
version_h.write('#define OGDF_VERSION_H\n\n')
version_h.write('#define OGDF_VERSION "' + date.today().strftime('%Y.%m') + '"\n\n')
version_h.write('#endif // OGDF_VERSION_H\n')
version_h.close()
