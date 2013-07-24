#!/usr/bin/env python
# extractPackage
#
# Part of OGDF
# (c) 2007
# Author: Markus Chimani, markus.chimani@cs.uni-dortmund.de
#########################################################

import os, sys, shutil

if len(sys.argv) != 1:
	print '===================================================='
	print 'indentifyFilesWithoutLIC'
	print '--------------'
	print '   Use to generate a list of files without set license'
	print '   property. This tool is part of OGDF'
	print ' '
	print '   Usage: python indentifyFilesWithoutLIC.py'
	print '===================================================='
	sys.exit(1)

svnCommand = 'svn propget '
tags = [ 'lic_gpl ', 'lic_com ', 'lic_no' ]

def DirWalk( curdir ):
	names = os.listdir( curdir )
	names.sort()

	for name in names:
		if name.startswith('.') or name.startswith('_') or (name=='Debug') or (name=='Release') or (name=='html'):
			continue
		fullname = os.path.normpath(os.path.join(curdir, name))
		if os.path.isdir(fullname) and not os.path.islink(fullname):
			objs = DirWalk( fullname )
		else:
			isLic = False
			for t in tags:
				with os.popen(svnCommand + ' ' + t + ' ' + fullname, 'r') as f:
					isLic = f.read(1)
				if isLic:
					break
			if not isLic:
				print fullname

# Call recursive function
DirWalk( '.' )

