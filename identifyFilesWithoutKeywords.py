# indentifyFilesWithoutKeywords
#
# Part of OGDF
# (c) 2012
# Author: Markus Chimani, markus.chimani@uni-jena.de
#         Carsten Gutwenger, carsten.gutwenger@cs.tu-dortmund.de
################################################################

import os, sys, shutil

if len(sys.argv) != 1:
	print '===================================================='
	print 'indentifyFilesWithoutKeywords'
	print '--------------'
	print '   Use to generate a list of files without set license'
	print '   property. This tool is part of OGDF'
	print ' '
	print '   Usage: python indentifyFilesWithoutKeywords.py'
	print '===================================================='
	sys.exit(1)

svnCommand = 'svn propget svn:keywords'

def DirWalk( curdir ):
	names = os.listdir( curdir )
	names.sort()

	for name in names:
		if name.startswith('.'):
			continue

		fullname = os.path.normpath(os.path.join(curdir, name))
		if (fullname=='include/coin') or (fullname=='src/coin'):
			continue

		if os.path.isdir(fullname) and not os.path.islink(fullname):
			objs = DirWalk( fullname )
		else:
			f = os.popen(svnCommand + ' ' + fullname, 'r')
			ch = f.read(1)

			if not ch:
				print fullname

# Call recursive function
DirWalk( 'include' )
DirWalk( 'src' )

