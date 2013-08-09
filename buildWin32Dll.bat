:: 07-08-2009
:: created by Sebastian Stein, sebastian.stein[at]tu-dortmund.de


@echo off
:: configure and create vcproj file to use dll as output
python makeVCProj.py config=makeVCProjDll.config.default

:: upgrade vcproj file if necessary
devenv ogdf.vcproj /Upgrade
devenv ogdf.vcproj /Clean
devenv ogdf.vcproj /Build Release

:: copy ogdf.dll into gryphon folder, which is necessary to
:: a) build the gryphon setup correctly
:: b) execute gryphon_edit.exe from the gryphon project folder (without installing it)

if exist Win32\Release\ogdf.dll (
	if exist ..\Gryphon\App copy "Win32\Release\ogdf.dll" "..\Gryphon\App"
)