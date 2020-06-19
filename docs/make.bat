@ECHO OFF 

pushd %~dp0 

REM Command file for Sphinx documentation 

set SPHINXBUILD = sphinx-build
set SOURCEDIR= source
set BUILDDIR= build

if "%1" == "" goto help


%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O% 
goto end 

:help 
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O% 

:end 

popd

