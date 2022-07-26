#!/bin/csh -f
# tar qnifft and move to pub directory
#
#/bin/rm /hosts/calvin/usr1/people/ftp/pub/qnifft.tar /hosts/calvin/usr1/people/ftp/pub/Untar_qnifft.com
#tar cvf /hosts/calvin/usr1/people/ftp/pub/qnifft.tar \
tar cvf qnifft.tar README Install.com Tar.com Untar.com \
*.lis delphi.def data examples source/qnifft/v2.2 source/utility source/qnifft/v1.4 \
source/Makefile*
# cp Untar.com /hosts/calvin/usr1/people/ftp/pub/Untar_qnifft.com
# cp qnifft.tar /hosts/calvin/usr1/people/ftp/pub/Untar_qnifft.com
