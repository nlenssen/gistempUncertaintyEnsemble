#!/bin/zsh
RHOST=ftp.ncdc.noaa.gov
RUSER=anonymous
RPASSWORD=nlenssen@nasa.gov


gftp -inv $RHOST <<EOF
user $RUSER $RPASSWORD
binary
cd $1
get $2 $3
bye
EOF
