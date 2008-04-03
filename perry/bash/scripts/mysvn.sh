#!/bin/sh
#Inspired by tinyurl.com/yt4ssq
 
# Scrub weird characters
MSG=`echo $@|tr ' ' '+'`
 
# Send twitter request
curl --basic --user "samesense:janie" --data-ascii "status=SVN CI: $MSG" "http://twitter.com/statuses/update.json"
 
# Send SVN request
svn ci -m $MSG