#!/bin/sh

# Scrub weird characters
MSG=`echo $@|tr ' ' '+'`

# Send twitter request
curl --basic --user "samesense:janie" --data-ascii "status=$MSG" "http://twitter.com/statuses/update.json"