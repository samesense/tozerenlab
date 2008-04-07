#!/usr/bin/env python

# txtIt.py
# given a subject line &
# body,
# send a txt message to my phone
# from my gmail account

import smtplib, sys

def txt(subject, body):
    sender = 'perry'
    
    # me
    to = '2152790262@txt.att.net'
    # dad
    # to = '2147080561@txt.att.net'
    # janie
    # christin
        
    s = smtplib.SMTP('smtp.gmail.com', 25)
    s.ehlo()
    s.starttls()
    s.ehlo()
    # very insecure
    s.login('samesense', 'PKg9J1=')
    headers = "From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n" % \
	      (sender, to, subject)
    message = headers + body
    s.sendmail(sender, to, message)
    s.quit()

def main():
    txt(sys.argv[1], sys.argv[2])
    sys.exit(0)

if __name__=='__main__':
    if len(sys.argv) != 3:
        print 'enter subject & body'
        sys.exit(0)
    else:
	main()


