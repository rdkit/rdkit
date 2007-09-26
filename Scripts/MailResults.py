#!/usr/bin/env python
import sys, smtplib, base64, cStringIO, time, os.path
from email import Message, Utils

SMTP_HOSTNAME=''
SMTP_USER=''
SMTP_DEST="rdkit-devel@lists.sourceforge.net"
def sendLog(fileName,summaryName=""):

  msg = Message.Message()
  msg["To"]=SMTP_DEST
  msg["From"]=SMTP_USER
  msg["Subject"]='RDKitBuild: Nightly Build Results for %s'%time.strftime("%d/%m/%Y")
  msg["Date"] = Utils.formatdate(localtime=1)
  msg["Message-ID"] = Utils.make_msgid()
  msg["Mime-version"] = "1.0"
  msg["Content-type"]= "Multipart/mixed"
  msg.preamble="Mime message\n"
  msg.epilogue=""

  subMsg = Message.Message()
  subMsg["Content-type"]= "text/plain"
  subMsg["Content-transfer-encoding"]="7bit"
  summaryText="Automatically generated email"
  if summaryName:
    try:
      summaryData=open(summaryName,'r').read()
    except IOError:
      summaryText += "\n Could not open summary file"
    else:
      summaryText += "\n\n    TEST SUMMARY\n"+summaryData
  subMsg.set_payload(summaryText)

  msg.attach(subMsg)

  subMsg = Message.Message()
  subMsg.add_header("Content-type","application/x-gzip",name=os.path.basename(fileName))
  subMsg.add_header("Content-transfer-encoding","base64")
  body=cStringIO.StringIO()
  base64.encode(open(fileName, 'rb'), body)
  subMsg.set_payload(body.getvalue())

  msg.attach(subMsg)
  
  smtp = smtplib.SMTP(SMTP_HOSTNAME)
  smtp.sendmail(msg['From'],
                [msg['To']],
                msg.as_string())
  smtp.quit()

if __name__=="__main__":
  import sys
  fName = sys.argv[1]
  if len(sys.argv)>2:
    summaryName=sys.argv[2]
  else:
    summaryName = ""
  sendLog(fName,summaryName)
