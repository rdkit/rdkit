#!/usr/bin/env python
import os, sys, re

from roundup import instance, date, hyperdb

def create(db, **issue):
  # first we need to create a 'msg' node containing the message
  msg = {}
  # check if the message was passed on the command line
  if issue.has_key('messages'):
    # from command line
    msg['content'] = issue['messages']
  else:
    msg['content'] = ' No Message Provided '
  msg['author'] = db.getuid()
  msg['date']   = date.Date('.')

  # create the 'msg' node
  issue['messages'] = db.msg.create(**msg) # will be put in list in next part

  # resolve linked and multilinked properties
  properties = db.issue.getprops()
  for key in issue.keys():
    if properties.has_key(key):
      if isinstance(properties[key], hyperdb.Link):
        if not issue[key].isdigit():
          # resolve linked property names into node ids
          klass = db.getclass(properties[key].classname)
          issue[key] = klass.lookup(issue[key])
      elif isinstance(properties[key], hyperdb.Multilink):
        # resolve multilinked property names into node ids
        klass = db.getclass(properties[key].classname)
        links = []
        for item in issue[key].split(','):
          if not item.isdigit():
            nodeid = klass.lookup(item)
          else:
            nodeid = item

          if not nodeid in links:
            links.append(nodeid)

        issue[key] = links

  # next we create the issue itself
  issue_id = db.issue.create(**issue)

  print "Created issue%s"%issue_id


if __name__ == '__main__':
  if len(sys.argv) > 1:
    sys.argv.pop(0)
    fname = sys.argv.pop(0)
    failureInfo = file(fname,'r').read()
    if not failureInfo:
      sys.exit(0)
    # convert script arguments to issue properties
    issue = {}
    tracker_home = None
    user = None
    issue['messages']=failureInfo
    for arg in sys.argv:
      arg = arg.split('=')

      if len(arg) == 2:
        if arg[0].lower() == 'tracker_home':
          # get tracker home from argument list
          tracker_home = arg[1]
        elif arg[0].lower() == 'user':
          # get user from argument list
          user = arg[1]
        else:
          # store argument in issue dictionary
          issue[arg[0]] = '='.join(arg[1:])

    # if no tracker_home on argument list, then try to get it from
    # the environment list
    if not tracker_home:
      tracker_home = os.getenv('TRACKER_HOME')

    # if no user on argument list, then try to get it from
    # the environment list
    if not user:
      user = os.getenv('USER')

      if not user:
        user = 'anonymous'

    # chevck if we have a valid tracker home
    if not os.path.exists(os.path.normpath('%s/config.ini'%tracker_home)):
      sys.stderr.write('ERROR: TRACKER_HOME (%s) is not a valid path\n'%tracker_home)
      sys.exit(-1)

    # open a tracker instance
    tracker = instance.open(tracker_home)
    # open the database
    db = tracker.open(user)

    try:
      # create issue
      create(db, **issue)

      # write changes to database
      db.commit()
    finally:
      # make sure we close the database in all conditions
      db.close()
