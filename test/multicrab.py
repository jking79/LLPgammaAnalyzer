#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
from httplib import HTTPException
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand


def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: %prog -c CMD -d DIR [-o OPT]\nThe multicrab command'
                   ' executes "crab CMD OPT" for each task contained in DIR\nUse'
                   ' multicrab -h for help"')

    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--crabCmd", dest="crabCmd",
         help=("The crab command you want to execute for each task in "
               "the DIR"), metavar="CMD")
    parser.add_option("-d", "--projDir", dest="projDir",
         help="The directory where the tasks are located", metavar="DIR")
    parser.add_option("-o", "--crabCmdOptions", dest="crabCmdOptions",
         help=("The options you want to pass to the crab command CMD"
               "tasklistFile"), metavar="OPT", default="")

    (options, args) = parser.parse_args()

    if args:
        parser.error("Found positional argument(s) %s." % args)
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided")
    if not options.projDir:
        parser.error("(-d DIR, --projDir=DIR) option not provided")
    if not os.path.isdir(options.projDir):
        parser.error("Directory %s does not exist" % options.projDir)

    return options


def main():
    """
    Main
    """
    options = getOptions()

    # Execute the command with its arguments for each task.
    for task in os.listdir(options.projDir):
        task = os.path.join(options.projDir, task)
        if not os.path.isdir(task):
            continue
        if str(options.crabCmd) != 'resubmit':
            try :
                print ("Executing (the equivalent of): crab %s %s %s" %
                      (options.crabCmd, task, options.crabCmdOptions))
                crabCommand(options.crabCmd, task, *options.crabCmdOptions.split())
		print (">> End Report >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            except HTTPException, hte :
                    print 'Command not executed'
        else:
            print ("Executing (the equivalent of): crab %s %s %s" %
                  (options.crabCmd, task, options.crabCmdOptions))
            os.system('crab resubmit ' + task) 
if __name__ == '__main__':
    main()   
