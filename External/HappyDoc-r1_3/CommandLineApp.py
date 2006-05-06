#!/usr/bin/env python
#
# Time-stamp: <00/09/04 07:56:44 dhellmann>
#
# COPYRIGHT
#
#   Permission to use, copy, modify, and distribute this software and
#   its documentation for any purpose and without fee is hereby
#   granted, provided that the above copyright notice appear in all
#   copies and that both that copyright notice and this permission
#   notice appear in supporting documentation, and that the name of Doug
#   Hellmann not be used in advertising or publicity pertaining to
#   distribution of the software without specific, written prior
#   permission.
# 
# DISCLAIMER
#
#   DOUG HELLMANN DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
#   INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN
#   NO EVENT SHALL DOUG HELLMANN BE LIABLE FOR ANY SPECIAL, INDIRECT OR
#   CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
#   OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
#   NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
#   CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 


"""Base class for building command line applications.

  The CommandLineApp class makes creating command line applications as
  simple as defining callbacks to handle options when they appear in
  'sys.argv'.

  To do

    - enhance intelligence of option handling

        - boolean options should not need to be implemented as functions

        - enable/disable with/without options

        - type checking for option arguments

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: CommandLineApp.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'Open Source',
    'created':'Tue, 23-May-2000 07:11:43 EDT',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.1 $',
    'date':'$Date: 2001/03/21 20:59:26 $',
    }

#
# Import system modules
#
import getopt
import sys
import string


#
# Import Local modules
#


#
# Module
#

class CommandLineApp:
    """Base class for building command line applications.
    
    Define a __doc__ string for the class to explain what the
    program does.

    When the argumentsDescription field is not empty,
    it will be printed appropriately in the help message.
    
    When the examplesDescription field is not empty,
    it will be printed last in the help message when
    the user asks for help.
    """

    argumentsDescription = ''
    
    examplesDescription = ''

    #
    # Exception names
    #
    ReservedOptionName = 'Reserved option name'
    HelpRequested='Help requested'
    InvalidOptionValue='Invalid value for option'
    InvalidArgument='Invalid argument to program'

    #
    # Globally useful configuration stuff.
    #
    optionHandlerPrefix = 'optionHandler_'
    optTypeLong = 'optTypeLong'
    optTypeShort = 'optTypeShort'
    optTakesParam = 'optTakesParam'
    optNoParam = 'optNoParam'

    verboseLevel = 1

    def __init__(self, commandLineOptions=sys.argv[1:]):
        "Initialize CommandLineApp."
        self.commandLineOptions = commandLineOptions
        self.__learnValidOpts__()
        self.__init_getopt__()
        self.appInit()

    def appInit(self):
        """Override this method to perform application initialization.

        This hook may be easier to override than the __init__ method,
        since it takes no parameters.
        """
        pass

    def shortOptionsStringGet(self):
        """Given the information learned through self.__learnValidOpts__,
        construct a string to be passed to getopt to set the valid
        single character option syntax.
        """
        sas = ''
        for (optName, optTakesValue,
             optLongOrShort, ignore, ignore) in self.shortOptions:
            sas = sas + optName
            if optTakesValue == self.optTakesParam:
                sas = sas + ':'
        return sas

    def longOptionsListGet(self):
        """Given the information learned through self.__learnValidOpts__,
        construct a list to be passed to getopt to set the valid
        long form option syntax.
        """
        lal = []
        for (optName, optTakesValue,
             optLongOrShort, ignore, ignore) in self.longOptions:
            if optTakesValue == self.optTakesParam:
                opt = '%s=' % optName
            else:
                opt = optName
            lal.append(opt)
        return lal

    def __init_getopt__(self):
        "Parse the command line options."
        shortOptionsString = self.shortOptionsStringGet() + 'h'
        longOptionList = self.longOptionsListGet()
        longOptionList.append('help')
        try:
            self.parsedOptions, self.remainingOpts = getopt.getopt(
                self.commandLineOptions, 
                shortOptionsString,
                longOptionList)
        except getopt.error, message:
            self.showHelp(message)
            raise getopt.error, message

    def constructOptionInfo(self, methodName, methodRef):
        """Return an info tuple for an option handler method.

        Given a method name, return the information tuple
        for that option to the program.  The tuple contains:

          (option name,
          flag showing whether the option takes a value,
          flag indicating long or short form option,
          help string for option)
        """
        optionName = methodName[len(self.optionHandlerPrefix):]
        if len(methodRef.func_code.co_varnames) > 1:
            optionTakesValue = self.optTakesParam
        else:
            optionTakesValue = self.optNoParam
        if len(optionName) > 1:
            optionLongForm = self.optTypeLong
        else:
            optionLongForm = self.optTypeShort
        return (optionName, optionTakesValue,
                optionLongForm, methodName, methodRef)

    def methodNameFromOption(self, option):
        """Given the name of an option, construct 
        and return the name of the method to handle it.
        """
        methodName = '%s%s' % (self.optionHandlerPrefix, option)
        return methodName
        
    def scanClassForOptions(self, cRef):
        "Scan through the inheritence hierarchy to find option handlers."
        for parentClass in cRef.__bases__:
            self.scanClassForOptions(parentClass)
        for componentName in cRef.__dict__.keys():
            component = cRef.__dict__[componentName]
            if componentName[:len(self.optionHandlerPrefix)] == self.optionHandlerPrefix and \
               type(component).__name__ == 'function':
                optionInfo = self.constructOptionInfo(componentName, component)
                if optionInfo[0] == 'h':
                    raise CommandLineApp.ReservedOptionName, 'h'
                self.allOptions[ optionInfo[0] ] = optionInfo
                self.allMethods[ componentName ] = optionInfo

    def __learnValidOpts__(self):
        """Derive the options which are valid for this application.

        Examine the methods defined for this class to
        learn what options the developer wants to use.  Options
        can be added by defining an optionHandler method with a
        name like optionHandler_<option name>.  If the option
        handler method takes an argument, the option will require
        an argument as well.
        """
        self.shortOptions = []
        self.longOptions = []
        self.allOptions = {}
        self.allMethods = {}
        self.scanClassForOptions(self.__class__)
        for optionName in self.allOptions.keys():
            optionInfo = self.allOptions[optionName]
            if optionInfo[2] == self.optTypeShort:
                self.shortOptions.append( optionInfo )
            else:
                self.longOptions.append( optionInfo )

    def handleHelpOption(self):
        #
        # Look for -h in self.optList.
        # if found, call help function 
        # then raise HelpRequested
        #
        for option in self.parsedOptions:
            if option[0] == '-h':
                self.showHelp()
                raise CommandLineApp.HelpRequested, 'Help message was printed.'
            if option[0] == '--help':
                self.showVerboseHelp()
                raise CommandLineApp.HelpRequested, 'Help message was printed.'
        
    def main(self, args):
        """Main body of your application.

        This is the main portion of the app, and is run after all of
        the arguments are processed.  Override this method to implment
        the primary processing section of your application.
        """
        pass

    def infoForOption(self, option):
        "Get the stored information about an option."
        optionBase = option
        if optionBase[:2] == '--':
            optionBase = option[2:]
        elif optionBase[0] == '-':
            optionBase = optionBase[1:]
        return self.allOptions[ optionBase ]
        
    def methodReferenceFromName(self, methodName):
        "Return a reference to the method with the given name."
        return self.allMethods[methodName][4]
    
    def run(self):
        """Entry point.

        Process options and execute callback functions as needed.
        This method should not need to be overridden, if the main()
        method is defined.
        """
        try:
            self.handleHelpOption()
            for option, optValue in self.parsedOptions:
                cleanOption = option
                while cleanOption[0] == '-':
                    cleanOption = cleanOption[1:]
                methName = self.methodNameFromOption(cleanOption)
                method = self.methodReferenceFromName(methName)
                optInfo = self.infoForOption(option)
                if optInfo[1] == self.optTakesParam:
                    optArg = string.split(optValue, ',')
                    if len(optArg) <= 1:
                        optArg = optValue
                    method(self, optArg)
                else:
                    method(self)
            apply(self.main, tuple(self.remainingOpts))
        except KeyboardInterrupt:
            try:
                self.interruptHandler()
            except AttributeError:
                sys.stderr.write('Cancelled by user.\n')
                pass
    
    def getSimpleSyntaxHelpString(self):    
        """Return syntax statement.
        
        Return a simplified form of help including only the 
        syntax of the command.
        """
        #
        # Initialize the return value
        #
        helpMsg = ''
        #
        # Initialize some lists of options to organize the
        # printing.
        #
        shortOptionNoArgs = 'h'
        shortOptionArgs = []
        longOptionNoArgs = ['help']
        longOptionArgs = []
        allOptionNames = self.allOptions.keys()
        allOptionNames.sort()
        #
        # Sort out each option into the correct group.
        #
        for option in allOptionNames:
            optName, optTakesValue, optLongOrShort, ignore, ignore = \
                 self.infoForOption(option)
            methodName = self.methodNameFromOption(option)
            method = self.methodReferenceFromName(methodName)
            if optTakesValue == self.optTakesParam:
                if optLongOrShort == self.optTypeLong:
                    longOptionArgs.append(optName)
                else:
                    shortOptionArgs.append(optName)
            else:
                if optLongOrShort == self.optTypeLong:
                    longOptionNoArgs.append(optName)
                else:
                    shortOptionNoArgs = shortOptionNoArgs + optName
        #
        # Print the name of the command, and the short options
        # which take no arguments.
        #
        helpMsg = '%s%s [-%s] ' % (helpMsg, sys.argv[0], shortOptionNoArgs)
        #
        # Print the short options which take arguments.
        #
        for option in shortOptionArgs:
            methodName = self.methodNameFromOption(option)
            method = self.methodReferenceFromName(methodName)
            valueName = method.func_code.co_varnames[1]
            helpMsg = '%s\n\t\t[-%s %s] ' % (helpMsg, option, valueName)
        #
        # Print the long options which take no arguments.
        #
        for option in longOptionNoArgs:
            helpMsg = '%s\n\t\t[--%s] ' % (helpMsg, option)
        #
        # Print the long options which take arguments.
        #
        for option in longOptionArgs:
            methodName = self.methodNameFromOption(option)
            method = self.methodReferenceFromName(methodName)
            valueName = method.func_code.co_varnames[1]
            helpMsg = '%s\n\t\t[--%s=%s] ' % (helpMsg, option, valueName)
        #
        helpMsg = '%s\n\n' % helpMsg
        return helpMsg

    def showSimpleSyntaxHelp(self):
        "Show basic syntax message."
        txt = self.getSimpleSyntaxHelpString()
        print txt

    def getOptionHelpString(self, dashInsert, option, 
                            valueInsert, docString):
        "Build the help string for an option."
        baseString = '\t%s%s%s\n\n' % (dashInsert, 
                                       option, 
                                       valueInsert, 
                                       )
        docStringLines = string.split(docString, '\n')
        optionString = baseString
        for line in docStringLines:
            line = string.strip(line)
            if not line: continue 
            optionString = '%s\t\t\t%s\n' % (optionString, line)
        optionString = optionString + '\n\n'
        return optionString
        
    def showOptionHelp(self, dashInsert, option, valueInsert, docString):
        'Format and print the help message for a single option.'
        #
        #print '\t%s%s%s' % (dashInsert, option, valueInsert)
        #print ''
        #print '\t\t\t%s' % docString
        #print ''
        print self.getOptionHelpString(dashInsert, option,
                                       valueInsert, docString)

    def showVerboseSyntaxHelp(self):
        "Show a full description of all options and arguments."
        txt = self.getVerboseSyntaxHelpString()
        print txt
        

    def getVerboseSyntaxHelpString(self):
        """Return the full description of the options and arguments.

        Show a full description of the options and arguments to the
        command in something like UNIX man page format. This includes
        
          - a description of each option and argument, taken from the
                __doc__ string for the optionHandler method for
                the option
                
          - a description of what additional arguments will be processed,
                taken from the class member argumentsDescription
                
        """
        #
        # Show them how to use it.
        #
        helpMsg = '\nSYNTAX:\n\n\t'
        #
        # Show the options.
        #
        helpMsg = '%s%sOPTIONS:\n\n' % (helpMsg, 
                                        self.getSimpleSyntaxHelpString())
        helpMsg = string.join( [helpMsg,
        
                    self.getOptionHelpString('-', 'h', '', 
                    'Displays abbreviated help message.'),
                    
                    self.getOptionHelpString('--', 'help', '', 
                    'Displays complete usage information.'),
                    ],
                    ''
                )
        #
        allOptionNames = self.allOptions.keys()
        allOptionNames.sort()
        for option in allOptionNames:
            optName, optTakesValue, optLongOrShort, ignore, ignore = \
                 self.infoForOption(option)
            methodName = self.methodNameFromOption(option)
            method = self.methodReferenceFromName(methodName)
            docString = method.__doc__
            if optTakesValue == self.optTakesParam:
                valueInsert = method.func_code.co_varnames[1]
            else:
                valueInsert = ''
            if optLongOrShort == self.optTypeLong:
                dashInsert = '--'
                if valueInsert:
                    valueInsert = '=%s' % valueInsert
            else:
                dashInsert = '-'
                if valueInsert:
                    valueInsert = ' %s' % valueInsert
            helpMsg = string.join( [ helpMsg,
                                     self.getOptionHelpString(dashInsert,
                                                              option,
                                                              valueInsert,
                                                              docString) ],
                                    '')
            #self.showOptionHelp(dashInsert, option, valueInsert, docString)
        helpMsg = helpMsg + '\n'
        if self.argumentsDescription:
            helpMsg = '%sARGUMENTS:\n\n%s\n\n' % (helpMsg,
                                                  self.argumentsDescription)
            #print 'ARGUMENTS:'
            #print ''
            #print self.argumentsDescription
            #print ''
        return helpMsg

    def showVerboseHelp(self):
        """Show a verbose help message explaining how to use the program.

        This includes:
        
           * a verbose description of the program, taken from the __doc__
             string for the class
             
           * an explanation of each option, produced by
             showVerboseSyntaxHelp()
             
           * examples of how to use the program for specific tasks,
             taken from the class member examplesDescription
             
        """
        #
        # Show the program name and
        # a description of what it is for.
        #
        print ''
        print '%s\n' % sys.argv[0]
        if self.__class__.__doc__:
            print ''
            for line in string.split(self.__class__.__doc__, '\n'):
                line = string.strip(line)
                if not line: continue
                print '\t%s' % line
            print ''
        self.showVerboseSyntaxHelp()
        if self.examplesDescription:
            print 'EXAMPLES:'
            print ''
            print self.examplesDescription
            print ''

    def showHelp(self, errorMessage=None):
        "Display help message when error occurs."
        #
        # If they made a syntax mistake, just
        # show them how to use the program.  Otherwise,
        # show the full help message.
        #
        if errorMessage:
            print ''
            print 'ERROR: ', errorMessage
            print ''
            print ''
            print '%s\n' % sys.argv[0]
            print ''
        self.showSimpleSyntaxHelp()

    def statusMessage(self, msg='', verboseLevel=1, error=None):
        """Print a status message to output.
        
        Arguments
        
            msg=''            -- The status message string to be printed.
            
            verboseLevel=1    -- The verbose level to use.  The message
                              will only be printed if the current verbose
                              level is >= this number.
                              
            error=None        -- If true, the message is considered an error and
                              printed as such.
                              
        """
        if self.verboseLevel >= verboseLevel:
            if error:
                output = sys.stderr
            else:
                output = sys.stdout
            output.write('%s\n' % msg)
            output.flush()
        return

    def debugMethodEntry(self, methodName, debugLevel=3, **nargs):
        "Print a method entry debug statement with named arguments."
        arg_str = string.join( map( lambda pair: '%s=%s' % pair, nargs.items() ), ', ')
        return self.statusMessage('%s(%s)' % (methodName, arg_str), debugLevel)
    
    def errorMessage(self, msg=''):
        'Print a message as an error.'
        self.statusMessage('ERROR: %s\n' % msg, 0)
                
    def optionHandler_v(self):
        """Increment the verbose level.
        Higher levels are more verbose.
        The default is 1.
        """
        self.verboseLevel = self.verboseLevel + 1
        self.statusMessage('New verbose level is %d' % self.verboseLevel,
                           4)
        return

    def optionHandler_q(self):
        'Turn on quiet mode.'
        self.verboseLevel = 0

        
class TestApp(CommandLineApp):
    """    This is a simple test application.

    It defines several optionHandler methods to handle
    some example options.  One option of each type is
    handled by an example.

    The __doc__ string for the class should contain
    the info about how to use the application. """

    examplesDescription = \
"""    a description of how to use the program 
    in various ways goes here.
"""                

    argumentsDescription = \
"""    a description of other arguments goes here
"""

    defaultLongFormOption='<default value>'

    def optionHandler_a(self, optValue):
        'get a value for a'
        print '%s: handling a: %s' % (self.__class__.__name__, optValue)

    def optionHandler_b(self):
        'toggle the value of b'
        print '%s: handling b' % (self.__class__.__name__,)

    def optionHandler_long_form_option(self):
        'boolean option'
        print '%s: handling long-form-option' % (self.__class__.__name__,)

    def optionHandler_long_form_with_value(self, optValue):
        """First line of help.
            get a value for long form option

            Default:<manually inserted>"""
        print '%s: handling long-form-with-value: %s' % (self.__class__.__name__,
                                 optValue)

    def main(self):
        'main loop'
        print '%s: LEFT OVERS: ' % (self.__class__.__name__,), self.remainingOpts

        
class SubClassTestApp(TestApp):
    'new doc string'
    
    def optionHandler_a(self, newA):
        'Doc string for SubClassTestApp'
        print 'New A:', newA
        TestApp.optionHandler_a(self, newA)
    
    def optionHandler_z(self):
        'Doc string for SubClassTestApp'
        print '%s -z' % self.__class__.__name__
        
    def optionHandler_option_list(self, optionList):
        'Doc string for SubClassTestApp'
        if type(optionList) == type([]):
            print '%s -z list: ' % self.__class__.__name__, optionList
        else:
            print '%s -z string: %s' % (self.__class__.__name__, optionList)
        
    def main(self, args):
        TestApp.main(self, args)
        raise 'not an error'

if __name__ == '__main__':

    testArgs = [
        #['-aaval', '-z'],
        ['-h'],
        ['--help'],
        #[ '--option_list=a,b,c', '--option_list', 'a,b,c', '--option_list=a']
        ]
    for t in testArgs:
        print '\n\n'
        print 'RUNNING: ', t
        try:
            myApp = SubClassTestApp(t)
            myApp.run()
        except CommandLineApp.HelpRequested:
            pass
        except getopt.error, msg:
            pass

    print 'Press <return>...'
    blah = sys.stdin.readline()
