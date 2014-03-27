#!/usr/bin/env python

# TODO: tab completion

import sys, os, fcntl, time, readline
from subprocess import Popen, PIPE
from math import log10

HIST_FILE = os.path.join(os.getenv("HOME"), ".sewrap.hist")
# Right now, unlimited.  Set with readline.set_history_length( <lines> )

def sewrapper(args, pipeout=None):
    flags = [a for a in args if a.startswith('-')]
    filenames = [a for a in args if not a.startswith('-')]
    commands = ['', 'pid ::= %i'%os.getpid()]
    title = Title('SE', 'P')
    if len(filenames) > 1:
        print "Warning: Received more than one filename:"
        print "\t", filenames
        print "\tUsing only the first."
    if filenames:
        commands.append('load "%s"'%filenames[0])
        title(title=filenames[0])

    p = Popen('evolver ' + ' '.join(flags), shell=True,
                stdin=PIPE, stdout=PIPE, close_fds=True)
    p.stdin.flush()
    fcntl.fcntl(p.stdout, fcntl.F_SETFL, os.O_NONBLOCK)
    try:
        readline.read_history_file(HIST_FILE)
    except IOError:
        pass

    sync = False

    while 1:
        try:
            try:
                line = p.stdout.read()
            except IOError:
                if not sync:
                    time.sleep(0.1)
                    continue

            sync = False
            if line == '':
    #            print "EOF"
                break
            lines = line.split('\n')
            for l in lines:
                if l.startswith("Difference:"):
                    try:
                        diff = log10(float(l.split()[1]))
                    except ValueError:
                        sys.stdout.write("SEWRAP: Can't figure out difference from this line.")
                    else:
                        title(running='%.1f'%diff)
            if len(lines) > 1:
                sys.stdout.write('\n'.join(lines[:-1]) + '\n')
                line = lines[-1]
            if line:
                title(running='Stop')
                if pipeout is not None:
                    pipeout.write("Ready\n")
                    pipeout.flush()
                if commands:
                    linein = commands[0]
                    commands = commands[1:]
                    sys.stdout.write(line + linein + '\n')
                    writePID = True
                else:
                    try:
                        linein = handle_input(line, title)
                    except EOFError:
                        sys.stdout.write('\n')
                        break # p.stdin.close() should send EOF
                p.stdin.write(linein + '\n')
                p.stdin.flush()
                title(running='Run')
                if pipeout is not None:
                    pipeout.write("Running\n")
                    pipeout.flush()
        except KeyboardInterrupt:
            pass
#            # Note that ^C also goes through to child automagically.  But in case
#            # evolver is really stuck, we'll give an option to force a quit
#            response = ''
#            while response not in ('y', 'n', 'q'):
#                response = raw_input(' detected.  Do you wish to exit? (y/n): ')[0]
#            if response in ('y', 'q'):
#                break
        except SyncInterrupt:
            sync = True

    p.stdout.close()
    p.stdin.close()
    try:
        readline.write_history_file(HIST_FILE)
    except IOError:
        print "Warning: Could not save history file to", HIST_FILE
    title(running='Quit')
    return 0

def handle_input(line, title):
    while True:
        linein = raw_input(line)
        if linein.startswith('#'):
            cmds = linein[1:].split()
            if cmds[0] == 'title':
                title(title=' '.join(cmds[1:]))
            elif cmds[0] == 'sync':
                raise SyncInterrupt # Kicks out back to the reading pipe loop
            else:
                sys.stdout.write('Did not recongnize command: #%s\n'%cmds[0])
        else:
            return linein

class Title():
    """Keep track of title info."""
    def __init__(self, title="", running=""):
        self.title = title
        self.running = running
        self.in_screen = os.environ.get('STY')

    def __call__(self, title=None, running=None):
        if title is not None:
            self.title = title
        if running is not None:
            self.running = running

        titlestring = "[%s]%s"%(self.running, self.title)
        if self.in_screen:
            sys.stdout.write("%sk%s%s\\"%(chr(27), titlestring, chr(27)))
        else:
            sys.stdout.write("%s]0;%s%s"%(chr(27), titlestring, chr(7)))

class SyncInterrupt(Exception):
    pass

if __name__ == "__main__":
    sys.exit(sewrapper(sys.argv[1:]))
