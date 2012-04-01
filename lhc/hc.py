#!/usr/bin/env python
# :exec set tabstop=4 softtab expandtab encoding=utf8 :

"""
Copyright (c) 2009, Don Peterson
Copyright (c) 2011, Vernon Mauery
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.
* The names of the contributors may not be used to endorse or
promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

#----------------------------------
# Python library stuff
from __future__ import division
import sys, getopt, os, time, readline
from socket import htonl
from atexit import register as atexit
from string import strip
import traceback
import re as regex
from tempfile import mkstemp
from debug import *
import config

try: from pdb import xx  # pdb.set_trace is xx; easy to find for debugging
except: pass

#----------------------------------
# Modules we are dependent on
try:
    import mpmath as m
    from simpleparse.stt import TextTools
    from simpleparse import generator
except ImportError:
    print """
This is a complex program that requires several external python
libraries.  Please install mpmath, simpleparse

apt-get install python-mpmath python-simpleparse

"""
    sys.exit(1)

#----------------------------------
# Modules needed in our package
from numeric import *
from stack import Stack
from mpformat import mpFormat
import constants
import console

# You may create your own display (GUI, curses, etc.) by derivation.  The
# default Display object just prints to stdout and should work with any
# console.
from display import Display

out = sys.stdout.write
err = sys.stderr.write
nl = "\n"
# Status numbers that can be returned
status_ok               = 0
status_quit             = 1
status_error            = 2
status_unknown_command  = 3
status_ok_no_display    = 4
status_interrupted      = 5
JULIAN_UNIX_EPOCH = Julian("1Jan1970:00:00:00")


class ParseError(Exception):
    pass

def nop(*args):
    """
    unimplimented
    """
    return None


class Calculator(object):
    def __init__(self, arguments, options):
        self.errors = []
        self.stack = Stack()
        self.stack_index = True
        self.constants = constants.ParseRawData()
        self.display = Display()     # Used to display messages to user
        self.fp = mpFormat()         # For formatting floating point numbers
        self.ap = mpFormat()         # For formatting arguments of complex numbers
        self.number = Number()
        self.registers = {}          # Keeps all stored registers
        self.commands_dict = {
            # Values are
            # [
            #   implementation function to call,
            #   number of stack arguments consumed,
            #   Optional dictionary of other needed things, such as:
            #       "pre" : [func, (args,)]  Function to execute before calling
            #                                the implementation function.
            #       "post": [func, (args,)]  Function to execute after calling
            #                                the implementation function.
            # ]

            # Binary functions
            "+"        : [self.add, 2],
            "-"        : [self.subtract, 2],
            "*"        : [self.multiply, 2],
            "/"        : [self.divide, 2],
            "div"      : [self.integer_divide, 2],
            "%"        : [self.Mod, 2],
            "mod"      : [self.Mod, 2],
            "and"      : [self.bit_and, 2],
            "&"        : [self.bit_and, 2],
            "or"       : [self.bit_or, 2],
            "|"        : [self.bit_or, 2],
            "xor"      : [self.bit_xor, 2],
            "<<"       : [self.bit_leftshift, 2],
            ">>"       : [self.bit_rightshift, 2],
            "%ch"      : [self.percent_change, 2],
            "comb"     : [self.combination, 2],  # Combinations of y taken x at a time
            "perm"     : [self.permutation, 2],  # Permutations of y taken x at a time
            "pow"      : [self.power, 2],  # Raise y to the power of x
            "^"        : [self.power, 2],  # Raise y to the power of x
            "atan2"    : [self.atan2, 2], # {"post" : self.Conv2Deg}], #
            "hypot"    : [self.hypot, 2],  # sqrt(x*x + y*y)
            "round"    : [self.Round, 2],  # Round y to nearest x
            "in"       : [self.In, 2],     # True if x is in interval y
            "=="       : [self.Equal, 2],     # True if x == y
            "!="       : [self.NotEqual, 2],  # True if x != y
            "<"        : [self.LessThan, 2],  # True if x < y
            "<="       : [self.LessThanEqual, 2], # True if x <= y
            ">"        : [self.GreaterThan, 2],      # True if x > y
            ">="       : [self.GreaterThanEqual, 2], # True if x >= y
            "="        : [self.DisplayEqual, 2],  # True if displayed strings of x & y are equal
            "iv"       : [self.ToIV, 2],   # Convert to [y,x] interval number
            "gcf"      : [self.gcf, 2],  # find the greatest common factor
            "lcd"      : [self.lcd, 2],  # find the lowest common denominator

            # Unary functions
            "I"        : [self.Cast_i, 1],  # Convert to integer
            "Q"        : [self.Cast_q, 1],  # Convert to rational at display resolution
            "QQ"       : [self.Cast_qq, 1], # Convert to rational at full precision
            "R"        : [self.Cast_r, 1],  # Convert to real number
            "C"        : [self.Cast_c, 1],  # Convert to complex number
            "T"        : [self.Cast_t, 1],  # Convert to time/date
            "V"        : [self.Cast_v, 1],  # Convert to interval number
            "cast"     : [self.cast, 1],  # Convert integer to current C int type
            "IP"       : [self.IP, 1],  # Convert to ip address
            "2deg"     : [self.ToDegrees, 1],  # Convert x to radians
            "2rad"     : [self.ToRadians, 1],  # Convert x to degrees
            "unix"     : [self.ToUnix, 1],  # Convert julian to unix timestamp
            "julian"   : [self.ToJulian, 1], # Convert unix timestamp to julian
            "2hr"      : [self.hr, 1],    # Convert to decimal hour format
            "2hms"     : [self.hms, 1],   # Convert to hour/minute/second format
            "fp"       : [self.first_part, 1],    # Integer part of x
            "sp"       : [self.second_part, 1],    # Fractional part of x

            "inv"      : [self.reciprocal, 1], # reciprocal of x
            "~"        : [self.bit_negate, 1],   # Flip all the bits of x
            "split"    : [self.split, 1], # Take rational, complex, or interval apart
            "chop"     : [self.Chop, 1],  # Convert x to its displayed value
            "conj"     : [self.conj, 1],  # Complex conjugate of x
            "sqrt"     : [self.sqrt, 1],  # Square root of x
            "cbrt"     : [self.cbrt, 1],  # Cube root of x
            "root"     : [self.root, 2],  # nth root of x
            "roots"    : [self.roots, 2],  # nth roots of x
            "sqr"      : [self.square, 1],# Square x
            "neg"      : [self.negate, 1], # negative of x
            "mid"      : [self.mid, 1],   # Take midpoint of interval number
            "sum"      : [self.Sum, 'x'],  # sum of top x values (depth sum for all)
            "!"        : [self.Factorial, 1],  # factorial
            "floor"    : [self.floor, 1], # Largest integer <= x
            "ceil"     : [self.ceil, 1],  # Smallest integer >= x
            "abs"      : [self.abs, 1],   # Absolute value of x
            "arg"      : [self.arg, 1],# {"post" : Conv2Deg}],  # Argument of complex
            "ln"       : [self.ln, 1],    # Natural logarithm
            "log2"     : [self.log2, 1],     # Base 2 logarithm
            "log"      : [self.log10, 1], # Base 10 logarithm
            "exp"      : [self.exp, 1], # Exponential function
            "bits"     : [self.bits, 1], # calculate the number of bits required for this integer
            "db"       : [self.db, 1], # take a number in db and express it as a std ratio
            "bd"       : [self.bd, 1], # take a ratio and express it in db

            # 0-nary functions
            "rand"     : [self.rand, 0],  # Uniform random number
            "ts"       : [self.unix_ts, 0], # return unix timestamp

            # trig functions
            "sin"      : [self.sin, 1],   # {"pre"  : self.Conv2Rad}],
            "cos"      : [self.cos, 1],   # {"pre"  : self.Conv2Rad}],
            "tan"      : [self.tan, 1],   # {"pre"  : self.Conv2Rad}],
            "asin"     : [self.asin, 1],  # {"post" : self.Conv2Deg}],
            "acos"     : [self.acos, 1],  # {"post" : self.Conv2Deg}],
            "atan"     : [self.atan, 1],  # {"post" : self.Conv2Deg}],
            "sec"      : [self.sec, 1],   # {"pre"  : self.Conv2Rad}],
            "csc"      : [self.csc, 1],   # {"pre"  : self.Conv2Rad}],
            "cot"      : [self.cot, 1],   # {"pre"  : self.Conv2Rad}],
            "asec"     : [self.asec, 1],  # {"post" : self.Conv2Deg}],
            "acsc"     : [self.acsc, 1],  # {"post" : self.Conv2Deg}],
            "acot"     : [self.acot, 1],  # {"post" : self.Conv2Deg}],
            "sinh"     : [self.sinh, 1],  # {"pre"  : self.Conv2Rad}],
            "cosh"     : [self.cosh, 1],  # {"pre"  : self.Conv2Rad}],
            "tanh"     : [self.tanh, 1],  # {"pre"  : self.Conv2Rad}],
            "asinh"    : [self.asinh, 1], # {"post" : self.Conv2Deg}],
            "acosh"    : [self.acosh, 1], # {"post" : self.Conv2Deg}],
            "atanh"    : [self.atanh, 1], # {"post" : self.Conv2Deg}],
            "sech"     : [self.sech, 1],  # {"pre"  : self.Conv2Rad}],
            "csch"     : [self.csch, 1],  # {"pre"  : self.Conv2Rad}],
            "coth"     : [self.coth, 1],  # {"pre"  : self.Conv2Rad}],
            "asech"    : [self.asech, 1], # {"post" : self.Conv2Deg}],
            "acsch"    : [self.acsch, 1], # {"post" : self.Conv2Deg}],
            "acoth"    : [self.acoth, 1], # {"post" : self.Conv2Deg}],

            # Stack functions
            "clr"      : [self.ClearStack, 0],
            "clear"    : [self.Reset, 0], # Reset the calculator state
            "stack"    : [self.SetStackDisplay, 1],
            "lastx"    : [self.lastx, 0], # Recall last x used
            "swap"     : [self.swap, 0],   # swap x and y
            "roll"     : [self.roll, 0],  # Roll stack
            "rolld"    : [self.rolld, 0],  # Roll stack down
            "over"     : [self.over, 0],  # push y onto the stack at the top
            "pick"     : [self.pick, 1],  # pick stack[x] off the stack and push it at the top
            "drop"     : [self.drop, 1],   # Pop x off the stack
            "drop2"    : [self.drop2, 2],   # Pop x and y off the stack
            "dropn"    : [self.dropn, 'x'],   # Pop x items off the stack
            "dup"      : [self.dup, 1],   # Push a copy of x onto the stack
            "dup2"     : [self.dup2, 2],   # Push a copy of x and y onto the stack
            "dupn"     : [self.dupn, 'x'],  # duplicate top x values on stack
            "depth"    : [self.depth, 0],  # Push stack depth onto stack

            # constants
            "phi"      : [self.Phi, 0],   # Golden ratio
            "pi"       : [self.Pi, 0],
            "e"        : [self.E, 0],
            "i"        : [self.I, 0],
            "j"        : [self.I, 0],
            "const"    : [self.const, 0],  # grab a list of constants

            # network functions
            # same net - 3 args -- 2 ips and a netmask
            # broadcast - 2 args - ip and a netmask
            "netmask"  : [self.netmask, 2],  # apply a given netmask to an ip address
            "cidr"     : [self.cidr, 2],  # apply a given netmask to an ip address
            # net match
            "le"       : [nop, 0],  # set little-endian integer mode
            "be"       : [nop, 0],  # set big-endian integer mode
            "htonl"    : [self.ntohl, 1],  # return htonl x
            "ntohl"    : [self.ntohl, 1],  # return ntohl x
            "=net"     : [self.samenet, 2],  # check to see if y and x are on the same subnet

            # Other stuff
            self.help.__name__: [self.help, 'match',
                            {
                                'regex': regex.compile(r"(?:help|[?])\s+(.+)?"),
                                'grammar': "(('?' / 'help'),(ws,delimited_func)?)",
                            }, ], # help
            "warranty" : [self.warranty, 0],  # Show warranty
            "quit"     : [self.quit, 0],  # Exit the program
            "deg"      : [self.deg, 0],  # Set degrees for angle mode
            "rad"      : [self.rad, 0],  # Set radians for angle mode
            "regs"     : [self.PrintRegisters, 0],
            "cfg"      : [self.ShowConfig, 0], # Show configuration
            "modulo"   : [self.Modulus, 1], # All answers displayed with this modulus
            "clrg"     : [self.ClearRegisters, 0],
            ">>."      : [self.display.logoff, 0],  # Turn off logging

            # Display functions
            "mixed"    : [self.mixed, 1], # Toggle mixed fraction display
            "rat"      : [self.Rationals, 1], # Toggle whether to use rationals
            "down"     : [self.ToggleDowncasting, 1],
            "on"       : [self.display.on, 0],  # Turn display of answers on
            "off"      : [self.display.off, 0],  # Turn display of answers off
            "prec"     : [self.Prec, 1],  # Set calculation precision
            "digits"   : [self.digits, 1],# Set significant figures for display
            "width"    : [self.width, 1], # Set line width
            "comma"    : [self.comma, 1], # Toggle comma decorating
            "fix"      : [self.fix, 0],  # Fixed number of places after decimal point
            "sig"      : [self.sig, 0],  # Display signification figures
            "sci"      : [self.sci, 0],  # Scientific notation display
            "eng"      : [self.eng, 0],  # Engineering display
            "engsi"    : [self.engsi, 0],  # Engineering display with SI prefix
            "raw"      : [self.raw, 0],  # raw fp mode
            "brief"    : [self.brief, 1],  # Fit number on one line
            "iva"      : [self.iva, 0],  # Interval display
            "ivb"      : [self.ivb, 0],  # Interval display
            "ivc"      : [self.ivc, 0],  # Interval display
            "show"     : [self.Show, 0],  # Show full precision of x register
            "debug"    : [self.Debug, 1], # Toggle the debug variable
            # angle modes
            "polar"    : [self.Polar, 0],  # Complex number display
            "rect"     : [self.Rectangular, 0],  # Complex number display
            # integer modes
            "sx"       : [self.C_sX, 1],  # Unsigned n-bit integer mode
            "ux"       : [self.C_uX, 1],  # Signed n-bit integer mode
            self.C_int.__name__: [self.C_int, 'match',
                            {
                                'regex': regex.compile(r"([su])([0-9]+)"),
                                'grammar': "([su],[0-9]+)"
                            }, ], # c int u32 style
            "dec"      : [self.dec, 0],  # Decimal display for integers
            "hex"      : [self.hex, 0],  # Hex display for integers
            "oct"      : [self.oct, 0],  # Octal for integers
            "bin"      : [self.bin, 0],  # Binary display for integers
            "roman"    : [self.roman, 0],  # roman numeral display for integers
            # The none display mode is primarily intended for debugging.  It
            # displays makes the mpmath numbers display in their native formats.

            # Some other math functions
            "gamma"    : [self.gamma, 1],
            "zeta"     : [self.zeta, 1],
            "ncdf"     : [self.Ncdf, 1],
            "invn"     : [self.Incdf, 1],

        }
        self.commands_dict['?'] = self.commands_dict['help']
        #t = datetime.now()
        #M.rand('init', 64)
        #M.rand('seed', (t.year+t.month+t.day)/(t.microsecond+1)+
        #    (((t.hour*60)+t.minute)*60+t.second*1000000)+t.microsecond)
        # set up readline stuff
        # check for dir and file
        # consider adding tab completion
        try:
            os.makedirs(os.path.expanduser('~')+'/.pycalc')
        except OSError:
            pass
        if hasattr(readline, "read_history_file"):
            try:
                readline.read_history_file(os.path.expanduser('~')+'/.pycalc/history')
            except IOError:
                pass
            atexit(self.cleanup)

        defined_functions = ["'nop'"]
        funcs = self.commands_dict.keys()
        funcs.sort(reverse=True)
        grammar_funcs = {}
        for f in funcs:
            funcinf = self.commands_dict[f]
            name = funcinf[0].__name__
            if name not in defined_functions:
                defined_functions.append(name)
            if len(funcinf) == 3 and 'grammar' in funcinf[2]:
                if name not in grammar_funcs:
                    grammar_funcs[name] = [funcinf[2]['grammar']]
                else:
                    grammar_funcs[name].append(funcinf[2]['grammar'])
            else:
                if name not in grammar_funcs:
                    grammar_funcs[name] = ["'%s'" %f,]
                else:
                    grammar_funcs[name].append("'%s'" %f)
        functions_grammar = []
        for name,func_vals in grammar_funcs.iteritems():
            functions_grammar.append('        %s := %s' % (name, ' / '.join(func_vals)))
        grammar_funcs = '\n'.join(functions_grammar)
        # print defined_functions
        grammar = ''.join(["""
        calculator_grammar := statement / ws
        statement := simple_statement / (simple_statement, ws, statement)
        simple_statement := delimited_func / constant / ipaddr / number
        constant := 'const'
        ipaddr := ipv6cidr / ipv4cidr / ipv6 / ipv4
        #ipv6 := (((hex_chars)?),':')+,((hex_chars)?),(':',((hex_chars)?))+
        ipv6cidr := ipv6,'/',[0-9],[0-9]?,[0-9]?
        ipv4cidr := ipv4,'/',[0-9],[0-9]?
        ipv6 := '::1' / '::' / ((hex_chars,':')+,(':'?,hex_chars)+)
        ipv4 := [0-9],[0-9]?,[0-9]?,'.',[0-9],[0-9]?,[0-9]?,'.',[0-9],[0-9]?,[0-9]?,'.',[0-9],[0-9]?,[0-9]?
        number := roman_number / rational_number / scaler_number / compound_number
        roman_number := roman_numeral / roman_numeral,roman_numeral
        roman_numeral := [Mm] / [Dd] / [Cc] / [Ll] / [Xx] / [Vv] / [Ii]
        rational_number := dec_whole , '/' , dec_whole
        compound_number := vector / array
        scaler_number := julian / complex_number / imag_number / real_number
        complex_number := (real_number_ns,('+'/'-'),imag_number) / ('(',real_number,',',real_number,')') / ('(', real_number, (',', ows)?, '<', real_number, ')')
        imag_number := real_number_ns,[ij]
        array := '[', vector_list, ']'
        vector_list := (vector, ',', vector_list) / vector
        vector := '[', real_number_list, ']'
        real_number_list := real_number, (','?, real_number)*
        real_number := ows, real_number_ns, ows
        real_number_ns := bin_number / oct_number / hex_number / dec_number
        # now and today are 'numbers' interpreted by Julian class
        julian := 'now' / 'today' / datetime
        datetime := [0-9],[0-9]?,month,[0-9],[0-9],[0-9],[0-9],[-:],[0-9],[0-9],':',[0-9],[0-9],':',[0-9],[0-9],('.',[0-9]+)?
        month := ([Jj],[Aa],[Nn]) / ([Ff],[Ee],[Bb]) / ([Mm],[Aa],[Rr]) / ([Aa],[Pp],[Rr]) / ([Mm],[Aa],[Yy]) / ([Jj],[Uu],[Nn]) / ([Jj],[Uu],[Ll]) / ([Aa],[Uu],[Gg]) / ([Ss],[Ee],[Pp]) / ([Oo],[Cc],[Tt]) / ([Nn],[Oo],[Vv]) / ([Dd],[Ee],[Cc])
        hex_number := '-'?,'0x',hex_chars
        dec_number := dec_float / dec_whole
        dec_whole := '-'?,dec_chars
        dec_float := '-'?,((dec_chars,'.',dec_chars)/(dec_chars,'.')/('.',dec_chars)/dec_chars),('e','-'?,dec_chars)?
        oct_number := '-'?,'0o',oct_chars
        bin_number := '-'?,'0b',bin_chars
        hex_chars := [0-9A-Fa-f]+
        dec_chars := [0-9]+
        oct_chars := [0-7]+
        bin_chars := [01]+
        ows := [ \n\t]*
        ws := [ \n\t],ows
        delimited_func := (ws,func,ws) / (ws,func) / (func,ws) / func
        """,
            "func := %s\n" % ' / '.join(defined_functions),
            grammar_funcs,
        ])
        # print grammar
        try:
            self.parser = generator.buildParser(grammar).parserbyname('calculator_grammar')
        except:
            print "Parser failed to build.  This may not work at all..."
            type,value,tb = sys.exc_info()
            traceback.print_exception(type, value, tb, None, sys.stdout)
        self.chomppre = regex.compile(r"^\s*")
        self.chomppost = regex.compile(r"\s*$")

        #---------------------------------------------------------------------------
        #---------------------------------------------------------------------------
        # Global variables
        self.stdin_finished = False  # Flags when stdin has reached EOF
        self.comment_line = "\xec\xeb"
        self.eof = "\xed\xee"
        self.argument_types = "%sThe two arguments must be the same type"
        self.factorial_cache = {0:1, 1:1, 2:2}
        self.process_stdin = False   # -s If true, our input comes from stdin
        self.run_checks = False      # -c Run checks
        self.quiet = False           # -q If true, don't print initial message
        self.testing = False         # -t If true, exit with nonzero status if x!=y
        self.tee_is_on = False       # True when logging to a file
        self.use_default_config_only = False

        # Used for binary conversions
        self.hexdigits = {
            "0" : "0000", "1" : "0001", "2" : "0010", "3" : "0011", "4" : "0100",
            "5" : "0101", "6" : "0110", "7" : "0111", "8" : "1000", "9" : "1001",
            "a" : "1010", "b" : "1011", "c" : "1100", "d" : "1101", "e" : "1110",
            "f" : "1111"}

        self.RunChecks()
        config.load()
        self.CheckEnvironment()
        self.GetConfiguration()

        if options.default_config:
            self.display.msg("Using default configuration only")
        if options.version:
            self.display.msg("hc version 7 (29 Mar 2012)")

    #---------------------------------------------------------------------------
    # Utility functions

    def use_modular_arithmetic(self, x, y):
        return (isint(x) and isint(y) and abs(config.cfg["modulus"]) > 1)

    def TypeCheck(self, x, y):
        if (not config.cfg["coerce"]) and (type(x) != type(y)):
            raise ValueError(self.argument_types % fln())

    def DownCast(self, x):
        """
        If x can be converted to an integer with no loss of information,
        do so.  If its a complex that can be converted to a real, do so.
        """
        if config.cfg["downcasting"] == False:
            return x
        if x == inf or x == -inf:
            return x
        elif isinstance(x, Rational):
            if x.d == 1:
                return x.n
        elif isinstance(x, m.mpf):
            if int(x) == x:
                return int(x)
        elif isinstance(x, m.mpc):
            if x.imag == 0:
                x = x.real
                if int(x) == x:
                    return int(x)
        return x

    def Conv2Deg(self, x):
        """
        Routine to convert the top of the stack element to degrees.  This
        is typically done after calling inverse trig functions.
        """
        try:
            if config.cfg["angle_mode"] == "deg":
                if isinstance(x, m.mpc):  # Don't change complex numbers
                    return x
                if isinstance(x, Zn): x = int(x)
                return m.degrees(x)
            return x
        except:
            raise ValueError("%sx can't be converted from radians to degrees" % fln())

    def Conv2Rad(self, x):
        """
        Routine to convert the top of the stack element to radians.  This
        is typically done before calling trig functions.
        """
        try:
            if config.cfg["angle_mode"] == "deg":
                if isinstance(x, m.mpc):  # Don't change complex numbers
                    return x
                if isinstance(x, Zn): x = int(x)
                return m.radians(x)
            return x
        except:
            raise ValueError("%sx can't be converted from degrees to radians" % fln())

    #---------------------------------------------------------------------------
    # Binary functions

    def add(self, y, x):
        """
    Usage: y x +

    Return the sum of the bottom two items on the stack (y + x)
        """
        if self.use_modular_arithmetic(x, y):
            return (x + y) % config.cfg["modulus"]
        self.TypeCheck(x, y)
        try:
            return y + x
        except TypeError, e:
            raise e
        except Exception, e:
            self.errors.append(str(e))
            return x + y

    def subtract(self, y, x):
        """
    Usage: y x -

    Return the difference of the bottom two items on the stack (y - x)
        """
        if self.use_modular_arithmetic(x, y):
            return (x - y) % config.cfg["modulus"]
        self.TypeCheck(x, y)
        try:
            return y - x
        except:
            return -(y - x)

    def multiply(self, y, x):
        """
    Usage: y x *

    Return the product of the bottom two items on the stack (y * x)
        """
        if self.use_modular_arithmetic(y, x):
            return (y*x) % config.cfg["modulus"]
        self.TypeCheck(y, x)
        try:
            return y*x
        except:
            return x*y

    def divide(self, y, x):
        """
    Usage: y x /

    Return the quotient of the bottom two items on the stack (y / x)
        """
        if self.use_modular_arithmetic(x, y):
            return (y//x) % config.cfg["modulus"]
        self.TypeCheck(y, x)
        if x == 0:
            if config.cfg["allow_divide_by_zero"]:
                if y > 0:
                    return m.inf
                elif y < 0:
                    return -m.inf
                else:
                    raise ValueError("%s0/0 is ambiguous" % fln())
            else:
                raise ValueError("%sCan't divide by zero" % fln())
        if isint(y) and isint(x):
            if y.C_division or x.C_division:
                return y // x
            q = Rational(int(y), int(x))
            if q.d == 1:
                return q.n
            return q
        try:
            return y/x
        except:
            return self.reciprocal(x/y)

    def Mod(self, n, d):
        """
    Usage: y x %

    Return the modulus of the bottom two items on the stack (y mod x)
        """
        self.TypeCheck(n, d)
        if isint(n) and isint(d):
            return Zn(n) % Zn(d)
        if isinstance(n, Zn): n = n.value
        if isinstance(d, Zn): d = d.value
        n = Convert(n, MPF)
        d = Convert(d, MPF)
        result = m.modf(n, d)
        if int(result) == result:
            result = Zn(int(result))
        return result

    def integer_divide(self, n, d):
        """
    Usage: y x div

    Return the integer division quotient of the bottom two items on the stack (y // x)
        """
        if self.use_modular_arithmetic(n, d):
            return (Zn(n)//Zn(d)) % config.cfg["modulus"]
        self.TypeCheck(n, d)
        if isint(n) and isint(d):
            if not isinstance(n, Zn): n = Zn(n)
            if not isinstance(d, Zn): d = Zn(d)
            return n // d
        n = Convert(n, MPF)
        d = Convert(d, MPF)
        return int(m.floor(n/d))

    def bit_and(self, y, x):
        """
    Usage: y x &

    Return the bitwise AND of the bottom two items on the stack (y & x)
        """
        self.TypeCheck(y, x)
        if isint(y) and isint(x):
            if not isinstance(y, Zn): y = Zn(y)
            if not isinstance(x, Zn): x = Zn(x)
            return y & x
        y = Convert(y, INT)
        x = Convert(x, INT)
        return y & x

    def bit_or(self, y, x):
        """
    Usage: y x |

    Return the bitwise OR of the bottom two items on the stack (y | x)
        """
        self.TypeCheck(y, x)
        if isint(y) and isint(x):
            if not isinstance(y, Zn): y = Zn(y)
            if not isinstance(x, Zn): x = Zn(x)
            return y | x
        y = Convert(y, INT)
        x = Convert(x, INT)
        return y | x

    def bit_xor(self, y, x):
        """
    Usage: y x xor

    Return the bitwise XOR of the bottom two items on the stack (y XOR x)
        """
        self.TypeCheck(y, x)
        if isint(y) and isint(x):
            if not isinstance(y, Zn): y = Zn(y)
            if not isinstance(x, Zn): x = Zn(x)
            return y ^ x
        y = Convert(y, INT)
        x = Convert(x, INT)
        return y ^ x

    def bit_leftshift(self, y, x):
        """
    Usage: y x <<

    Return the bitwise left shift of the bottom two items on the stack (y << x)
        """
        self.TypeCheck(y, x)
        if isint(y) and isint(x):
            if not isinstance(y, Zn): y = Zn(y)
            if not isinstance(x, Zn): x = Zn(x)
            return y << x
        y = Convert(y, INT)
        x = Convert(x, INT)
        return y << x

    def bit_rightshift(self, y, x):
        """
    Usage: y x >>

    Return the bitwise right shift of the bottom two items on the stack (y >> x)
        """
        self.TypeCheck(y, x)
        if isint(y) and isint(x):
            if not isinstance(y, Zn): y = Zn(y)
            if not isinstance(x, Zn): x = Zn(x)
            return y >> x
        y = Convert(y, INT)
        x = Convert(x, INT)
        return y >> x

    def percent_change(self, y, x):
        """
    Usage: y x %ch

    Return the percent change between the bottom two items on the stack
        """
        y = Convert(y, MPF)
        x = Convert(x, MPF)
        if y == 0:
            raise ValueError("%sBase is zero for %ch" % fln())
        return 100*(x - y)/y

    def combination(self, y, x):
        """
    Usage: y x comb

    Return the statistical combination of the bottom two items on the stack
        """
        if (not config.cfg["coerce"]) and \
           (not isint(y)) and (not isint(x)):
            raise ValueError(self.argument_types % fln())
        y = Convert(y, INT)
        x = Convert(x, INT)
        return int(self.permutation(y, x)//self.Factorial(x))

    def permutation(self, y, x):
        """
    Usage: y x perm

    Return the statistical permutation of the bottom two items on the stack
        """
        if (not config.cfg["coerce"]) and \
           (not isint(y)) and (not isint(x)):
            raise ValueError(self.argument_types % fln())
        y = Convert(y, INT)
        x = Convert(x, INT)
        return int(self.Factorial(y)//self.Factorial(y - x))

    def power(self, y, x):
        """
    Usage: y x ^

    Return the value of the pow() function applied to the bottom two items on the stack (y^x)
        """
        return pow(y, x)

    #---------------------------------------------------------------------------
    # Unary functions

    def reciprocal(self, x):
        """
    Usage: x inv

    Returns the reciprocal of x (1/x)
        """
        if x == 0:
            if config.cfg["allow_divide_by_zero"]:
                return inf
            else:
                raise ValueError("%sDivision by zero" % fln())
        if isint(x):
            return Rational(1, x)
        elif isinstance(x, Rational):
            return Rational(x.d, x.n)
        return m.mpf(1)/x

    def bit_negate(self, x):
        """
    Usage: x ~

    Returns the bit-negated version of x (x gets cast to an int)
        """
        return ~Convert(x, INT)

    def negate(self, x):
        """
    Usage: x neg

    Returns negative x (-(x))
        """
        return -x

    def conj(self, x):
        """
    Usage: x conj

    Returns the complex conjugate of complex number x
        """
        if isinstance(x, m.mpc):
            return x.conjugate()
        else:
            return x

    def sqrt(self, x):
        """
    Usage: x sqrt

    Returns the square root of x
        """
        if isinstance(x, Zn): x = int(x)
        return m.sqrt(x)

    def cbrt(self, x):
        """
    Usage: x cbrt

    Returns the cube root of x
        """
        if isinstance(x, Zn): x = int(x)
        return m.cbrt(x)

    def root(self, y, x, k=0):
        """
    Usage: y x root

    Returns the xth root of y
        """
        if isinstance(x, Zn): x = int(x)
        if isinstance(y, Zn): y = int(y)
        return m.root(y, x, k)

    def roots(self, y, x):
        """
    Usage: y x roots

    Returns all the xth roots of y
        """
        return [ self.root(y, x, k) for k in xrange(x) ]

    def square(self, x):
        """
    Usage: x sqr

    Returns the square of x
        """
        return x*x

    def mid(self, x):
        """
    Usage: x mid

    Returns the midpoint for interval number x
        """
        if isinstance(x, m.ctx_iv.ivmpf):
            return x.mid
        else:
            raise ValueError("%sNeed an interval number for mid" % fln())

    def Factorial(self, x):
        """
    Usage: x !

    Returns the factorial of x.  This returns the exact factorial up to
    cfg['factorial_limit'] and a floating point approximation beyond that.
        """
        def ExactIntegerFactorial(x):
            if x in self.factorial_cache:
                return self.factorial_cache[x]
            else:
                if x > 2:
                    y = 1
                    for i in xrange(2, x+1):
                        y *= i
                    self.factorial_cache[x] = y
                    return y
        limit = config.cfg["factorial_limit"]
        if limit < 0 or not isint(limit):
            raise SyntaxError("%sFactorial limit needs to be an integer >= 0" % fln())
        if isint(x) and x >= 0:
            if limit == 0 or (limit > 0 and x < limit):
                return ExactIntegerFactorial(x)
        return m.factorial(x)

    def Sum(self, *args):
        """
    Usage: x sum

    Returns the sum of the bottom x items on the stack
    """
        s = 0
        try:
            for x in args:
                s = self.add(s, x)
        except Exception, e:
            self.display.msg("%sStack is not large enough, %s" % (fln(), e))
            return None
        return s

    def floor(self, x):
        """
    Usage: x floor

    Returns the next integer less than or equal to x
    """
        if isinstance(x, Zn): x = int(x)
        return m.floor(x)

    def ceil(self, x):
        """
    Usage: x ceil

    Returns the next integer greater than or equal to x
    """
        if isinstance(x, Zn): x = int(x)
        return m.ceil(x)

    def atan2(self, x, y):
        """
    Usage: y x atan2

    Returns the arc tangent of the angle with legs y and x (gets angle sign correct)
        """
        if isinstance(x, Zn): x = int(x)
        return m.atan2(x, y)

    def hypot(self, x, y):
        """
    Usage: y x hypot

    Returns the hypotenuse of the right triangle with legs x and y
        """
        if isinstance(x, Zn): x = int(x)
        return m.hypot(x, y)

    def sin(self, x):
        """
    Usage: x sin

    Returns the sine of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.sin(self.Conv2Rad(x))

    def cos(self, x):
        """
    Usage: x cos

    Returns the cosine of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.cos(self.Conv2Rad(x))

    def tan(self, x):
        """
    Usage: x tan

    Returns the tangent of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.tan(self.Conv2Rad(x))

    def asin(self, x):
        """
    Usage: x asin

    Returns the arc-sine of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.asin(x))

    def acos(self, x):
        """
    Usage: x acos

    Returns the arc-cosine of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.acos(x))

    def atan(self, x):
        """
    Usage: x atan

    Returns the arctangent of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.atan(x))

    def sec(self, x):
        """
    Usage: x sec

    Returns the secant of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.sec(self.Conv2Rad(x))

    def csc(self, x):
        """
    Usage: x csc

    Returns the cosecant of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.csc(self.Conv2Rad(x))

    def cot(self, x):
        """
    Usage: x cot

    Returns the cotangent of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.cot(self.Conv2Rad(x))

    def asec(self, x):
        """
    Usage: x asec

    Returns the arc-secant of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.asec(x))

    def acsc(self, x):
        """
    Usage: x acsc

    Returns the arc-cosecant of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.acsc(x))

    def acot(self, x):
        """
    Usage: x acot

    Returns the arc-cotangent of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.acot(x))

    def sinh(self, x):
        """
    Usage: x sinh

    Returns the hypebolic sine of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.sinh(self.Conv2Rad(x))

    def cosh(self, x):
        """
    Usage: x cosh

    Returns the hypebolic cosine of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.cosh(self.Conv2Rad(x))

    def tanh(self, x):
        """
    Usage: x tanh

    Returns the hypebolic tangent of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.tanh(self.Conv2Rad(x))

    def asinh(self, x):
        """
    Usage: x asinh

    Returns the hypebolic arc-sine of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.asinh(x))

    def acosh(self, x):
        """
    Usage: x acosh

    Returns the hypebolic arc-cosine of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.acosh(x))

    def atanh(self, x):
        """
    Usage: x atanh

    Returns the hypebolic arctangent of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.atanh(x))

    def sech(self, x):
        """
    Usage: x sech

    Returns the hyperbolic secant of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.sech(self.Conv2Rad(x))

    def csch(self, x):
        """
    Usage: x csch

    Returns the hyperbolic cosecant of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.csch(self.Conv2Rad(x))

    def coth(self, x):
        """
    Usage: x coth

    Returns the hyperbolic cotangent of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.coth(self.Conv2Rad(x))

    def asech(self, x):
        """
    Usage: x asech

    Returns the hyperbolic arc-secant of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.asech(x))

    def acsch(self, x):
        """
    Usage: x acsch

    Returns the hyperbolic arc-cosecant of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.acsch(x))

    def acoth(self, x):
        """
    Usage: x acoth

    Returns the hyperbolic arc-cotangent of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return self.Conv2Deg(m.acoth(x))

    def ln(self, x):
        """
    Usage: x ln

    Returns the natural log of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.ln(x)

    def log2(self, x):
        """
    Usage: x log2

    Returns the log2 of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.ln(x)/m.ln(2)

    def log10(self, x):
        """
    Usage: x log

    Returns the log10 of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.log10(x)

    def exp(self, x):
        """
    Usage: x exp

    Returns e raised to the power of top item on the stack (e^x)
        """
        if isinstance(x, Zn): x = int(x)
        return m.exp(x)

    def bits(self, x):
        """
    Usage: x bits

    Returns the number of bits required to represent x as an integer
        """
        if isinstance(x, Zn): x = int(x)
        if x < 0: x = self.abs(x)
        return int(self.ceil(self.log2(x)))

    def db(self, x):
        """
    Usage: x db

    Returns the x (in dB) as a standard ratio
    x = 10*log(y) ==> x db -> y
    or a shortcut for '10 x 10 / ^'
        """
        return self.power(mpf('10.0'), x / mpf('10.0'))

    def bd(self, x):
        """
    Usage: x bd

    Returns the x as dB
    y = 10*log(x) ==> x bd -> y
    or a shortcut for 'x log10 10 *'
        """
        return 10 * self.log10(x)

    def abs(self, x):
        """
    Usage: x abs

    Returns the absolute value of the top item on the stack
        """
        return abs(x)

    def arg(self, x):
        """
    Usage: x arg

    Returns the complex argument of the top item on the stack
        """
        if isinstance(x, Zn): x = int(x)
        return m.arg(x)

    def gamma(self, x):
        """
    Usage: x gamma

    Returns the gamma function at x
        """
        if isinstance(x, Zn): x = int(x)
        return m.gamma(x)

    def zeta(self, x):
        """
    Usage: x zeta

    Returns the zeta function at x
        """
        if isinstance(x, Zn): x = int(x)
        return m.zeta(x)

    def Ncdf(self, x):
        'Normal probability CDF'
        return ncdf(x, 0, 1)

    def Incdf(self, x):
        'Inverse of normal probability CDF'
        if not (0 < x < 1):
            raise ValueError("%sInverse normal CDF requires 0 < argument < 1" % \
                fln())

        if   x < 0.01:  start = mpf("-2.3")
        elif x < 0.05:  start = mpf("-1.64")
        elif x < 0.10:  start = mpf("-1.28")
        elif x < 0.25:  start = mpf("-0.67")
        elif x < 0.50:  start = mpf("0")
        elif x < 0.75:  start = mpf("0.67")
        elif x < 0.90:  start = mpf("1.28")
        elif x < 0.95:  start = mpf("1.64")
        elif x < 0.99:  start = mpf("2.3")
        else:           start = mpf("3")
        y = m.findroot(lambda z: ncdf(z) - x, 0)
        return y

    def rand(self):
        """
    Usage: rand

    Return a uniformly-distributed random number in [0, 1).  We use
    the os.urandom function to return a group of random bytes, then convert
    the bytes to a binary fraction expressed in decimal.
        """
        numbytes = ceil(mp.prec/mpf(8)) + 1
        bytes = os.urandom(numbytes)
        number = self.sum([ord(b)*mpf(256)**(-(i+1)) for i, b in enumerate(list(bytes))])
        return number

    def unix_ts(self):
        """
    Usage: ts

    Return the current Unix date/time as a float
    (use 2Jul to transform to a julian date)
        """
        return m.mpf(time.time())

    ############################################################################
    # constants.  Should these be handled differently?
    ############################################################################

    def Phi(self):
        """
    Usage: phi

    Returns Phi (the golden ratio)
        """
        return m.mpf(m.phi)

    def Pi(self):
        """
    Usage: pi

    Returns Pi
        """
        return m.mpf(m.mp.pi)

    def E(self):
        """
    Usage: e

    Returns e
        """
        return m.mpf(m.mp.e)

    def I(self):
        """
    Usage: i (or j)

    Returns i (or j, depending on how you see it)
        """
        return m.mpc(0,1)

    def choose_a_const(self, clist, keys):
        # read names out in groups of 24
        name = ''
        offset = 0
        prompt = "Type number to choose, enter for (%d) more, q to quit: "
        winnowed_size = len(clist)
        while name == '' and offset < winnowed_size:
            options = []
            top = min(offset + 24, winnowed_size) - offset
            for i in xrange(1, top+1):
                ck = keys[clist[offset+i-1][1]]
                options.append("% 2d: %s" % (i, ck))
            self.display.msg('\n'.join(options))
            name = self.chomp(raw_input(prompt % (winnowed_size - (offset + top))))
            if name != '':
                if name == 'q':
                    return None
                try:
                    name = int(name)
                    if name > 0 and name < top:
                        lcname, idx = clist[offset + name - 1]
                        # print "choosing idx=%d -> %s (%d)" % (idx, lcname, name)
                        return self.constants[keys[idx]]
                    self.errors.append("Invalid selection: %d"%name)
                    return None
                except ValueError, e:
                    self.errors.append("Invalid selection: %s"%name)
                    return None
            offset += 24
        return None

    def const(self, line=''):
        """
    Usage: const

    Presents a list of constants to choose from
        """
        ord_names = self.constants.keys()
        ord_names.sort(key=str.lower)
        lc_names = [ n.lower() for n in ord_names ]

        line = self.chomp(line)
        if line == '':
            name = self.chomp(raw_input("Type constant name (or part) or * to list: "))
            if name == '': return None
        else:
            name = line
        if name == '*':
            winnowed = [ (lc_names[i], i) for i in xrange(len(lc_names)) ]
            return self.choose_a_const(winnowed, ord_names)
        else:
            lc_name = name.lower()
            winnowed = [ (lc_names[i], i) for i in xrange(len(lc_names)) if lc_name in lc_names[i] ]
            if len(winnowed) == 0:
                self.errors.append("No constant like '%s'"%name)
                return None
            if len(winnowed) == 1:
                return self.constants[ord_names[winnowed[0][1]]]
            return self.choose_a_const(winnowed, ord_names)

    ############################################################################
    # networking functions
    ############################################################################

    def ntohl(self, x):
        """
    Usage: x ntohl
           x htonl

    Changes int from net (big-endian) to host order and back again
        """
        if not isint(x):
            raise TypeError("ntohl requires an integer argument")
        return htonl(int(x))

    def netmask(self, y, x):
        """
    Usage: y x netmask

    Adds netmask information to an IPv4 address.  In the netmask form,
    use something like '192.168.2.50 255.255.255.0 netmask'.  Do not
    use this form for IPv6 addresses, use cidr.

    An alternative is to enter the IP address with the netmask
    already appended in cidr format such as '192.168.2.50/24'.
        """
        if not isint(x):
            raise TypeError("netmask requires an integer argument")
        if not isint(y):
            raise TypeError("cidr requires an integer (or IP) argument for IP address")
        if not isinstance(x, ipaddr):
            x1 = ipaddr(x, ipvn='ipv4')
        x1 = ~x
        z = self.bits(x1.value)
        if (1 << z)-1 != x1.value:
            raise ValueError("Invalid netmask %s" % str(x))
        return ipaddr(y, (32-z))

    def cidr(self, y, x):
        """
    Usage: y x cidr

    Adds netmask information to an IP address.  In the netmask form,
    use something like '192.168.2.50 24 cidr'.

    An alternative is to enter the IP address with the netmask
    already appended in cidr format such as '192.168.2.50/24'.
        """
        if not isint(x):
            raise TypeError("cidr requires an integer argument for netmask")
        if not isint(y):
            raise TypeError("cidr requires an integer (or IP) argument for IP address")
        return ipaddr(y, x)

    def samenet(self, y, x):
        pass


    ############################################################################
    # Stack callback functions
    ############################################################################

    def ClearStack(self):
        """
    Usage: clear

    Clear the stack, but keep register and other settings
        """
        self.stack.clear_stack()

    def Reset(self):
        """
    Usage: reset

    Reset the calculator to initial settings, clear stack, reset registers, etc.
        """
        self.ClearRegisters()
        self.ClearStack()
        config.load()
        self.ConfigChanged()

    def SetStackDisplay(self, x):
        """
    Usage: n stack

    Set the display size of the stack to n items.  Items beyond n are still
    saved as part of the stack, but are not displayed on a refresh.
        """
        msg = "Stack display size be an integer >= 0"
        if int(x) == x:
            if x >= 0:
                config.cfg["stack_display"] = int(x)
                return None
            else:
                self.display.msg(msg)
                return x
        else:
            self.display.msg(msg)
            return x

    def lastx(self):
        """
    Usage: lastx

    Push the saved last x back onto the stack
        """
        assert self.stack.lastx != None, "Bug:  stack.lastx is None"
        return self.stack.lastx


    def swap(self):
        """
    Usage: swap

    Swap the bottom two stack items
        """
        try:
            self.stack.swap()
        except:
            self.display.msg("%sStack is not large enough" % fln())

    def roll(self):
        """
    Usage: roll

    Roll the stack up (1 => 2, 2 => 3, n => 1)
        """
        self.stack.roll(0)

    def rolld(self):
        """
    Usage: rolld

    Roll the stack down (1 => n, 2 => 1, n => n -  1)
        """
        self.stack.roll(-1)

    def over(self):
        """
    Usage: over

    Pushes the second to bottom item onto the stack as the bottom
        """
        return self.stack[1]

    def pick(self, x):
        """
    Usage: n pick

    Pushes the nth to bottom item onto the stack as the bottom
        """
        x = int(x)
        return self.stack[x-1]

    def drop(self, x):
        """
    Usage: drop

    Drops the bottom item off the stack
        """
        return None

    def drop2(self, y, x):
        """
    Usage: drop2

    Drops the bottom two items off the stack
        """
        return None

    def dropn(self, *args):
        """
    Usage: n dropn

    Drops the bottom n items off the stack
        """
        return None

    def dup(self, x):
        """
    Usage: dup

    Duplicates the bottom item on the stack
        """
        return x, x

    def dup2(self, y, x):
        """
    Usage: dup2

    Duplicates the bottom two items on the stack
        """
        return y, x, y, x

    def dupn(self, *args):
        """
    Usage: n dupn

    Duplicates the bottom n items on the stack
        """
        return args+args

    def depth(self):
        """
    Usage: depth

    Pushes the stack depth onto the bottom of the stack
        """
        return len(self.stack)

    ############################################################################
    # Casting and converting functions
    ############################################################################

    def Cast(self, x, newtype, use_prec=False):
        try:
            try:
                if use_prec == True:
                    digits = 0
                else:
                    digits = max(1, self.fp.num_digits)
            except:
                pass
            return Convert(x, newtype, digits)
        except Exception, e:
            self.display.msg("%sCouldn't perform conversion" % fln())
            raise e

    def Cast_i(self, x):
        """
    Usage: x I

    Returns x casted as a integer value
        """
        return self.Cast(x, INT)

    def Cast_qq(self, x):
        """
    Usage: x QQ

    Returns x casted as a rational value (using displayed precision)
        """
        return self.Cast(x, RAT, use_prec=True)

    def Cast_q(self, x):
        """
    Usage: x Q

    Returns x casted as a rational value
        """
        return self.Cast(x, RAT, use_prec=False)

    def Cast_r(self, x):
        """
    Usage: x R

    Returns x casted as a real value
        """
        return self.Cast(x, MPF)

    def Cast_c(self, x):
        """
    Usage: x C

    Returns x casted as a complex value
        """
        return self.Cast(x, MPC)

    def Cast_t(self, x):
        """
    Usage: x T

    Returns x casted as a Julian value
        """
        return self.Cast(x, JUL)

    def Cast_v(self, x):
        """
    Usage: x V

    Returns x casted as an interval value
        """
        return self.Cast(x, MPI)

    def cast(self, x):
        """
    Usage: x cast

    Returns x casted as an int of the current "C" stype integer type
        """
        return Zn(int(x))

    def IP(self, x):
        """
    Usage: IP

    Convert integer to an IP address
        """
        return ipaddr(x)

    def ToDegrees(self, x):
        """
    Usage: x 2deg

    Returns x converted from radians to degrees
        """
        if x == 0: return x
        if isinstance(x, m.mpc):
            raise ValueError("%sNot an appropriate operation for a complex number" % fln())
        return degrees(x)

    def ToRadians(self, x):
        """
    Usage: x 2rad

    Returns x converted from degrees to radians
        """
        if x == 0: return x
        if isinstance(x, m.mpc):
            raise ValueError("%sNot an appropriate operation for a complex number" % fln())
        return radians(x)

    def ToUnix(self, x):
        """
    Usage: x unix

    Returns x (which must be a Julian date) as a Unix timestamp
        """
        if not isinstance(x, Julian):
            raise ValueError("%sThis function requires a Julian date (use T?)" % fln())
        utc_offset = time.mktime(time.localtime()) - time.mktime(time.gmtime())
        if time.daylight:
            utc_offset += 3600
        return (self.Cast_r(x-JULIAN_UNIX_EPOCH))*86400-utc_offset

    def ToJulian(self, x):
        """
    Usage: x julian

    Returns x (interpreted as a Unix timestamp) as a Julian date
        """
        utc_offset = time.mktime(time.localtime()) - time.mktime(time.gmtime())
        if time.daylight:
            utc_offset += 3600
        return Julian((self.Cast_r(x)+utc_offset)/86400)+JULIAN_UNIX_EPOCH

    def hr(self, x):
        """
    Usage: x hr

    Convert hms to decimal hours
        """
        x = Convert(x, MPF)
        hours = int(x)
        x -= hours
        x *= 100
        minutes = int(x)
        x -= minutes
        x *= 100
        return hours + minutes/mpf(60) + x/3600

    def hms(self, x):
        """
    Usage: x hms

    Convert decimal hours to hours.MMSSss
        """
        x = Convert(x, MPF)
        hours = int(x)
        x -= hours
        x *= 60
        minutes = int(x)
        if minutes == 60:
            hours += 1
            minutes = 0
        x -= minutes
        seconds = 60*x
        if seconds == 60:
            minutes += 1
            seconds = 0
        if minutes == 60:
            hours += 1
            minutes = 0
        return hours + minutes/mpf(100) + seconds/mpf(10000)

    def split(self, x):
        """
    Usage: x split

    Returns the respective parts of various compound numbers:
    complex, rational, interval, Julian, vectors and floats
        """
        if isinstance(x, mpf):
            ip = Convert(x, INT)
            fp = Convert(x, MPF)
            return ip, fp - ip.value
        if isinstance(x, mpc):
            return x.real, x.imag
        elif isinstance(x, Rational):
            return x.n, x.d
        elif isinstance(x, ctx_iv.ivmpf):
            return mpf(x.a), mpf(x.b)
        elif isinstance(x, Julian):
            return mpf(x.value.a), mpf(x.value.b)
        else:
            msg = "%sapart requires rational, complex, or interval number"
            raise TypeError(msg % fln())

    def first_part(self, x):
        """
    Usage: x fp

    Returns the first part of any compound number:
        Floating point => integer part
        Complex        => real part
        Rational       => numerator
        Interval       => lower boundary
        Julian         => first date
        """
        a,b = self.split(x)
        return a

    def second_part(self, x):
        """
    Usage: x sp

    Returns the second part of any compound number:
        Floating point => float part
        Complex        => imaginary part
        Rational       => denominator
        Interval       => upper boundary
        Julian         => last date
        """
        a,b = self.split(x)
        return b

    def ToIV(self, y, x):
        """
    Usage: y x iv

    Convert to interval number [y,x]
        """
        if y > x:
            msg = "%sy must be <= x"
            raise ValueError(msg % fln())
        y = Convert(y, MPF)
        x = Convert(x, MPF)
        return mpi(y, x)

    def gcf(self, y, x):
        """
    Usage: y x gcf

    Returns the greatest common factor of x and y
        """
        def subgcf(a, b):
            if b == 0: return a
            return subgcf(b, a % b)

        if not isint(x) or not isint(y):
            raise TypeError("operands to gcf must be integers")
        if y > x:
            return subgcf(y, x)
        return subgcf(x, y)

    def lcd(self, y, x):
        """
    Usage: y x lcd

    Returns the lowest common denominator of x and y
        """
        if not isint(x) or not isint(y):
            raise TypeError("operands to lcd must be integers")
        return self.multiply(y, x)/self.gcf(y, x)

    def Chop(self, x):
        """
    Usage: x chop

    Returns the the value of x as displayed
        """
        return self.number(self.Format(x).replace(" ", ""))

    def Prec(self, x):
        """
    Usage: x prec

    Set floating point precision to x digits
        """
        if isint(x) and x > 0:
            mp.dps = int(x)
            config.cfg["prec"] = int(x)
            if config.cfg["fp_digits"] > mp.dps:
                config.cfg["fp_digits"] = mp.dps
            if self.fp.num_digits > mp.dps:
                self.fp.digits(mp.dps)
            return None
        else:
            self.display.msg("You must supply an integer > 0")

    def digits(self, x):
        """
    Usage: x digits

    Set floating point display to x digits
        """
        if int(x) == x:
            if x >= 0:
                d = min(int(x), mp.dps)
                config.cfg["fp_digits"] = d
                self.fp.digits(min(int(x), mp.dps))
                return None
            else:
                self.display.msg("Use an integer >= 0")
        else:
            self.display.msg("You must supply an integer >= 0")

    def Round(self, y, x):
        """
    Usage: y x round

    Round y to the nearest x.  Algorithm from PC Magazine, 31Oct1988, pg 435.
        """
        y = Convert(y, MPF)
        x = Convert(x, MPF)
        sgn = 1
        if y < 0: sgn = -1
        return sgn*int(mpf("0.5") + abs(y)/x)*x

    def In(self, y, x):
        """
    Usage: y x in

    Both x and y are expected to be interval numbers or x a number
    and y and interval number.  Returns the boolean 'x in y'.
        """
        msg = "%sy needs to be an interval number or Julian interval"
        if isinstance(y, m.ctx_iv.ivmpf):
            if self.testing:
                if x not in y:
                    s = "x was not in y:" + nl
                    s += "  x = " + repr(x) + nl
                    s += "  y = " + repr(y)
                    self.display.msg(s)
                    exit(1)
            return x in y
        elif isinstance(y, Julian):
            if not isinstance(y.value, m.ctx_iv.ivmpf):
                raise ValueError(msg % fln())
            if self.testing:
                if x not in y.value:
                    s = "x was not in y:" + nl
                    s += "  x = " + repr(x) + nl
                    s += "  y = " + repr(y.value)
                    self.display.msg(s)
                    exit(1)
            if isinstance(x, Julian):
                return x.value in y.value
            else:
                return x in y.value
        else:
            raise ValueError(msg % fln())

    ############################################################################
    # Display modification functions
    ############################################################################

    def mixed(self, x):
        """
    Usage: x mixed

    Show the rationals as mixed fractions or not
        """
        if x != 0:
            config.cfg["mixed_fractions"] = True
            Rational.mixed = True
        else:
            config.cfg["mixed_fractions"] = False
            Rational.mixed = False

    def Debug(self, x):
        """
    Usage: x debug

    Set or clear the debug flag based on x
        """
        if x != 0:
            debug(True)
        else:
            debug(False)

    def Show(self):
        """
    Usage: x show

    Show the full precision of the bottom value on the stack
        """
        def showx(x, prefix=""):
            if mp.dps < 2:
                return
            sign, mant, exponent = to_digits_exp(x._mpf_, mp.dps)
            s = mant[0] + "." + mant[1:] + "e" + str(exponent)
            if sign ==  -1:
                s = "-" + s
            self.display.msg(" " + prefix + s)
        from mpmath.libmp.libmpf import to_digits_exp
        x = self.stack[0]
        if isinstance(x, m.mpf):
            showx(x)
        elif isinstance(x, m.mpc):
            self.display.msg(" x is complex")
            showx(x.real, "  x.real:  ")
            showx(x.imag, "  x.imag:  ")
        elif isinstance(x, m.ctx_iv.ivmpf):
            self.display.msg(" x is an interval number")
            showx(x.a, "  x.a:  ")
            showx(x.b, "  x.b:  ")

    def comma(self, x):
        """
    Usage: x comma

    If x, use commas to decorate displayed values
        """
        if x != 0:
            config.cfg["fp_comma_decorate"] = True
        else:
            config.cfg["fp_comma_decorate"] = False
        mpFormat.comma_decorate = config.cfg["fp_comma_decorate"]

    def width(self, x):
        """
    Usage: x width

    Set display width to x (x must be > 20)
        """
        if isint(x) and x > 20:
            config.cfg["line_width"] = int(x)
        else:
            self.display.msg("width command requires an integer > 20")

    def Rectangular(self):
        """
    Usage: rec

    Set rectangular mode for display of complex numbers and vectors
        """
        config.cfg["imaginary_mode"] = "rect"

    def Polar(self):
        """
    Usage: polar

    Set polar mode for display of complex numbers and vectors
        """
        config.cfg["imaginary_mode"] = "polar"

    def fix(self):
        """
    Usage: fix

    Set fixed-point mode for display of floating point numbers
        """
        config.cfg["fp_format"] = "fix"

    def sig(self):
        """
    Usage: sig

    Set significant digits mode for display of floating point numbers
        """
        config.cfg["fp_format"] = "sig"

    def sci(self):
        """
    Usage: sci

    Set scientific mode for display of floating point numbers
        """
        config.cfg["fp_format"] = "sci"

    def eng(self):
        """
    Usage: eng

    Set engineering mode for display of floating point numbers
        """
        config.cfg["fp_format"] = "eng"

    def engsi(self):
        """
    Usage: eng

    Set engineering mode for display of floating point numbers
        """
        config.cfg["fp_format"] = "engsi"

    def raw(self):
        """
    Usage: raw

    Set raw mode for display of floating point numbers
        """
        config.cfg["fp_format"] = "raw"

    def dec(self):
        """
    Usage: dec

    Set decimal mode for display of integers
        """
        config.cfg["integer_mode"] = "dec"

    def hex(self):
        """
    Usage: hex

    Set hexadecimal mode for display of integers
        """
        config.cfg["integer_mode"] = "hex"

    def oct(self):
        """
    Usage: oct

    Set octal mode for display of integers
        """
        config.cfg["integer_mode"] = "oct"

    def bin(self):
        """
    Usage: bin

    Set binary mode for display of integers
        """
        config.cfg["integer_mode"] = "bin"

    def roman(self):
        """
    Usage: roman

    Set roman numeral mode for display of integers
        """
        config.cfg["integer_mode"] = "roman"

    def iva(self):
        """
    Usage: iva

    Set interval mode A for display of intervals
        """
        config.cfg["iv_mode"] = "a"
        Julian.interval_representation = "a"

    def ivb(self):
        """
    Usage: ivb

    Set interval mode B for display of intervals
        """
        config.cfg["iv_mode"] = "b"
        Julian.interval_representation = "b"

    def ivc(self):
        """
    Usage: ivc

    Set interval mode C for display of intervals
        """
        config.cfg["iv_mode"] = "c"
        Julian.interval_representation = "c"

    def on(self):
        """
    Usage: on

    Set display output on
        """
        self.display.on()
        return status_ok_no_display

    def off(self):
        """
    Usage: off

    Set display output off
        """
        self.display.off()
        return status_ok_no_display

    def deg(self):
        """
    Usage: deg

    Set angle mode to degrees.  This means that all values
    used by the various trigonometric functions are assumed
    to be already expressed in degrees.  To do this, all
    values are converted behind the scenes to radians before
    passing them to the functions
        """
        config.cfg["angle_mode"] = "deg"

    def rad(self):
        """
    Usage: rad

    Set angle mode to radians.  This means that all values
    used by the various trigonometric functions are assumed
    to be already expressed in radians.
        """
        config.cfg["angle_mode"] = "rad"

    def Rationals(self, x):
        """
    Usage: x rat

    If x, show rationals as rationals instead of decimals
        """
        if x != 0:
            config.cfg["no_rationals"] = True
        else:
            config.cfg["no_rationals"] = False

    def ToggleDowncasting(self, x):
        """
    Usage: x down

    Toggle downcasting: if X, downcast floats to ints if precision permits
        """
        if x != 0:
            config.cfg["downcasting"] = True
        else:
            config.cfg["downcasting"] = False

    #---------------------------------------------------------------------------
    # Other functions

    def Modulus(self, x):
        """
    Usage: modulo

    Set up modulus arithmetic with X as the modulus (1 or 0 to cancel)
        """
        if isinstance(x, m.mpc) or isinstance(x, m.ctx_iv.ivmpf):
            raise ValueError("%sModulus cannot be a complex or interval number" % fln())
        if x == 0:
            config.cfg["modulus"] = 1
        else:
            config.cfg["modulus"] = x
        return None

    def ClearRegisters(self):
        """
    Usage: clrg

    Clears all registers
        """
        self.registers = {}

    def ShowConfig(self):
        """
    Usage: cfg

    Shows the current config
        """
        d = {True:"on", False:"off"}
        per = d[config.cfg["persist"]]
        st = str(config.cfg["stack_display"])
        lw = config.cfg["line_width"]
        mf = str(config.cfg["mixed_fractions"])
        dc = d[config.cfg["downcasting"]]
        sps = d[config.cfg["fp_show_plus_sign"]]
        am = config.cfg["angle_mode"]
        im = config.cfg["integer_mode"]
        imm = config.cfg["imaginary_mode"]
        sd = str(config.cfg["stack_display"])
        fmt = config.cfg["fp_format"]
        dig = str(config.cfg["fp_digits"])
        ad = str(config.cfg["arg_digits"])
        af = config.cfg["arg_format"]
        pr = str(mp.dps)
        br = d[config.cfg["brief"]]
        nr = d[config.cfg["no_rationals"]]
        cd = d[config.cfg["fp_comma_decorate"]]
        adz = d[config.cfg["allow_divide_by_zero"]]
        iv = config.cfg["iv_mode"]
        cdiv = d[config.cfg["C_division"]]
        dbg = d[get_debug()]
        if 1:
            s = '''Configuration:
      Stack:%(st)s    Commas:%(cd)s   +sign:%(sps)s   Allow divide by zero:%(adz)s
      iv%(iv)s    brief:%(br)s  C-type integer division:%(cdiv)s   Rationals:%(nr)s
      Line width:%(lw)s    Mixed fractions:%(mf)s     Downcasting:%(dc)s
      Complex numbers:  %(imm)s    arguments:%(af)s %(ad)s digits Debug:%(dbg)s
      Display: %(fmt)s %(dig)s digits   prec:%(pr)s  Integers:%(im)s  Angles:%(am)s''' \
        % locals()

        self.display.msg(s)

    def brief(self, x):
        """
    Usage: x brief

    Set display to truncate long numbers to one line (shown with ...)
        """
        if x != 0:
            config.cfg["brief"] = True
        else:
            config.cfg["brief"] = False

    ############################################################################
    # End of callback functions
    ############################################################################

    def Flatten(self, n=0):
        '''The top of the stack contains a sequence of values.  Remove them
        and put them onto the stack.  If n is 0, it means do this with all the
        elements; otherwise, just do it with the left-most n elements and
        discard the rest.
        '''
        assert n >= 0
        if not len(self.stack):
            return
        items = list(self.stack.pop())
        if not n:
            n = len(items)
        self.stack.stack += items[:n]

    def ConfigChanged(self):
        try:
            self.fp.digits(config.cfg["fp_digits"])
        except:
            raise ValueError("%s'fp_digits' value in configuration is bad" % fln())
        try:
            self.ap.digits(config.cfg["arg_digits"])
        except:
            raise ValueError("%s'arg_digits' value in configuration is bad" % fln())
        mpFormat.comma_decorate = config.cfg["fp_comma_decorate"]
        mpFormat.cuddle_si = config.cfg["fp_cuddle_si"]
        mpFormat.explicit_plus_sign = config.cfg["fp_show_plus_sign"]
        Rational.mixed = config.cfg["mixed_fractions"]
        Zn().C_division = config.cfg["C_division"]
        if isint(config.cfg["prec"]) and int(config.cfg["prec"]) > 0:
            mp.dps = config.cfg["prec"]
        else:
            raise ValueError("%s'prec' value in configuration is bad" % fln())

    def GetFullPath(self, s):
        '''If s doesn't have a slash in it, prepend it with the directory where
        our executable is.  If it does have a slash, verify it's usable.
        '''
        if "/" not in s:
            path, file = os.path.split(sys.argv[0])
            return os.path.join(path, s)
        else:
            return os.normalize(s)

    def SaveConfiguration(self):
        if config.cfg["persist"]:
            c = os.path.expanduser(os.path.join("~", ".config", "hc", "config"))
            msg = "%sCould not write %s to:\n  %s"
            try:
                WriteSettings(c, config)
            except:
                self.display.msg(msg % (fln(), "config", c))
        if config.cfg["persist_registers"]:
            r = os.path.expanduser(os.path.join("~", ".config", "hc", "registers"))
            try:
                WriteSettings(r, self.registers)
            except:
                self.display.msg(msg % (fln(), "registers", r))
        if config.cfg["persist_stack"]:
            s = os.path.expanduser(os.path.join("~", ".config", "hc", "stack"))
            try:
                WriteList(s, self.stack.stack)
            except:
                self.display.msg(msg % (fln(), "stack", s))

    def GetLineWidth(self):
        config.cfg["line_width"],height = console.size()

    def GetConfiguration(self):
        import config
        self.GetLineWidth()
        self.ConfigChanged()
        us = "Using default configuration"
        if config.cfg["persist"]:
            c, r, s = config.cfg["config_file"], config.cfg["config_save_registers"], \
                      config.cfg["config_save_stack"]
            if c and not self.use_default_config_only:
                try:
                    d = {}
                    p = GetFullPath(c)
                    execfile(p, d, d)
                    config = d["cfg"]
                    self.ConfigChanged()
                except:
                    msg = "%sCould not read and execute configuration file:" % fln() + \
                          nl + "  " + c
                    self.display.msg(msg)
                    self.display.msg(us)
            if r:
                try:
                    d = {}
                    p = GetFullPath(r)
                    execfile(p, d, d)
                    self.registers = d["registers"]
                except:
                    msg = "%sCould not read and execute register file:"  % fln() + \
                          nl + "  " + r
                    self.display.msg(msg)
            if s:
                try:
                    d = {}
                    p = GetFullPath(s)
                    execfile(p, d, d)
                    global stack
                    self.stack.stack = d["mystack"]
                except:
                    msg = "%sCould not read and execute stack file:" % fln() + \
                          nl + "  " + s
                    self.display.msg(msg)

    def DisplayStack(self):
        size = config.cfg["stack_display"]
        assert size >= 0 and isint(size)
        stack = self.stack._string(self.Format, size)
        if len(stack) > 0:
            self.display.msg(stack)
        if config.cfg["modulus"] != 1:
            self.display.msg(" (mod " + self.Format(config.cfg["modulus"])+ ")")
        if len(self.errors) > 0:
            self.display.msg("\n".join(self.errors))
            self.errors = []

    def EllipsizeString(self, s, desired_length, ellipsis):
        '''Remove characters out of the middle of s until the length of s
        with the ellipsis inserted is <= the desired length.  Note:  the
        string returned will be of length desired_length or one less.
        '''
        had_exponent = "e" in s or "E" in s
        had_dp = "." in s
        if len(s) <= desired_length:
            return s
        if len(s) < desired_length - len(ellipsis) + 3:
            raise Exception("%sProgram bug:  string too short" % fln())
        left, right = s[:len(s)//2], s[len(s)//2:]
        while len(left) + len(right) + len(ellipsis) > desired_length:
            try:
                left = left[:-1]
                right = right[1:]
            except:
                raise Exception("%sProgram bug:  string too short" % fln())
        new_s = left + ellipsis + right
        chopped_exponent = had_exponent and ("E" not in new_s and "e" not in new_s)
        chopped_dp = had_dp and "." not in new_s
        if chopped_exponent or chopped_dp:
            self.display.msg("Warning:  floating point number was 'damaged' by inserting ellipsis")
        return new_s

    def Format(self, x, item_is_x=True):
        '''Format the four different types of numbers.  Return a string in
        the proper format.  The item_is_x arg is because we never want to
        ellipsize x, no matter the size.  This is passed in by the stack display
        function as it processes the stack.
        '''
        width = abs(config.cfg["line_width"])
        brief = config.cfg["brief"] and not item_is_x
        e = config.cfg["ellipsis"]
        im = config.cfg["integer_mode"]
        stack_header_allowance = 5
        if isinstance(x, ipaddr):
            s = str(x)
        elif isint(x):
            if isint_native(x):
                x = Zn(x)
            if im == "dec":
                s = str(x)
            elif im == "hex":
                s = hex(x)
            elif im == "oct":
                s = oct(x)
            elif im == "bin":
                s = x.bin()
            elif im == "roman":
                s = x.roman()
            else:
                raise Exception("%s'%s' integer mode is unrecognized" % (im, fln()))
            # Prepend a space or + if this is being done in the mpFormat
            # object.  This is a hack; eventually, there will be a single
            # number object where the formatting is handled.
            if x >= Zn(0):
                if mpFormat.implicit_plus_sign == True:  sign = " "
                if mpFormat.explicit_plus_sign == True:  sign = "+"
                s = sign + s
            if s[-1] == "L": s = s[:-1]  # Handle old python longs
            if brief:
                s = self.EllipsizeString(s, width - stack_header_allowance, e)
            return s
        elif isinstance(x, Rational):
            if config.cfg["no_rationals"]:
                x = Zn(x.n)/Zn(x.d)
                s = self.fp.format(x, config.cfg["fp_format"])
            else:
                s = str(x)
                if x >= Rational(0):
                    if mpFormat.implicit_plus_sign == True:  sign = " "
                    if mpFormat.explicit_plus_sign == True:  sign = "+"
                    s = sign + s
            if len(s) > width//2:
                s = s.replace("/", " / ") # Makes / easier to see
            if brief:
                size = (width - stack_header_allowance)//2 - 1
                s = self.EllipsizeString(s, size, e)
            return s
        elif isinstance(x, mpf):
            s = self.fp.format(x, config.cfg["fp_format"])
            if s[-1] == ".": s = s[:-1]  # Remove a trailing dot
            if brief:
                s = self.EllipsizeString(s, width - stack_header_allowance, e)
            return s
        elif isinstance(x, mpc):
            space = config.cfg["imaginary_space"]
            s = ""
            if space:
                s = " "
            sre = self.fp.format(x.real, config.cfg["fp_format"])
            sim = self.fp.format(abs(x.imag), config.cfg["fp_format"]).strip()
            if config.cfg["ordered_pair"]:
                if brief:
                    size = (width - stack_header_allowance)//2 - 4
                    sre = self.EllipsizeString(sre, size, e).strip()
                    sim = self.EllipsizeString(sim, size, e)
                s = "(" + sre + "," + s + sim + ")"
            else:
                mode = config.cfg["imaginary_mode"]
                first = config.cfg["imaginary_unit_first"]
                unit = config.cfg["imaginary_unit"]
                if mode == "polar":
                    # Polar mode
                    sep = config.cfg["polar_separator"]
                    angle_mode = config.cfg["angle_mode"]
                    mag = abs(x)
                    ang = arg(x)
                    if angle_mode == "deg":
                        ang_sym = config.cfg["degree_symbol"]
                        ang *= 180/pi
                    else:
                        ang_sym = "rad"
                    if angle_mode != "deg" and angle_mode != "rad":
                        self.display.msg("Warning:  bad angle_mode('%s') in configuration" \
                            % angle_mode)
                    m = str(mag)
                    a = str(ang)
                    if config.cfg["fp_format"] != "raw":
                        m = self.fp.format(mag, config.cfg["fp_format"])
                        a = self.ap.format(ang, config.cfg["arg_format"])
                    if brief:
                        size = (width - stack_header_allowance)//2 - \
                               len(sep) - 4
                        mag = self.EllipsizeString(m, size, e)
                        ang = self.EllipsizeString(a, size, e)
                    s = m + config.cfg["polar_separator"] + a + " " + ang_sym
                else:
                    # Rectangular mode
                    if brief:
                        size = (width - stack_header_allowance)//2 - 1
                        sre = self.EllipsizeString(sre, size, e)
                        sim = self.EllipsizeString(sim, size, e)
                    if x.real == 0:
                        # Pure imaginary
                        sign = ""
                        if mpFormat.implicit_plus_sign == True:  sign = " "
                        if mpFormat.explicit_plus_sign == True:  sign = "+"
                        if x.imag < 0:
                            if x.imag == -1:
                                s = "-" + unit
                            else:
                                if first:
                                    s = "-" + unit + sim
                                else:
                                    s = "-" + sim + unit
                        elif x.imag == 0:
                            s = sign + sre
                        else:
                            if x.imag == 1:
                                s = sign + unit
                            else:
                                if first:
                                    s = sign + unit + sim
                                else:
                                    s = sign + sim + unit
                    else:
                        if x.imag < 0:
                            if space:
                                if first:
                                    s = sre + s + "-" + s + unit + sim
                                else:
                                    s = sre + s + "-" + s + sim + unit
                            else:
                                if first:
                                    s = sre + "-" + unit + sim
                                else:
                                    s = sre + "-" + sim + unit
                        elif x.imag == 0:
                            s = sre
                        else:
                            if space:
                                if first:
                                    s = sre + s + "+" + s + unit + sim
                                else:
                                    s = sre + s + "+" + s + sim + unit
                            else:
                                if first:
                                    s = sre + "+" + unit + sim
                                else:
                                    s = sre + "+" + sim + unit
                if mode != "rect" and mode != "polar":
                    self.display.msg("Warning:  bad imaginary_mode('%s') in configuration" \
                        % mode)
            return s
        elif isinstance(x, ctx_iv.ivmpf):
            a = mpf(x.a)
            b = mpf(x.b)
            mid = mpf(x.mid)
            delta = mpf(x.delta)/2
            f = config.cfg["fp_format"]
            mode = config.cfg["iv_mode"]
            sp = ""
            if config.cfg["iv_space"]:
                sp = " "
            if mode == "a":
                mid = self.fp.format(mid, f)
                delta = self.fp.format(delta, f).strip()
                s = mid + sp + "+-" + sp + delta
            elif mode == "b":
                if mid != 0:
                    pct = 100*delta/mid
                else:
                    pct = mpf(0)
                mid = self.fp.format(mid, f)
                pct = self.fp.format(pct, f).strip()
                s = mid + sp + "(" + pct + "%)"
            elif mode == "c":
                a = self.fp.format(a, f)
                b = self.fp.format(b, f).strip()
                br1, br2 = config.cfg["iv_brackets"]
                s = br1 + a.strip() + "," + sp + b + br2
            else:
                raise ValueError("%s'%s' is unknown iv_mode in configuration" % \
                    (fln(), mode))
        elif isinstance(x, Julian):
            return str(x)
        else:
            self.errors.append("%sError in Format():  Unknown number format" % fln())
            return str(x)
        return s

    def WriteList(self, filename, name, list):
        try:
            f = open(filename, "wb")
            p = f.write
            p("from mpmath import *" + nl)
            p("from rational import Rational" + nl)
            p("from integer import Zn" + nl)
            p("from julian import Julian" + nl)
            p("mp.dps = " + str(mp.dps) + nl + nl)
            p(name + " = [" + nl)
            indent = "  "
            for item in list:
                s = repr(item)
                if s == "<pi: 3.14159~>": s = "pi"
                p(indent + s + "," + nl)
            p("]" + nl)
            f.close()
        except Exception, e:
            msg = ("%sError trying to write list '%s':" % (fln(), name)) + nl + str(e)
            self.display.msg(msg)
            raise

    def WriteSettings(self, filename, name, dictionary):
        try:
            f = open(filename, "wb")
            p = f.write
            p("from mpmath import *" + nl)
            p("from rational import Rational" + nl)
            p("from integer import Zn" + nl)
            p("from julian import Julian" + nl)
            p("mp.dps = " + str(mp.dps) + nl + nl)
            p(name + " = {" + nl)
            keys = dictionary.keys()
            keys.sort()
            indent = "  "
            for key in keys:
                s = repr(dictionary[key])
                if s == "<pi: 3.14159~>": s = "pi"
                p(indent + '"' + key + '"' + " : " + s + "," + nl)
            p("}" + nl)
            f.close()
        except Exception, e:
            msg = ("%sError trying to write dictionary '%s':" % (fln(), name)) + nl + str(e)
            self.display.msg(msg)
            raise

    def PrintRegisters(self):
        """
    Usage: regs

    Displays the contents of the registers
        """
        if not self.registers:
            raise ValueError("%sThere are no registers defined" % fln())
        names = self.registers.keys()
        names.sort()
        lengths = [len(name) for name in names]
        fmt = "%%-%ds  %%s\n" % max(lengths)
        s = ""
        for name in names:
            s += fmt % (name, self.Format(self.registers[name]))
        self.display.msg(s)

    def CheckEnvironment(self):
        '''Look at the environment variables defined in
        cfg["environment"] and execute any commands in them.  Note we only
        do this if we're not using the default configuation (-d option).
        '''
        if self.use_default_config_only: return
        for var in config.cfg["environment"]:
            if var in os.environ:
                try:
                    finished = False
                    status = None
                    cmd_line = os.environ[var]
                    cmds = ParseCommandInput(cmd_line)
                    n = len(cmds) - 1
                    for i, cmd in enumerate(cmds):
                        status = ProcessCommand(cmd, self.commands_dict, i==n)
                        if status == status_quit:
                            finished = True
                    if finished:
                        raise Exception("%sGot a quit command" % fln())
                except Exception, e:
                    msg = "%sFor environment variable '%s', got exception:" + nl
                    self.display.msg(msg % (fln(), var) + str(e))


    def RunChecks(self):
        '''Run checks on various things to flag things that might need to be
        fixed.
        '''
        if not self.run_checks:  return
        # Look for commands that don't have associated help strings
        method = regex.compile(r"<bound method [_a-z][_a-z0-9]*[.]([^ ]*) of .*", regex.I)
        undocumented = []
        for f in self.commands_dict.keys():
            if self.commands_dict[f][0].__doc__ is None:
                name = method.match(self.commands_dict[f][0].__str__()).groups()[0]
                undocumented.append(name)
        if len(undocumented):
            print "undocumented functions: %s" % ' '.join(undocumented)

    def GetRegisterName(self, cmd):
        cmd = strip(cmd)
        if len(cmd) < 2:
            raise ValueError("%sYou must give a register name" % fln())
        return cmd[1:]

    def RecallRegister(self, cmd):
        name = GetRegisterName(cmd)
        if name not in self.registers:
            raise ValueError("%sRegister name '%s' is not defined" % (fln(), name))
        self.stack.push(self.registers[name])
        return status_ok

    def StoreRegister(self, cmd):
        name = GetRegisterName(cmd)
        if len(self.stack) == 0:
            raise Exception("%sStack is empty" % fln())
        self.registers[name] = stack[0]
        return status_ok

    def C_int(self, cmd, val):
        try:
            n = int(val)
        except:
            msg = "%%s'%s' is not a valid integer for int command" % val
            raise ValueError(msg)

        if n > 0:
            if n < 1:
                msg = "%sInteger for int command must be > 0"
                raise ValueError(msg % fln())
            Number.bits = n
            Zn.num_bits = n
        else:
            Number.bits = 0
            Zn.num_bits = 0
        # TODO This is ugly and needs refactoring...
        if cmd == 's':
            Number.signed = True
            Zn.is_signed = True
        else:
            Number.signed = False
            Zn.is_signed = False

    def C_sX(self, val):
        """
    Usage: sX where X is 'X' or X is an integer
           x sX -> set C signed integer mode with bits defined by the
                top value on the stack
           s5 -> set C signed integer mode with 5 bits
        """
        self.C_int('s', val)

    def C_uX(self, val):
        """
    Usage: uX where X is 'X' or X is an integer
           x uX -> set C unsigned integer mode with bits defined by the
                top value on the stack
           u5 -> set C unsigned integer mode with 5 bits
        """
        self.C_int('u', val)

    def GreaterThanEqual(self, x, y):
        """
    Usage: y x >=

    If x >= y, return True; otherwise, return False.
        """
        result = (x >= y)
        if not result and self.testing: exit(1)
        return result

    def GreaterThan(self, x, y):
        """
    Usage: y x >

    If x > y, return True; otherwise, return False.
        """
        result = (x > y)
        if not result and self.testing: exit(1)
        return result

    def LessThanEqual(self, x, y):
        """
    Usage: y x <=

    If x <= y, return True; otherwise, return False.
        """
        result = (x <= y)
        if not result and self.testing: exit(1)
        return result

    def LessThan(self, x, y):
        """
    Usage: y x <

    If x < y, return True; otherwise, return False.
        """
        result = (x < y)
        if not result and self.testing: exit(1)
        return result

    def Equal(self, x, y):
        """
    Usage: y x =

    If x and y are equal, return True; otherwise, return False.
        """
        result = (x == y)
        if not result and self.testing:
            exit(1)
        return result

    def NotEqual(self, x, y):
        """
    Usage: y x !=

    If x and y are not equal, return True; otherwise, return False.
        """
        result = (x != y)
        if not result and self.testing:
            exit(1)
        return result

    def DisplayEqual(self, x, y):
        """
    Usage: y x =

    If x and y are equal, return True; otherwise, return False.
    IMPORTANT:  the test of equality is whether the string
    representations in the current display mode match.  If they do,
    the numbers are equal.  If testing is True, then if the two
    numbers compare false, we exit with a status of 1.
        """
        sx, sy = self.Format(x), self.Format(y)
        result = (sx == sy)
        if not result and self.testing:
            exit(1)
        return result

    def cleanup(self):
        self.SaveConfiguration()
        readline.write_history_file(os.path.expanduser('~')+'/.pycalc/history')

    def push(self, val):
        if val is None:
            raise Exception('%sBAD! val is None'%fln())
        self.stack.push(val)

    def pop(self):
        return self.stack.pop()

    def read_line(self, stream=None):
        if stream:
            s = stream.readline()
        elif self.process_stdin:
            s = sys.stdin.readline()
        else:
            try:
                line = raw_input(config.cfg["prompt"])
            except KeyboardInterrupt:
                print
                sys.exit()
        # it looks like readline automatically adds stuff to history
        #readline.add_history(line)
        if line and line[-1] == nl:
            self.display.log('--> "%s"' % line, suppress_nl=True)
        else:
            self.display.log('--> "%s"' % line)
        pos = line.find("#")  # Delete comments
        if pos != -1:
            line = line[:pos]
            if pos == 0:
                # This line was nothing but a comment
                return ''
        return line

    def chomp(self, line):
        return self.chomppost.sub("", self.chomppre.sub("", line))

    def flatten_tags(self, tags):
        ft = []
        while type(tags) == list and len(tags):
            ft.append(tags[0][0])
            tags = tags[0][3]
        return ft

    def token(self):
        # snag the next token from the line
        line = self.read_line()
        # print "got new line: '%s'" % line
        success, taglist, next = TextTools.tag(line, self.parser)
        while self.chomp(line) != '':
            if not success:
                raise ParseError("Not a command or value: '%s'"%line)
            #print "'%s' yielding '%s', '%s'" % (line, self.chomp(line[:next]),taglist)
            yield self.chomp(line[:next]), self.flatten_tags(taglist), line[next:]
            # print "back from yield '%s'"%self.chomp(line[:next])
            line = self.chomp(line[next:])
            #line = line[next:]
            success, taglist, next = TextTools.tag(line, self.parser)

    def prepare_args(self, fn, inf):
        if debug(): print "prepare_args(%s,%s)"%(fn,str(inf))
        args = []
        v = None
        n = inf[1]
        if len(inf) == 3 and n == 'match' and 'regex' in inf[2]:
            matches = inf[2]['regex'].match(fn)
            if matches:
                return [ g for g in matches.groups() ]
            return []
        if n == 'x':
            v = self.pop()
            if not isint(v):
                self.push(v)
                raise TypeError("'%s' requires an integer as the first argument" % fn)
            n = int(v)
        #print "stack size is %d"%len(self.stack)
        l = len(self.stack)
        if n > l:
            if v is not None:
                self.push(v)
                raise IndexError("'%d %s' requires 1+%d args (stack size is %d)" %
                    (n, fn, n, l+1))
            else:
                raise IndexError("'%s' requires %d args (stack size is %d)" %
                    (fn, n, l))
        for i in range(n):
            val = self.pop()
            args.insert(0, val)
        return args

    def run(self):
        isiterable = lambda obj: getattr(obj, '__iter__', False)
        while True:
            arg = ''
            try:
                for arg,tag,line in self.token():
                    if debug(): print arg,line,tag
                    if 'func' in tag:
                        cmdidx = tag.index('func')
                    else:
                        cmdidx = False
                    if cmdidx and tag[cmdidx+1] in self.commands_dict:
                        command = tag[cmdidx+1]
                    else:
                        command = arg
                    if arg == "const":
                        cv = self.commands_dict['const'][0](line)
                        if cv is not None:
                            self.push(cv)
                        break
                    elif command in self.commands_dict:
                        try:
                            args = self.prepare_args(arg, self.commands_dict[command])
                            if debug(): print args
                            try:
                                retval = self.commands_dict[command][0](*args)
                            except (ValueError, TypeError), e:
                                retval = args
                                if debug():
                                    self.errors.append(traceback.format_exc())
                                else:
                                    self.errors.append(str(e))
                        except (IndexError, TypeError), e:
                            self.errors.append(str(e))
                            continue
                        if not isiterable(retval):
                            retval = [retval]
                        for v in retval:
                            if v is not None:
                                if isint_native(v):
                                    v = Zn(v)
                                self.push(v)
                    elif arg in ['null', 'nop']:
                        pass
                    else:
                        # this should be a number....
                        num = self.chomp(arg)
                        #print "num = '%s', arg = '%s'"%(num,arg)
                        if len(num) > 0:
                            try:
                                num = self.number(self.chomp(arg), tag)
                                if num is not None:
                                    self.push(num)
                            except ValueError:
                                self.errors.append("Invalid input: %s" % arg)
                if arg not in ['help', '?']:
                    self.DisplayStack()
            except EOFError:
                break
            except ParseError:
                type,value,tb = sys.exc_info()
                print value
            except SystemExit:
                raise
            except:
                print "Something bad happened.  Don't do that again!"
                type,value,tb = sys.exc_info()
                traceback.print_exception(type, value, tb, None, sys.stdout)
        readline.write_history_file()

    def help(self, args=None):
        """
    Usage: help [function]

    Lists the functions implemented or displays help for the
    requested function
        """
        if args is not None:
            args = args.split()
            if args:
                arg = args[0]
                if arg in self.commands_dict:
                    if self.commands_dict[arg][0].__doc__ is None:
                        print "No help for %s" % arg
                    else:
                        print self.commands_dict[arg][0].__doc__
                else:
                    print "unknown function:", arg
                return
        maxlen = 0
        functions = []
        for k in self.commands_dict.iterkeys():
            maxlen = max(maxlen, len(k))
            functions.append(k)
        functions.sort();
        maxlen += 1
        printed = 0
        for k in functions:
            s = k + " "*maxlen
            print s[0:maxlen],
            printed += maxlen
            if printed > (72-maxlen):
                print
                printed = 0
        print "\n"
        print "Delimiters are space, tab, and newline.\n"

    def list_constants(self):
        """
    Usage: constants

    Lists the constants available for use by name
        """
        for a,k in self.constants.iteritems():
            print "%s = %s" % (a, k.show(self.base, config.cfg["prec"], self.vector_mode, self.angle_mode))

    def warranty(self):
        """
    Usage: warranty

    Displays the license and warranty information
        """
        global __doc__
        print __doc__

    def todo(self):
        """
    Usage: todo

    Displays the things that still need to be fixed
        """
        print """
parse vector [2 3]
parsing of things that should break
    dup23
        """

    def quit(self):
        """
    Usage: quit

    Exits the program
        """
        sys.exit();

def ParseCommandLine():
    from optparse import OptionParser
    usage = "usage: %prog [options]"
    descr = "Command line RPN calculator"
    parser = OptionParser(usage, description=descr)
    c,d,s,r,t,v = ("Check that commands have help info",
                   "Use default configuration in hc.py file only",
                   "Take input from stdin",
                   "Read input from file",
                   "Exit with status 1 if = or == are False",
                   "Display program version")
    parser.add_option("-c", "--run-checks", action="store_true", help=c)
    parser.add_option("-d", "--default-config", action="store_true", help=d)
    parser.add_option("-s", "--read-stdin", action="store_true", help=s)
    parser.add_option("-r", "--read-file", dest="file", help=r)
    parser.add_option("-t", "--testing-mode", action="store_true", help=t)
    parser.add_option("-v", "--version", action="store_true", help=v)
    return parser.parse_args(args=None, values=None)

def main(argv):
    finished = False
    status = None
    opt, arg = ParseCommandLine()
    console.set_title("Vrnn's Smrt Clcltr")
    calculator = Calculator(arg, opt)
    try:
        calculator.run()
    except KeyboardInterrupt, e:
        pass
    except EOFError, e:
        pass
    print
    sys.exit(0)
