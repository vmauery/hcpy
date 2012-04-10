import os, sys
#---------------------------------------------------------------------------
# Configuration information.

defcfg = {
    # If any of these environment variables exist, execute the
    # commands in them.
    "environment" : ["HCPYINIT", "hcpyinit"],

    # Angle mode:  must be either 'deg' or 'rad'
    "angle_mode" : "rad",

    # Integer mode: must be 'dec', 'hex', 'oct', 'bin', or 'ip'.
    "integer_mode" : "dec",

    # title to display on console for xterms and stuff
    "console_title" : "%(tty)s: hc",

    # Prompt to be displayed for each input (may be empty)
    "prompt" : "> ",

    # If true, coerce means to change arguments' types as needed to
    # calculate a function.  Otherwise, a ValueError exception will be
    # thrown.
    "coerce" : True,

    # If true, we'll allow x/0 to be infinity as long as x != 0
    "allow_divide_by_zero" : True,

    # The calculator is capable of handling rational numbers.  If
    # you'd rather a division of two integers result in a real number,
    # set no_rationals to True.  The rat command toggles this setting.
    "no_rationals" : True,

    # For display of complex numbers.  Mode is either rect or polar.
    "imaginary_mode" : "rect",
    "imaginary_unit" : "i",
    "imaginary_unit_first" : False,  # If true, 1+i3
    "imaginary_space" : False,       # If true, 1 + 3i or (1, 3)
    "ordered_pair" : False,          # If true, (1,3)
    "polar_separator" : " <| ",      # Used in polar display
    "degree_symbol" : "deg",         # chr(248) might be OK for cygwin bash
    "infinity_symbol" : "inf",       # chr(236) might be OK for cygwin bash
    "arg_digits" : 2,                # Separate formatting num digits for args
    "arg_format" : "fix",            # Separate formatting type for arguments

    # Factorials of integers <= this number are calculated exactly.
    # Otherwise, the mpmath factorial function is used which returns either
    # an mpf or mpc.  Set to zero if you want all factorials to be
    # calculated exactly (warning:  this can result in long calculation
    # times and lots of digits being printed).
    "factorial_limit" : 20001,

    # The following string is used to separate commands on the command
    # input.  If this string is not in the command line, the command is
    # parsed into separate commands based on whitespace.
    "command_separator" : ";",

    # How many items of the stack to show.  Use 0 for all.
    "stack_display" : 0,

    # If the following variable is True, we will persist our settings from
    # run to run.  Otherwise, our configuration comes from this dictionary
    # and the stack and registers are empty when starting.
    "persist" : False,
    "persist_registers" : False,
    "persist_stack" : False,

    # The following variables determines how floating point numbers are
    # formatted. Legitimate values are:  fix for fixed, sig for significant
    # figures, sci for scientific, eng for engineering, engsi for
    # engineering with SI prefixes after the number, and "none".  If
    # "none", then the default mpmath string representation is used.
    # Change the mpformat.py file if you wish to change things such as
    # decimal point type, comma decoration, etc.
    #
    # NOTE:  there are mpFormat class variables in mpformat.py that you
    # might want to examine and set to your tastes.  They are not included
    # here because they will probably be set only once by a user.
    "fp_format" : "raw",
    "fp_digits" : 10,
    "fp_show_plus_sign" : False,  # If true, "+3.4" instead of " 3.4"
    "fp_comma_decorate" : True,   # If true, 1,234 instead of 1234
    "fp_cuddle_si" : False,       # If true, "12.3k" instead of "12.3 k"

    # Set how many digits of precision the mpmath library should use.
    "prec" : 30,

    # Settings for interval number display.  iva mode is a+-b; ivb mode is
    # a(b%); ivc mode is <a, b>.
    "iv_space" : True,              # a +- b vs. a+-b, a (b%) vs a(b%)
    "iv_mode" : "b",

    # If brief is set to true, truncate numbers to fit on one line.
    "brief" : True,

    # Set the line width for the display.  Set it to a negative integer to
    # instruct the program to try to first read the value from the COLUMNS
    # environment variable; if not present, it will use the absolute value.
    "line_width" : 75,

    # String to use when an ellipsis is needed (used by brief command)
    "ellipsis" : "."*3,

    # If true, display fractions as mixed fractions.
    "mixed_fractions" : True,

    # If this number is not 1, then it is used for modular arithmetic.
    "modulus" : 1,

    # If this is true, then results are downcast to simpler types.
    # Warning:  this may cause suprises.  For example, if prec is set to 15
    # and you execute 'pi 1e100 *', you'll see an integer result.  This
    # isn't a bug, but is a consequence of the finite number of digits of
    # precision -- the downcast happens because int(1e100*pi) == 1e100*pi.
    "downcasting" : False,

    # The integers used by this program exhibit python-style division.
    # This behavior can be surprising to a C programmer, because in
    # python, (-3) // 8 is equal to -1 (this is called floor division);
    # most C programmers would expect this to be zero.  Set the
    # following variable to True if you want (-3) // 8 to be zero.
    "C_division" : True,

    # Scripts that can be called using the ! command are in the following
    # directory.  Any python script in this directory will have its main()
    # function called and the value that is returned will be pushed on the
    # stack.  This lets you write auxiliary scripts that prompt you to help
    # you get a number you need without cluttering up the commands or
    # registers of this program.  Example:  an astronomy.py script could
    # prompt you for which astronomical constant you wanted to use.  Set
    # this entry to the empty string or None if you don't want this
    # behavior.  You can also give the name of the function you want to
    # be called.  This function will be called with the display object,
    # which you can use to send messages to the user.
    "helper_scripts" : "d:/p/math/hcpy/helpers",
    "helper_script_function_name" : "main",
}
cfg = {}

def load():
    global cfg
    global defcfg
    cfg = defcfg
    config_file = os.path.expanduser(os.path.join("~", ".config", "hc", "config"))
    try:
        os.stat(config_file)
    except OSError:
        return
    try:
        execfile(config_file, {}, cfg)
    except Exception, e:
        print "Error loading user config: %s" % str(e)
    for k in cfg:
        if k not in defcfg:
            del cfg[k]


