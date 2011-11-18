'''
Provide numeric objects that interoperate with mpmath objects.

---------------------------------------------------------------------------
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
'''

from mpmath import mpf, mpc, mpi, ctx_iv, eps, mp, pi
from mpformat import mpFormat, inf
from debug import *
import socket
import time
import re
from string import strip
from si import suffixes_ln

try: from pdb import xx  # pdb.set_trace is xx; easy to find for debugging
except: pass

def isint_native(x):
    return isinstance(x, (int, long))

def isint(x):
    return isinstance(x, (int, long, Zn))

def gcd(a, b):
    '''Determine the greatest common divisor of integers u and v.
    Euclid's algorithm from Knuth, vol 2, pg 320.
    '''
    if not isint(a) or not isint(b):
        raise ValueError("Arguments must be integers")
    a, b = int(a), int(b)
    if a == 0:  return b
    if b == 0:  return a
    while b != 0:
        a, b = b, a % b
    return a

class Rational(object):
    mixed = False  # If set to true, str() returns mixed form
    def __init__(self, a=0, b=1):
        if b == 1:
            if isinstance(a, mpf):
                r = self.frac(a)
                self.n, self.d = r.n, r.d
                return
            else:
                try:
                    # C++ copy constructor behavior
                    self.n = a.numer()
                    self.d = a.denom()
                    return
                except:
                    pass
        if b == 0:
            raise ZeroDivisionError("Denominator is zero")
        else:
            g = gcd(a, b)
            self.n = int(a)//g
            self.d = int(b)//g

    def __abs__(self):
        return Rational(abs(self.n), abs(self.d))

    def __pos__(self):
        return Rational(self)

    def __neg__(self):
        return Rational(-self.n, self.d)

    def __radd__(self, other):
        return self.__add__(other)

    def __add__(self, other):
        if isinstance(other, Rational):
            y = Rational(self.n*other.d + other.n*self.d, self.d*other.d)
            if y.d == 1:
                return y.n
            else:
                return y
        elif isinstance(other, float):
            raise ValueError("float addition not supported")
        else:
            assert isinstance(other, int) or isinstance(other, mpf) or \
                   isinstance(other, mpc) or isinstance(other, ctx_iv.ivmpf)
            n = other*self.d + self.n
            if isinstance(n, int):
                y = Rational(n, self.d)
                if y.d == 1:
                    return y.n
                else:
                    return y
            else:
                return n/self.d

    def __rsub__(self, other):
        return -self.__sub__(other)

    def __sub__(self, other):
        if isinstance(other, Rational):
            return Rational(self.n*other.d - other.n*self.d, self.d*other.d)
        elif isinstance(other, float):
            raise ValueError("float subtraction not supported")
        else:
            assert isinstance(other, int) or isinstance(other, mpf) or \
                   isinstance(other, mpc) or isinstance(other, ctx_iv.ivmpf)
            n = self.n - other*self.d
            if isinstance(n, int):
                y = Rational(n, self.d)
                if y.d == 1:
                    return y.n
                else:
                    return y
            else:
                return n/self.d

    def __rmul__(self, other):
        return self.__mul__(other)

    def __mul__(self, other):
        if isinstance(other, Rational):
            return Rational(self.n*other.n, self.d*other.d)
        elif isinstance(other, float):
            raise ValueError("float multiplication not supported")
        else:
            assert isint(other) or isinstance(other, mpf) or \
                   isinstance(other, mpc) or isinstance(other, ctx_iv.ivmpf)
            n = other*self.n
            if isinstance(n, int):
                y = Rational(n, self.d)
                if y.d == 1:
                    return y.n
                else:
                    return y
            else:
                return n/self.d

    def __rdiv__(self, other):
        one = Rational(1)
        return one/(self*(one/other))

    def __div__(self, other):
        if other == 0:
            raise ZeroDivisionError("Divisor is zero")
        if isinstance(other, Rational):
            return Rational(self.n*other.d, self.d*other.n)
        elif isint(other):
            return Rational(self.n, self.d*int(other))
        elif isinstance(other, float):
            raise ValueError("float division not supported")
        else:
            assert isinstance(other, mpf) or \
                   isinstance(other, mpc) or \
                   isinstance(other, ctx_iv.ivmpf)
            return (self.n/other)/self.d

    def __truediv__(self, other):
        return self.__div__(other)

    def _mixed(self):
        sign = ""
        if (self.n < 0 and self.d > 0) or (self.n > 0 and self.d < 0):
            sign = "-"
        whole, remainder = divmod(abs(self.n), abs(self.d))
        if  whole == 0:
            if remainder == 0:
                return "0"
            else:
                return(sign + "%s/%s" % (remainder, self.d))
        else:
            if remainder == 0:
                return sign + str(whole)
            else:
                return(sign + "%s %s/%s" % (whole, remainder, self.d))

    def __str__(self):
        if Rational.mixed:
            return self._mixed()
        else:
            return "%d/%d" % (self.n, self.d)

    def __repr__(self):
        return "Rational(%d, %d)" % (self.n, self.d)

    def __float__(self):
        return float(self.n)/self.d

    def __getattr__(self, key):
        '''This function is necessary when a Rational gets compared to an
        mpf.  The mpf comparison routine looks for the _mpf_ attribute; if
        it finds it, then it knows it has an mpf.  We fake it out by doing
        an mpf conversion at the point the information is needed.  Thus,
        the comparison will be done with the proper number of digits; this
        wouldn't necessarily be true if we cached the _mpf_ data earlier.
        '''
        if key == "_mpf_":
            return (mpf(self.n)/mpf(self.d))._mpf_
        else:
            raise AttributeError("'%s' not an attribute" % key)

    def __cmp__(self, other):
        if other == None:
            return -1
        if isinstance(other, Rational):
            if (self.n == other.n) and (self.d == other.d):
                return 0
            else:
                if self.mpf() < other.mpf(): return -1
                else: return 1
        elif isinstance(other, mpf):
            a, b = self.mpf(), other.mpf()
            if a < b: return -1
            elif a > b: return 1
            else: return 0
        elif isinstance(other, int):
            a, b = self.mpf(), mpf(other)
            if a < b: return -1
            elif a > b: return 1
            else: return 0
        else:
            raise ValueError("Second argument is unsupported type")

    def mpf(self):
        return mpf(self.n)/mpf(self.d)

    def mpc(self):
        return mpc(mpf(self.n)/mpf(self.d), 0)

    def mpi(self):
        n = mpf(self.n)/mpf(self.d)
        return mpi(n, n)

    def denom(self):
        return self.d

    def numer(self):
        return self.n

    def frac(self, x, digits=0, max_iterations=0):
        '''Converts an mpf to a Rational approximation and returns a
        Rational object.  The iteration goal is to return a rational
        approximation whose difference from x is less than 1/10**mp.dps;
        i.e., the approximation as good as the current number of digits
        used in the mpmath library.

        x will be converted to an mpf type.  This routine is converted from
        a C routine I downloaded from NETLIB on Tue Nov 11 18:08:26 1997
        from NIST Guide to Available Math Software, source for module FRAC
        from package C.

        digits sets the precision of the conversion at 10**(-digits).  If
        it is 0, then mp.dps is used.

        Set max_iterations to a postive nonzero value to limit the number
        of iterations.

        See http://gams.nist.gov/,
        http://gams.nist.gov/serve.cgi/Package/C/,
        http://www.netlib.org/c/frac (location of actual source)
        '''
        if isinstance(x, mpc) or isinstance(x, complex):
            x = abs(x)
        elif isinstance(x, ctx_iv.ivmpf):
            x = x.mid
        elif isinstance(x, int) or isinstance(x, long):
            x = mpf(x)
        elif isinstance(x, Zn):
            x = mpf(int(x))
        else:
            if not isinstance(x, mpf):
                # Note we explicitly do not handle floats
                raise SyntaxError("Unsupported type")
        # Handle a special case that results in no convergence.
        if x == mpf("1")/mpf("2"):
            return Rational(1, 2)
        digits_increase = 2
        if digits == 0:
            digits = mp.dps
        error = mpf(10)**(-digits)
        # We'll increase the precision to ensure the error is greater
        # than the minimal fractional digits.
        mp.dps += digits_increase
        sign = 1
        if x < 0:
            sign = -1
            x = abs(x)
        if x == 0:
            mp.dps -= digits_increase
            return Rational(0, 1)
        d, D, n, r, epsilon, iterations = 1, 1, int(x), mpf(1), mpf(0), 0
        N = n + 1
        one = True
        while True:
            iterations += 1
            if max_iterations and iterations > max_iterations:
                raise Exception("Too many iterations")
            if not one:
                if r > 1.0:
                    N += n*int(r)
                    D += d*int(r)
                    n += N
                    d += D
                else:
                    r = 1.0/r
            else:
                one = False
            r = 0.0
            if x*d != n:
                r = (N - x*D)/(x*d - n)
                if r <= 1.0:
                    t = N
                    N = n
                    n = t
                    t = D
                    D = d
                    d = t
            epsilon = abs(1.0 - n/(x*d))
            if epsilon <= error:
                mp.dps -= digits_increase
                return Rational(sign*n, d)
            if r != 0.0:
                continue
            m = 10.0
            while m*epsilon < 1.0:
                m *= 10.0
            epsilon = (1.0/m)*(int(0.5 + m*epsilon))
            if r == 0.0:
                mp.dps -= digits_increase
                raise Exception("Error in frac() routine")

if __name__ == "__main__":
    # Test code

    def gcd_tests():
        assert(gcd(8, 12) == 4)
        assert(gcd(8L, 12) == 4)
        assert(gcd(8, 12L) == 4)
        assert(gcd(8L, 12L) == 4)

    def generalTests():
        three = Rational(3)
        assert(str(three) == "3/1")
        third = Rational(1,3)
        assert(str(third) == "1/3")
        fifth = Rational(1,5)
        assert(str(fifth) == "1/5")
        assert(str(third + fifth) == "8/15")
        assert(str(third*fifth) == "1/15")
        assert(str(third-fifth) == "2/15")
        assert(str(fifth/third) == "3/5")
        assert(float(fifth) == 0.2)
        assert(float(third) == 1/3.)
        assert(third.numer() == 1)
        assert(third.denom() == 3)
        new_third = Rational(third)  # Show we support init from ourselves
        assert(str(new_third) == "1/3")

    def mixedTests():
        badPi = Rational(22,7)
        Rational.mixed = True
        assert(str(badPi) == "3 1/7")
        proper_fraction = Rational(3,5)
        assert(str(proper_fraction) == "3/5")
        whole_num  =  Rational(8,2)
        assert(str(whole_num) == "4")
        zero  =  Rational(0,1)
        assert(str(zero) == "0")

    def errorTests():
        try:
            r = Rational(5, 0)
            raise Exception("Fail: didn't detect zero denominator.")
        except ZeroDivisionError, detail:
            pass

    def no_diff(a, b):
        return abs(a) - abs(b) <= eps

    def mpfTests(): # Test with mpf's
        onep2 = mpf("1.2")
        third = Rational(1,3)
        # add
        assert no_diff(onep2 + third, onep2 + 1/mpf(3))
        # radd
        assert no_diff(third + onep2, onep2 + 1/mpf(3))
        # sub
        assert no_diff(onep2 - third, onep2 - 1/mpf(3))
        # rsub
        assert no_diff(third - onep2, -onep2 + 1/mpf(3))
        # mul
        assert no_diff(onep2 * third, onep2 * 1/mpf(3))
        # rmul
        assert no_diff(third * onep2, onep2 * 1/mpf(3))
        # div
        assert no_diff(onep2 / third, onep2 /(1/mpf(3)))
        # rdiv
        assert no_diff(third / onep2, (1/mpf(3)) / onep2)

    def mpcTests(): # Test with mpc's
        onep2 = mpc("1.2", "88")
        third = Rational(1,3)
        # add
        assert no_diff(onep2 + third, onep2 + 1/mpc(3))
        # radd
        assert no_diff(third + onep2, onep2 + 1/mpc(3))
        # sub
        assert no_diff(onep2 - third, onep2 - 1/mpc(3))
        # rsub
        assert no_diff(third - onep2, -onep2 + 1/mpc(3))
        # mul
        assert no_diff(onep2 * third, onep2 * 1/mpc(3))
        # rmul
        assert no_diff(third * onep2, onep2 * 1/mpc(3))
        # div
        assert no_diff(onep2 / third, onep2 /(1/mpc(3)))
        # rdiv
        assert no_diff(third / onep2, (1/mpc(3)) / onep2)

    def no_idiff(a, b):
        if not no_diff(a.mid, b.mid):
            return False
        if not no_diff(a.delta, b.delta):
            return False
        return True

    def mpiTests(): # Test with mpi's
        onep2 = mpi("1", "2")
        third = Rational(1,3)
        # add
        assert no_idiff(onep2 + third, onep2 + 1/mpi(3))
        # radd
        assert no_idiff(third + onep2, onep2 + 1/mpi(3))
        # sub
        assert no_idiff(onep2 - third, onep2 - 1/mpi(3))
        # rsub
        assert no_idiff(third - onep2, -onep2 + 1/mpi(3))
        # mul
        assert no_idiff(onep2 * third, onep2 * 1/mpi(3))
        # rmul
        assert no_idiff(third * onep2, onep2 * 1/mpi(3))
        # div
        assert no_idiff(onep2 / third, onep2 /(1/mpi(3)))
        # rdiv
        assert no_idiff(third / onep2, (1/mpi(3)) / onep2)

    def integerTests():
        one = 1
        third = Rational(1,3)
        # add
        assert (one + third) == (Rational(4,3))
        # radd
        assert (third + one) == (Rational(4,3))
        # sub
        assert (one - third) == (Rational(2,3))
        # rsub
        assert (third - one) == (Rational(-2,3))
        # mul
        assert (one * third) == (Rational(1,3))
        # rmul
        assert (third * one) == (Rational(1,3))
        # div
        assert (one / third) == (Rational(3,1))
        # rdiv
        assert (third / one) == (Rational(1,3))

    def comparisonTest():
        one = 1
        oned = mpf(1)
        third = Rational(1,3)
        # < and <=
        assert third < one
        assert third < oned
        assert third <= one
        assert third <= oned
        assert third <= third
        assert third <= third
        # > and >=
        assert one   >  third
        assert oned  >  third
        assert one   >= third
        assert oned  >= third
        assert third >= third
        assert third >= third
        # == and !=
        assert third != one
        assert one   != third
        assert third != oned
        assert oned  != third
        assert third == Rational(1, 3)
        assert not (third != Rational(1, 3))

    def fracTest():
        '''I experimentally determined that factor = 3.1 for the loop
        shown up to 1000 digits (but it only goes to 100 for speed).
        '''
        factor = mpf(3.1)
        for digits in xrange(5, 100, 5):
            mp.dps = digits
            r = Rational()
            approx = r.frac(pi)
            num = mpf(approx.n)/approx.d
            #assert abs(pi - num) < factor*mpf(10)**(-digits)
            if abs(pi - num) > factor*mpf(10)**(-digits):
                print "digits", digits,
                a = str(abs(pi - num))
                print a[:6], a[-6:]
    generalTests()
    mixedTests()
    errorTests()
    mpfTests()
    mpcTests()
    mpiTests()
    comparisonTest()
    fracTest()
    exit(0)
'''
$Id: integer.py 1.22 2009/02/10 05:24:01 donp Exp $

Zn objects can be set to behave as a signed or unsigned n n-bit integer
in addition to a regular python integer.  The properties you can set
are:

    bits        Must be >= 0.  If 0, then the object behaves like a
                regular python integer.  If > 0, then the object
                behaves like an integer with that number of bits.

    signed      If true, behaves as a signed integer; otherwise
                behaves as an unsigned integer.  Since the normal unlimited
                precision python integers cannot have an unsigned
                characteristic, this only applies when bits is
                nonzero.

    C_division  If True, integer division behaves like it does in C;
                i.e., 3/8 == -3/8 == 3/-8 == 0.  If False, then the
                behavior is pythonic:  floor division.  This means
                -3/8 == 3/-8 == -1.

Four bit two's complement representation:
Uns binary 2's comp
15  1111   -1
14  1110   -2
13  1101   -3
12  1100   -4
11  1011   -5
10  1010   -6
 9  1001   -7
 8  1000   -8
 7  0111    7
 6  0110    6
 5  0101    5
 4  0100    4
 3  0011    3
 2  0010    2
 1  0001    1
 0  0000    0
'''


def isint(x):
    return isinstance(x, int) or isinstance(x, long) or isinstance(x, Zn)

class Zn(object):
    # These characters are used in the str representation of Zn objects
    # Example:  a 4-bit signed value of -2 is given as '-2<4s>'.
    left  = "<"
    right = ">"
    space = ""  # Put a space between the number and its designator

    # These variables are private
    num_bits = 0
    is_signed = True
    use_C_division = False
    base = 0
    high_unsigned = 0
    high_signed = 0

    # This variable is used to hold 0, 1, or 2.  These settings have
    # to do with the subtleties of negating 2's complement numbers.
    # See the comments in the __neg__ method.
    negate_zero = 0

    def __init__(self, value=0, proto=None):
        self.n = 0
        if proto is not None:
            self.num_bits = proto.num_bits
            self.is_signed = proto.is_signed
        elif isinstance(value, Zn):
            self.num_bits = value.num_bits
            self.is_signed = value.is_signed
        else:
            self.num_bits = Zn.num_bits
            self.is_signed = Zn.is_signed
        if self.num_bits > 0:
            self.base = 2**self.num_bits
        if type(value) == str:
            self._set_from_string(value)
        elif isinstance(value, int) or isinstance(value, long):
            self.value = value
        elif isinstance(value, Zn):
            self.value = value.value
        else:
            raise TypeError("%sCan't set integer from value '%s'" % \
                (fln(), str(value)))
        self._update()

    def _set_from_string(self, value):
        '''We set our value from a string.  This can be either a
        regular string for an integer or long or a string gotten from
        our str() method.  In the second case, the value will be made
        to fit in the current representation, regardless of how many
        bits or whether it was signed or unsigned when str'd.  Note
        that the left character must match our current setting or an
        exception will be raised.
        '''
        assert type(value) == type("")
        try:
            if Zn.left in value:
                self.value = int(value.split(Zn.left)[0])
            else:
                self.value = int(value)
        except:
            msg = "%sCan't set integer from '%s'"
            raise ValueError(msg % (fln(), value))


    # Properties
    def get_C_division(self):
        return Zn.use_C_division

    def set_C_division(self, use_C_division):
        if use_C_division != True and use_C_division != False:
            msg = "%suse_C_division class variable must be True or False"
            raise ValueError(msg % fln())
        Zn.use_C_division = use_C_division

    C_division = property(get_C_division, set_C_division, \
        doc="Set to True for C type integer division")

    def get_bits(self):
        return self.num_bits

    def set_bits(self, bits):
        if not isinstance(bits, int) and \
           not isinstance(bits, long) and \
           not isinstance(bits, Zn):
            msg = "%sNumber of bits must be an integer"
            raise ValueError(msg % fln())
        if bits < 0:
            msg = "%sNumber of bits in integer must be >= 0"
            raise ValueError(msg % fln())
        if bits != self.num_bits:
            self.num_bits = bits
            self.num_bits = bits
            if bits == 0:
                self.base = 0 # Will cause 0 div if we use for % calc
            else:
                self.base = 2**bits
            self._update()

    bits = property(get_bits, set_bits, \
        doc="Number of bits in integer (0 for unlimited)")

    def get_signed(self):
        return self.is_signed

    def set_signed(self, signed):
        if signed != True and signed != False:
            raise ValueError("%ssigned must be True or False" % fln())
        self.is_signed = signed
        self._update()

    signed = property(get_signed, set_signed, doc="Signed if True")

    def set_value(self, value):
        if isinstance(value, int) or isinstance(value, long):
            self.n = value
        elif isinstance(value, Zn):
            self.n = value.n
        else:
            try:
                self.n = int(value)
            except:
                raise TypeError("%sCan't set integer from '%s'" % \
                    (fln(), str(value)))
        self._update()

    def get_value(self):
        return self.n

    value = property(get_value, set_value, doc="Integer's value")

    def set_negate_zero(self, value):
        msg = "%svalue must be 0, 1, or 2"
        if not isinstance(value, int) and not isinstance(value, long):
            raise ValueError(msg % fln())
        if value < 0 or value > 2:
            raise ValueError(msg % fln())
        Zn.negate_zero = value
        self._update()

    def get_negate_zero(self):
        return Zn.negate_zero

    negative_zero = property(get_negate_zero, set_negate_zero, \
        doc="If true, -Zn(0) == Zn(-(2**(n-1)))")

    def _update(self):
        "The object's state has changed."
        if self.num_bits != 0:
            assert self.base == 2**self.num_bits
        if self.num_bits == 0:
            self.is_signed = True
        else:
            if self.is_signed:
                self.n &= (self.base - 1)  # Mask off the desired bits
                # If high bit is on, then convert to negative in 2's
                # complement
                if self.n & 2**(self.num_bits - 1):
                    self.n -= self.base
            else:
                self.n &= (self.base-1)
        # Check our invariants
        if self.num_bits == 0:
            assert self.is_signed == True
        else:
            if self.is_signed:
                assert -(self.base >> 1) <= self.n < (self.base >> 1)
            else:
                assert 0 <= self.n < self.base

    def _auto_cast(self, y):
        '''y must be a Zn object for us to interoperate with.  We can
        convert regular integers.
        '''
        if isinstance(y, Zn):
            p = Zn()
            p.is_signed = y.is_signed and self.is_signed
            p.num_bits = max(y.num_bits, self.num_bits)
            y1 = Zn(y, proto=p)
            x1 = Zn(self, proto=p)
            return y1, x1
        elif isinstance(y, mpf):
            return y, mpf(self.value)
        elif isinstance(y, mpc):
            return y, mpc(self.value, 0)
        elif isinstance(y, Julian):
            return y, Julian(self.value)
        elif isinstance(y, ctx_iv.ivmpf):
            return y, mpi(self.value)
        elif isinstance(y, Rational):
            return y, Rational(self.value, 1)
        else:
            return y, self.value

    def __hex__(self):
        self._update()
        t = ""
        if self.num_bits != 0:
            t = self._suffix()
        sign = ""
        num_hex_digits, r = divmod(self.num_bits, 4)
        if r != 0: num_hex_digits += 1
        v = self.n
        if v < 0:
            sign = " "
            if self.is_signed and self.num_bits != 0:
                v = self.n
                v &= (self.base - 1)  # Mask off the desired bits
                v |= (2**(self.bits - 1))
        s = hex(v)[2:]
        if s[-1] == "L": s = s[:-1]     # Remove "L"
        if self.num_bits != 0:
            while len(s) < num_hex_digits:
                s = "0" + s
        if self.num_bits != 0:  assert len(s) == num_hex_digits
        return "%s0x%s%s" % (sign, s, t)

    def __oct__(self):
        self._update()
        t = ""
        if self.num_bits != 0:
            t = self._suffix()
        sign = ""
        num_oct_digits, r = divmod(self.num_bits, 3)
        if r != 0: num_oct_digits += 1
        v = self.n
        if v < 0:
            sign = " "
            if self.is_signed and self.num_bits != 0:
                v = self.n
                v &= (self.base - 1)  # Mask off the desired bits
                v |= (2**(self.bits - 1))
        s = oct(v)[1:]
        if s[0] == "o":  s = s[1:]  # Remove leading 'o' if present
        if self.num_bits != 0:
            while len(s) < num_oct_digits:
                s = "0" + s
        if self.num_bits != 0:  assert len(s) == num_oct_digits
        return "%s0o%s%s" % (sign, s, t)

    def bin(self):
        'Binary representation'
        self._update()
        t = ""
        if self.num_bits != 0:
            t = self._suffix()
        sign = ""
        v = self.n
        if v < 0:
            sign = " "
            if self.is_signed and self.num_bits != 0:
                v = self.n
                v &= (self.base - 1)  # Mask off the desired bits
                v |= (2**(self.bits - 1))
        s = bin(v)[2:]
        while len(s) > 1 and s[0] == "0":  # Remove leading 0's
            s = s[1:]
        if s[-1] == "L": s = s[:-1]     # Remove "L"
        if self.num_bits != 0:
            while len(s) > self.num_bits:  # Trim leading 0's to get num bits
                assert s[0] == "0", "s = '%s'" % s
                s = s[1:]
        if self.num_bits == 0:
            while len(s) > 1 and s[0] == "0":  # Remove leading 0's
                s = s[1:]
        else:
            if len(s) > self.num_bits:
                s = s[len(s) - self.num_bits:]
            else:
                # Add leading zeros if length is not == num bits
                while len(s) < self.num_bits:
                    s = "0" + s
        if self.num_bits != 0:
            assert len(s) == self.num_bits, "s='%s'  %d bits" % (s, self.num_bits)
        return "%s0b%s%s" % (sign, s, t)

    def roman(self):
        self._update()
        sign = " "
        v = self.n
        if v < 0:
            sign = "-"
        v = abs(v)
        s = []
        if v > 3999:
            return str(self)
        while v >= 1000:
            s.append('M')
            v -= 1000
        if v >= 900:
            s.append('CM')
            v -= 900
        if v >= 500:
            s.append('D')
            v -= 500
        if v >= 400:
            s.append('CD')
            v -= 400
        while v >= 100:
            s.append('C')
            v -= 100
        if v >= 90:
            s.append('XC')
            v -= 90
        if v >= 50:
            s.append('L')
            v -= 50
        if v >= 40:
            s.append('XL')
            v -= 40
        while v >= 10:
            s.append('X')
            v -= 10
        if v >= 9:
            s.append('IX')
            v -= 9
        if v >= 5:
            s.append('V')
            v -= 5
        if v >= 4:
            s.append('IV')
            v -= 4
        while v >= 1:
            s.append('I')
            v -= 1
        return "%s%s" % (sign, ''.join(s))

    def __int__(self):
        return self.n

    def __long__(self):
        return self.n

    def __float__(self):
        return mpf(self.n)

    def __complex__(self):
        return mpc(self.n, 0)

    def _suffix(self):
        fmt = Zn.space + Zn.left + "%s%d" + Zn.right
        if self.is_signed:
            return fmt % ("s", self.num_bits)
        else:
            return fmt % ("u", self.num_bits)

    def __str__(self):
        self._update()
        if self.num_bits == 0:
            s = str(self.value)
        else:
            t = self._suffix()
            if self.is_signed:
                s = str(self.value) + t
            else:
                if self.n < 0:
                    s = str(self.base + self.n)
                else:
                    s = str(self.value)
                s += t
        return s

    def __repr__(self):
        if debug():
            if self.signed:
                s = 's'
            else:
                s = 'u'
            return "Zn(%d<%c%d>)"%(self.n, s, self.bits)
        return "Zn(%d)" % self.n

    def _sgn(self, x):
        if x >= 0: return 1
        return -1

    def __abs__(self):
        'See comments under __neg__ for some subleties.'
        msg = "%sCan't take the absolute value of the most negative number"
        if self.is_signed == True and self.n == -(self.base >> 1):
            raise ValueError(msg % fln())
        return Zn(abs(self.n))

    def __neg__(self):
        '''Dealing with the subtleties of 2's complement arithmetic.
        For further details, see
        http://en.wikipedia.org/wiki/Two%27s_complement_arithmetic.

        The negate_zero flag allows three different behaviors with
        regards to negating zero.  I'll use 4-bit integer arithmetic
        to illustrate things.

            0   Normal 2's complement arithmetic.  -Zn(-8) == Zn(-8).
            1   -Zn(-8) == Zn(0), -Zn(0) == Zn(0).
            2   -Zn(-8) == Zn(0), -Zn(0) == Zn(-8).

        In 2's complement arithmetic, the negative of a number is
        gotten by inverting the bits and adding one.  This is handy
        from a hardware perspective because it's fast, but it does
        leave us with a "special" number:  the most negative number.
        It's special because it is its own negative.  This causes a
        problem if you wish to calculate its absolute value, since it
        doesn't exist (in an n-bit signed 2's complement system, the
        values represented go from -2**n to (2**n)-1; 2**n, the
        absolute value of -2**n, can't be represented.  Thus, the only
        logical way to handle abs(-2**n) is to raise an exception.

        Here are the bit representations:
             U        S           S   U
             7  0111  7    1111  -1  15
             6  0110  6    1110  -2  14
             5  0101  5    1101  -3  13
             4  0100  4    1100  -4  12
             3  0011  3    1011  -5  11
             2  0010  2    1010  -6  10
             1  0001  1    1001  -7   9
             0  0000  0    1000  -8   8
        The S are the signed numbers and the U are the unsigned
        numbers.

        The other two situations are provided should the user prefer
        to see that behavior.  All three behaviors are a bit
        counterintuitive, but negate_zero == 1 is probably closest to
        common experience in that -0 == 0.  negate_zero == 2 is
        perhaps logically more satisfying in that all bit patterns are
        now paired together and switch between each other on negation.

        This behavior is now obvious to me, but it wasn't obvious when
        I first wrote __neg__; I simply returned Zn(-self.value).  I
        was playing around with 64 bit integers.  In the calculator, I
        set int64, then entered -1 into x.  I changed to uint64 and
        this changed to 18446744073709551615<u64>.  I divided by 2 and
        cast it back to an I to get 9223372036854775808<u64>.
        Changing this to an int64 yielded -9223372036854775808<s64>.
        This is when I noticed chs was idempotent on this number.  At
        first I thought it was a bug and marked it down to
        investigate.  I reproduced the behavior on 4-bit integers and
        then understood what was going on.  I wish I had read the
        wikipedia article first, as it would have saved me a bunch of
        time.  :^)  For your final lesson in 2's complement arithmetic,
        predict what -9223372036854775808<s64> squared is.

        This also indicates a lesson I learned in a software course:
        always test at the boundaries.  I intuitively knew to mess
        around with numbers like 2**(n-1) and flip between the signed
        an unsigned values and see what happened.
        '''
        if self.n == 0 and Zn.negate_zero and self.is_signed:
            return Zn(-(self.base >> 1))
        return Zn(-self.n)

    def __coerce__(self, other):
        # hmmmmm.... the types we can encounter
        print "coersion: %s, %s" %(self, other)
        if isint_native(other):
            return self, Zn(other)
        ot = type(other)
        if ot == mpc:
            return mpc(self.value, 0), other
        if ot == mpf:
            return mpf(self.value), other
        if isint_native(other):
            return self, Zn(other)
        print "Zn could not coerce type '%s'"%str(ot)
        return None

    def __add__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value + y1.value, proto=x1)
        return x1 + y1

    def __iadd__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value += z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 + x1

    def __radd__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(y1.value + x1.value, proto=x1)
        return y1 + x1

    def __sub__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value - y1.value, proto=x1)
        return x1 - y1

    def __isub__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value -= z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 - x1

    def __rsub__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(y1.value - x1.value, proto=x1)
        return y1 - x1

    def __mul__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value * y1.value, proto=x1)
        return x1 * y1

    def __imul__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value *= z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 * x1

    def __rmul__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value * y1.value, proto=x1)
        return x1 * y1

    def __div__(self, y):
        #print "__div__(%d, %d)"%(self.value, y.value)
        y1, x1 = self._auto_cast(y)
        if Zn.use_C_division:
            if x1.is_signed == True:
                sign = x1._sgn(x1.n)*x1._sgn(y1.n)
                if x1.num_bits != 0:
                    assert -(x1.base >> 1) <= x1.n < (x1.base >> 1)
                    assert -(x1.base >> 1) <= y1.n    < (x1.base >> 1)
                    m = x1.base >> 1
                    return Zn(sign*((abs(x1.value) % m)//(abs(y1.value) % m)), proto=x1)
                else:
                    return Zn(sign*(abs(x1.value)//abs(y1.value)), proto=x1)
            else:
                if x1.num_bits != 0:
                    assert 0 <= x1.n < x1.base
                    assert 0 <= y1.n    < x1.base
                return Zn(x1.n//y1.n, proto=x1)
        else:
            return Zn(x1.value//y1.value, proto=x1)

    def __idiv__(self, y):
        y1, x1 = self._auto_cast(y)
        if not isinstance(y1, Zn) or not isinstance(x1, Zn):
            raise TypeError("Invalid types to integer __idiv__")
        if Zn.use_C_division:
            sign_x = x1._sgn(x1.value)
            sign_y = x1._sgn(y1.value)
            x1.value = sign_x*sign_y*(abs(x1.value)//abs(y1.value))
        else:
            x1.value //= y1.value
        self.is_signed = x1.is_signed
        self.num_bits = x1.num_bits
        self.value = x1.value
        return self

    def __rdiv__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return y1.__div__(x1)
        return y1 / x1

    __floordiv__ = __div__
    __ifloordiv__ = __idiv__
    __rfloordiv__ = __rdiv__

    def __mod__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value % y1.value, proto=x1)
        return x1 % y1

    def __imod__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value %= z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 % x1

    def __rmod__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return y1.__mod__(x1)
        return y1 % x1

    def __cmp__(self, y):
        if isinstance(y, Zn):
            if self.value < y.value:  return -1
            if self.value == y.value: return 0
            else:                     return 1
        elif isinstance(y, int) or isinstance(y, long):
            if self.value < y:  return -1
            if self.value == y: return 0
            else:               return 1
        elif isinstance(y, mpf):
            if mpf(self.value) < y:  return -1
            if mpf(self.value) == y: return 0
            else:                    return 1
        else:
            return -1

    def __and__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value & y1.value, proto=x1)
        return x1 & y1

    def __iand__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value &= z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 & x1

    def __rand__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(y1.value & x1.value, proto=x1)
        return y1 & x1

    def __or__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value | y1.value, proto=x1)
        return x1 | y1

    def __ior__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value |= z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 | x1

    def __ror__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(y1.value | x1.value, proto=x1)
        return y1 | x1

    def __xor__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value ^ y1.value, proto=x1)
        return x1 ^ y1

    def __ixor__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value ^= z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 ^ x1

    def __rxor__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(y1.value ^ x1.value, proto=x1)
        return y1 ^ x1

    def __lshift__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value << y1.value, proto=x1)
        return x1 << y1

    def __ilshift__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value <<= z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 << x1

    def __rlshift__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return y1.__lshift__(x1)
        return y1 << x1

    def __rshift__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value >> y1.value, proto=x1)
        return x1 >> y1

    def __irshift__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            self.value >>= z.value
            self.is_signed = x1.is_signed
            self.bits = x1.num_bits
            return self
        return y1 >> x1

    def __rrshift__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return y1.__rshift__(x1)
        return y1 >> x1

    def __invert__(self):
        return Zn(~self.value)

    def __truediv__(self, y):
        if isinstance(y, Zn): y = mpf(y.value)
        return self.value / y

    def __pow__(self, y):
        y1, x1 = self._auto_cast(y)
        if isinstance(y1, Zn) and isinstance(x1, Zn):
            return Zn(x1.value ** y1.value, proto=x1)
        return x1 ** y1

class ipaddr(Zn):
    def __init__(self, value=0, cidr=None, ipvn=None, prototype=None):
        Zn.__init__(self)
        self.is_signed = False
        self.cidr = None
        if prototype is not None:
            self.ipvn = prototype.ipvn
            self.cidr = prototype.cidr
        else:
            if isinstance(value, ipaddr):
                self.ipvn = value.ipvn
                self.cidr = value.cidr
                value = value.value
            else:
                value = int(value)
            if ipvn is None:
                if value < 0xffffffff:
                    self.ipvn = 'ipv4'
                else:
                    self.ipvn = 'ipv6'
            else:
                self.ipvn = ipvn
            if isinstance(cidr, str):
                cidr = int(cidr)
            if cidr is not None:
                if cidr < 0:
                    raise ValueError("CIDR must not be negative")
                if self.ipvn == 'ipv4' and cidr > 32:
                    raise ValueError("CIDR for a IPv4 address cannot be greater than 32")
                elif self.ipvn == 'ipv6' and cidr > 128:
                    raise ValueError("CIDR for a IPv6 address cannot be greater than 128")
            self.cidr = cidr
        if self.ipvn == 'ipv4':
            self.bits = 32
        if self.ipvn == 'ipv6':
            self.bits = 128
        self.value = value
        self._update()

    def __str__(self):
        'IP address representation'
        cidr = ''
        if self.cidr is not None:
            cidr = '/%d'%self.cidr
        if self.ipvn == 'ipv6':
            pack = lambda n: n>0 and pack(n>>8)+chr(n&0xff) or ''
            v = pack(self.value)
            # ipv6 address?
            if len(v) < 16:
                v = '\x00'*(16-len(v))+v
            return ' %s%s' % (socket.inet_ntop(socket.AF_INET6, v), cidr)
        else:
            v = self.value
            v = [ (v >> 24) & 0xff, (v >> 16) & 0xff, (v >> 8) & 0xff, v & 0xff ]
            v = [ '%u'%c for c in v ]
            return ' %s%s' % ('.'.join(v), cidr)

    def __repr__(self):
        return "ipaddr(%x/%d)" % (self.n, self.cidr)

    def _check_type(self, v):
        if not isint(v):
            raise TypeError("Operations on IP addresses are limited to integers")
        if self.ipvn == 'ipv6' or v.ipvn == 'ipv6':
            ipvn = 'ipv6'
        else:
            ipvn = 'ipv4'
        return v, ipaddr(0, max(self.cidr, v.cidr), ipvn)

    def __add__(self, y):
        y = self._check_type(y)
        return ipaddr(self.value + y.value, self.cidr)

    def __iadd__(self, y):
        y = self._check_type(y)
        self.value += y.value
        return self

    def __sub__(self, y):
        y = self._check_type(y)
        return ipaddr(self.value - y.value, self.cidr)

    def __isub__(self, y):
        y = self._check_type(y)
        self.value -= y.value
        return self

    def __mul__(self, y):
        y = self._check_type(y)
        return ipaddr(self.value*y.value, self.cidr)


if __name__ == "__main__":
    # Run unit tests
    # Signed:    0  1  2  3  4  5  6  7  -8  -7  -6  -5  -4  -3  -2  -1
    # Unsigned:  0  1  2  3  4  5  6  7   8   9  10  11  12  13  14  15
    def sgn(x):
        if x < 0: return -1
        return 1
    def TestSignedPythonArithmetic(n, step=1):
        Zn().bits = n
        Zn().signed = True
        Zn().C_division = False
        m = 2**n >> 1
        it = xrange(-m, m, step)
        for i in it:
            for j in it:
                x, y = Zn(i), Zn(j)
                assert x+y == Zn((i+j))
                assert x-y == Zn((i-j))
                assert x*y == Zn((i*j))
                try: assert x//y == Zn(i//j)
                except ZeroDivisionError: pass
    def TestUnsignedPythonArithmetic(n, step=1):
        Zn().bits = n
        Zn().signed = False
        Zn().C_division = False
        b = 2**n
        it = xrange(0, b, step)
        for i in it:
            for j in it:
                x, y = Zn(i), Zn(j)
                assert x+y == Zn((i+j) % b)
                assert x-y == Zn((i-j) % b)
                assert x*y == Zn((i*j) % b)
                try: assert x/y == Zn((i//j) % b)
                except ZeroDivisionError: pass
    def TestSignedCArithmetic(n, step=1):
        Zn().bits = n
        Zn().signed = True
        Zn().C_division = True
        m = 2**n >> 1
        it = xrange(-m, m, step)
        for i in it:
            for j in it:
                x, y = Zn(i), Zn(j)
                assert x+y == Zn(i+j)
                assert x-y == Zn(i-j)
                assert x*y == Zn(i*j)
                try:
                    sign = sgn(i)*sgn(j)
                    assert x//y == Zn(sign*((abs(i) % m)//(abs(j) % m)))
                except ZeroDivisionError: pass
    def TestUnsignedCArithmetic(n, step=1):
        Zn().bits = n
        Zn().signed = False
        Zn().C_division = True
        it = xrange(0, 2**n, step)
        for i in it:
            for j in it:
                x, y = Zn(i), Zn(j)
                assert x+y == Zn(i+j)
                assert x-y == Zn(i-j)
                assert x*y == Zn(i*j)
                try:
                    assert x//y == Zn(i//j)
                except ZeroDivisionError: pass
    def BitTwiddling():
        Zn().bits = 4
        def twiddle(x, y, signed):
            u, v = Zn(x), Zn(y)
            u.signed = signed
            assert u & v == Zn(x & y)
            assert u | v == Zn(x | y)
            assert u ^ v == Zn(x ^ y)
            assert u << v == (Zn((x << y) % (Zn.high_unsigned - 1)))
            assert u >> v == Zn(x >> y)
            if y: assert u % v == Zn(x % y)
            assert ~u == Zn(~x)
            w = Zn(u.value); w &= v;  assert w == Zn(x & y)
            w = Zn(u.value); w |= v;  assert w == Zn(x | y)
            w = Zn(u.value); w ^= v;  assert w == Zn(x ^ y)
            w = Zn(u.value); w <<= v; assert w == Zn((x << y) % (Zn.high_unsigned - 1))
            w = Zn(u.value); w >>= v; assert w == Zn(x >> y)
        n = Zn.high_unsigned
        for i in xrange(n):
            for j in xrange(n):
                twiddle(i, j, True)
                twiddle(i, j, False)
        n = Zn.high_signed
        for i in xrange(-n, n):
            for j in xrange(-n, n):
                twiddle(i, j, True)
                twiddle(i, j, False)
    def Disallowed():
        '''Note that using lambda functions for the arithmetic operations
        doesn't work correctly -- e.g., 3 + Zn is allowed.
        '''
        def plus(x, y):
            return x + y
        x = Zn(1)
        oplist = (
            (lambda x, y: x + y,  "+"),
            (lambda x, y: x - y,  "-"),
            (lambda x, y: x * y,  "*"),
            (lambda x, y: x / y,  "/"),
            (lambda x, y: x // y, "//"),
            (lambda x, y: x % y,  "%"),
            (lambda x, y: x & y,  "&"),
            (lambda x, y: x | y,  "|"),
            (lambda x, y: x ^ y,  "^"),
            (lambda x, y: x << y, "<<"),
            (lambda x, y: x >> y, ">>"),
        )
        t = 3
        for op, opname in oplist:
            for arg in (t, long(t), float(t), complex(t, t)):
                try: y = op(arg, x); assert False, opname
                except TypeError: pass
                if arg != t and arg != long(t):
                    try: y = op(x, arg); assert False, opname
                    except TypeError: pass
    def TestChangingNumberOfBits():
        # Test signed
        results = (
            (15, -1, -1, -1),
            (14, -2, -2, -2),
            (13, -3, -3,  1),
            (12, -4, -4,  0),
            (11, -5,  3, -1),
            (10, -6,  2, -2),
            ( 9, -7,  1,  1),
            ( 8, -8,  0,  0),
            ( 7,  7, -1, -1),
            ( 6,  6, -2, -2),
            ( 5,  5, -3,  1),
            ( 4,  4, -4,  0),
            ( 3,  3,  3, -1),
            ( 2,  2,  2, -2),
            ( 1,  1,  1,  1),
            ( 0,  0,  0,  0),
        )
        x = Zn()
        x.bits = 4
        x.signed = True
        for item in results:
            i, compl, bits3, bits2 = item
            x.bits = 4
            x.value = i
            assert x.n == compl
            x.bits = 3
            assert x.n == bits3
            x.bits = 2
            assert x.n == bits2
        # Test unsigned
        results = (
            (15,  7,  3),
            (14,  6,  2),
            (13,  5,  1),
            (12,  4,  0),
            (11,  3,  3),
            (10,  2,  2),
            ( 9,  1,  1),
            ( 8,  0,  0),
            ( 7,  7,  3),
            ( 6,  6,  2),
            ( 5,  5,  1),
            ( 4,  4,  0),
            ( 3,  3,  3),
            ( 2,  2,  2),
            ( 1,  1,  1),
            ( 0,  0,  0),
        )
        x = Zn()
        x.bits = 4
        x.signed = False
        for item in results:
            i, bits3, bits2 = item
            x.bits = 4
            x.value = i
            assert x.n == i
            x.bits = 3
            assert x.n == bits3
            x.bits = 2
            assert x.n == bits2
    for n in xrange(1, 6):
        TestSignedPythonArithmetic(n)
        TestUnsignedPythonArithmetic(n)
        TestSignedCArithmetic(n)
        TestUnsignedCArithmetic(n)
    BitTwiddling()
    Disallowed()
    TestChangingNumberOfBits()
    exit(0)
'''
$Id: julian.py 1.15 2009/02/11 02:39:22 donp Exp $

Provides a Julian object which behaves as either an mpf or mpi number.
It represents astronomical Julian days and can be arithmetically combined
with integers, reals, etc.

The representation is based on the algorithms given in "Astronomical
Algorithms", 2nd ed., by Jean Meeus, Willman-Bell, Inc., 1998.
'''

class Julian(object):
    '''Convert a date/time specification to a Julian day object
    represented by an mpmath mpf number.  A date is interpreted as the
    Julian calendar if it is on or before 4Oct1582 and the Gregorian
    calendar if the date is on or after 15Oct1582.  A date between
    those two dates is an error if be_strict is True.
    '''

    # The following class variable is used to offset the time from the
    # conventional astronomical Julian day, which is a pure integer
    # for noon of the day it represents.  This only affects the
    # displayed date/time, as the internal representation used is
    # always the astronomical Julian day.  Example:  if you want the
    # displayed date/time to represent the conventional day beginning
    # at midnight, set the offset to -12 hours, which needs to be
    # -0.5, as the units are (Julian) days.
    day_offset = mpf(0.0)

    # The following lowercase strings are used to identify the month and
    # translated it to an integer between 1 and 12 inclusive.  Hopefully,
    # I've isolated the routines well enough that you, the user, can modify
    # the code as needed to make things work with nearly any calendar.
    month_names = {
        "Jan" : 1, "Feb" : 2, "Mar" : 3, "Apr" : 4, "May" : 5, "Jun" : 6,
        "Jul" : 7, "Aug" : 8, "Sep" : 9, "Oct" : 10, "Nov" : 11, "Dec" : 12,
    }
    name_months = dict([[val, key] for key, val in month_names.items()])

    # This is a bit of a hack, but it is intended to allow the
    # calculator to communicate how Julian dates expressed as interval
    # numbers are supposed to be represented as strings.
    interval_representation = "a"

    # If be_strict is True, there are dates that are illegal and will
    # result in an exception.
    be_strict = False

    # ---------------------------------------------------------------
    # Class variables below here are private
    fp = mpFormat()

    def __init__(self, s="now"):
        '''Initialization can be done with numerous different objects.
        The basic need is for an mpf or mpi object that represents the
        standard astronomical Julian day.  mpi objects are allowed so
        as to model uncertainty in a date.

        string:
            "h:[m[:s]]"
                Hour, minute, second format that represents a time
                today.
            "dd.dd[m[y]]"
                dd.ddd is a real number representing the day number
                (integer part) and the time as a decimal fraction.
                m is a string representing the month; the number of
                characters are in the class variable month_name_size.
                y, if present, represents the year; if not present, the
                current year is taken; y must be an integer.
            "d[m[y]][:hh[:mm[:ss]]]"
                d is an integer between 1 and 31 representing the day.
                m is a string representing the month; the number of
                characters are in the class variable month_name_size.
                y, if present, represents the year; if not present, the
                current year is taken; y must be an integer.  hh and
                mm are the integer hour and minute and must be
                integers between 0 and 59. ss is a real number of
                seconds >= 0 and < 60.
            "today"
                Represents an integer that represents the current
                integral Julian date.
            "now"
                Represents the current Julian date that represents the
                current date and time.
            NOTE:  both "today" and "now" use the python time
            library's strftime and localtime functions to get the
            date/time information.
        integer:
            Is converted to an mpf.
        rational:
            Is converted to an mpf.
        mpf:
            Represents an astronomical Julian date and, by virtue of
            the decimal fractional part, a specific time on that date.
            A fractional part of 0 represents noon on the particular
            Julian date; the class variable Julian.day_offset is first
            subtracted from the number s to get the astronomical date
            and time.
        mpi:
            Same as an mpf, but allows the representation of a time
            interval.
        '''
        if isinstance(s, int) or \
           isinstance(s, long) or \
           isinstance(s, Zn):
            self.value = mpf(int(s))
        elif isinstance(s, Rational):
            self.value = mpf(s.n)/mpf(s.d)
        elif isinstance(s, mpf):
            self.value = s
        elif isinstance(s, ctx_iv.ivmpf):
            self.value = s
        elif isinstance(s, (str, unicode)):
            y, M, d, h, m, s = self._convert_string(s.strip())
            self.value = self._to_julian(y, M, d, h, m, s)
        else:
            raise ValueError("%sUnrecognized type for Julian day" % fln())

    def _convert_string(self, s):
        '''The forms we allow are 'today', 'now', or the following two:
        1.  d.d[m[y]] where m is a month name and y is an integer year.
        2.  d[m[y]][:h[:m[:s]]] where d is an integer, m is a month name,
            and y is an integer year.  h and m are integers for the hour
            and minutes and s can be a floating point string.
        3.  :h[:m[:s]] Note the first colon is mandatory.
        '''
        def form1(s):
            day, M, y = "", "", ""
            for month in Julian.month_names:
                if month.lower() in s:
                    M = Julian.month_names[month]
                    fields = s.split(month.lower())
                    if len(fields) == 1:
                        day = fields[0]
                        break
                    elif len(fields) == 2:
                        day, y = fields
                        break
                    else:
                        raise Exception("")
            if M == "" and y == "":
                y, M = [int(i) for i in time.strftime("%Y %m").split()]
            elif y == "":
                y = int(time.strftime("%Y"))
            if y:
                y = int(y)
            if day == "":
                day = s
            day = mpf(day)
            d = int(day)
            fractional_part = day - d
            h, m, s = self._to_hms(fractional_part)
            return self._check(y, M, d, h, m, s)
        def form2(date, time):
            y, M, d, h, m, s = form1(date)
            # Now parse the time
            h, m, s = 0, 0, 0
            try:
                fields = time.split(":")
                if len(fields) == 1:
                    # Hour only
                    h = int(fields[0])
                elif len(fields) == 2:
                    # Hour and minute
                    h = int(fields[0])
                    m = int(fields[1])
                elif len(fields) == 3:
                    # Hour, minute, and seconds
                    h = int(fields[0])
                    m = int(fields[1])
                    s = mpf(fields[2])
                else:
                    raise SyntaxError("")
            except:
                raise ValueError("'%s' is a bad h:m:s specification" % time)
            msg = ""
            if not (0 <= h < 24): msg = "Bad hour specification"
            if not (0 <= m < 60): msg = "Bad minute specification"
            if not (0 <= s < 60): msg = "Bad second specification"
            if msg:
                raise ValueError("%s" % fln() + msg)
            return self._check(y, M, d, h, m, s)
        def form3(st):
            assert st[0] == ":" and len(st) > 1
            d, M, y, h, m, s = [int(i) for i in
                time.strftime("%d %m %Y %H %M %S", time.localtime()).split()]
            h, m, s = 0, 0, 0
            fields = [i for i in st[1:].split(":") if i]
            if not (1 <= len(fields) <= 3):
                raise ValueError("%s'%s' is a bad time specification" % (fln(), st))
            if len(fields) == 1:
                h = int(fields[0])
            elif len(fields) == 2:
                h, m = [int(i) for i in fields]
            elif len(fields) == 3:
                h, m = [int(i) for i in fields[:2]]
                s = mpf(fields[2])
            return self._check(y, M, d, h, m, s)
        s = s.lower()
        if s == "today" or s == "now":
            return self._convert_now(s)
        loc = s.find(":")
        try:
            if loc != -1:
                if loc == 0:
                    y, M, d, h, m, s = form3(s)
                else:
                    date = s[:loc]
                    time = s[loc+1:]
                    y, M, d, h, m, s = form2(date, time)
            else:
                y, M, d, h, m, s = form1(s)
            return y, M, d, h, m, s
        except ValueError:
            raise
        except:
            raise ValueError("%sCould not convert '%s'" % (fln(), s))

    def _convert_now(self, now):
        nowtime = time.time()
        d, M, y, h, m, s = [int(i) for i in
            time.strftime("%d %m %Y %H %M %S", time.localtime(nowtime)).split()]
        if now == "today":
            return y, M, d, 0, 0, mpf("0")
        elif now == "now":
            fp = mpf(nowtime) - long(nowtime)
            s = mpf(s) + fp
            return y, M, d, h, m, s
        else:
            raise Exception("%sProgram bug:  bad string" % fln())

    def _to_hms(self, fractional_part_of_day):
        'Return hours, minutes, seconds'
        assert 0 <= fractional_part_of_day < 1
        if 0: # Old method
            second = 24*3600*fractional_part_of_day
            fractional_seconds = second - int(second)
            second = int(second)
            minute, second = divmod(second, 60)
            hour, minute = divmod(minute, 60)
            second += fractional_seconds
        else: # New method.  See comments in _to_date()
            hours = mpf("24")*fractional_part_of_day
            # Round this value to the nearest microhour
            onemillion = mpf("1e6")
            hours += mpf("5e-6")
            hours *= onemillion
            hours = int(hours)/onemillion
            h = int(hours)
            minutes = mpf("60")*(hours - h)
            m = int(minutes)
            s = mpf("60")*(minutes - m)
            # We'll truncate to the nearest tenth of a second
            s = int(10*s)/mpf("10")
        msg = ""
        if not (0 <= s < 60):
            msg = "Program bug in calculation of seconds"
        if not (0 <= m < 60):
            msg = "Program bug in calculation of minutes"
        if not (0 <= h < 24):
            msg = "Program bug in calculation of hours"
        if msg:
            raise ValueError("%s" % fln() + msg)
        return h, m, s

    def _check(self, year, month, day, hour, minute, second):
        '''Check the arguments and return them all as integers, except
        for seconds.
        '''
        try:
            msg = ""
            if int(year) != year: msg = "Year needs to be an integer"
            if int(month) != month: msg = "Month needs to be an integer"
            if int(month) < 1 or int(month) > 12: msg = "Bad month"

            # If day is an mpf, then we must have the time being 0
            if isinstance(day, mpf):
                if hour != 0 or minute != 0 or second != 0:
                    msg = "hr, min, sec must be 0 if day is real number"
                fractional_part = day - int(day)
                day = int(day)
                hour, minute, second = self._to_hms(fractional_part)
            else:
                if int(day) != day: msg = "Day needs to be an integer"
            if int(hour) != hour: msg = "Hour needs to be an integer"
            if int(minute) != minute: msg = "Minute needs to be an integer"
            if not isinstance(second, int) and \
               not isinstance(second, long) and \
               not isinstance(second, Zn) and \
               not isinstance(second, mpf):
                msg = "Second not of proper type"
            M, d, h, m = [int(i) for i in (month, day, hour, minute)]
            biggest_day = self._days_in_month(M, year)
            if d < 1 or d > biggest_day: msg = "Bad day"
            if h < 0 or h > 23: msg = "Bad hour"
            if m < 0 or m > 59: msg = "Bad minute"
            if second < 0 or second >= 60: msg = "Bad second"
            if year == 1582 and month == 10:
                if d > 4 and d < 15:
                    if Julian.be_strict:
                        msg = "Illegal date for Gregorian and Julian calendars"
                    else:
                        # Set to one or the other
                        if d <= 9: d = 4
                        else: d = 15
            if msg:
                raise ValueError("%s" % fln() + msg)
            return year, M, d, h, m, mpf(second)
        except ValueError:
            raise
        except:
            raise Exception("%sBad date/time" % fln())

    def _to_julian(self, year, month, day, hour=0, minute=0, second=0):
        'Algorithm is on page 61 of Meeus'
        Y, M, D, h, m, s = self._check(year, month, day, hour, minute, second)
        # Add the time to D; D will become a real number
        D += (mpf(h) + mpf(m)/mpf(60) + mpf(s)/mpf(3600))/mpf(24)
        if M == 1 or M == 2:
            Y -= 1
            M += 12
        A = Y//100
        B = 2 - A + A//4
        if Y <= 1582:
            # Set B = 0 if we're in the Julian calendar
            if Y == 1582:
                if M < 10 or (M == 10 and D < 5):
                    B = 0
            else:
                B = 0
        jd = int(mpf("365.25")*(Y + 4716)) + int(mpf(30.6001)*(M + 1)) + \
             D + B - mpf("1524.5")
        return jd

    def _is_leap_year(self, y):
        '''Rule is on page 62 of Meeus.  In the Julian calendar, a
        year is a leap year if it is divisible by 4.
        '''
        if not (isinstance(y, int) or isinstance(y, long) or isinstance(y, Zn)):
            raise ValueError("%sYear must be an integer" % fln())
        if y <= 1582:
            return (y % 4 == 0)
        else:
            # For Gregorian calendar
            if (y % 400 == 0) or ((y % 4 == 0) and (y % 100 != 0)):
                return True
            return False

    def _days_in_month(self, month, year):
        '''Assumes month and year are integers.
        '''
        if month == 2:
            if self._is_leap_year(year): return 29
            else: return 28
        elif month in (4, 6, 9, 11): return 30
        elif month in (1, 3, 5, 7, 8, 10, 12): return 31
        else: raise ValueError("%sBad month" % fln())

    def _to_date(self, jd):
        'Algorithm on page 63 of Meeus.'
        if not isinstance(jd, mpf):
            raise ValueError("%sJulian day must be a real number" % fln())
        if jd < 0:
            raise ValueError("%sJulian day cannot be negative" % fln())
        jd += mpf("0.5")
        Z = int(jd)
        F = jd - Z
        if Z < 2299161:
            A = Z
        else:
            alpha = int((Z - mpf("1867216.25"))/mpf("36524.25"))
            A = Z + 1 + alpha - alpha//4
        B = A + 1524
        C = int((B - mpf("122.1"))/mpf("365.25"))
        D = int(mpf("365.25")*C)
        E = int((B - D)/mpf("30.6001"))
        D = B - D - int(mpf("30.6001")*E) + F
        if E < 14:
            M = E - 1
        else:
            M = E - 13
        if M > 2:
            y = C - 4716
        else:
            y = C - 4715
        d = int(D)
        fractional_part = D - d
        # Note:  a test case of ":22" for jd = 2454871.4166666665
        # showed a problem with roundoff, in that the decimal hours
        # calculate to 21.99999 99960 to 10 places.  This value
        # results in the following tuple being returned:
        # (2009, 2, 8, 21, 59, mpf('59.999986588954926')).  Because of
        # this test case, I decided to round hours to the nearest
        # microhour by adding 5e-7, etc. in _to_hms().
        h, m, s = self._to_hms(fractional_part)
        return y, M, d, h, m, s

    def _st(self, val):
        assert isinstance(val, mpf)
        Julian.fp.digits(1)  # We'll display to the nearest tenth second
        if val < 250000:
            return self._units(val)
        y, M, d, h, m, s = self._to_date(val)
        try:
            t = ["%d %s %d" % (d, Julian.name_months[M], y)]
            if max(h,m,s) != 0:
                t += [" %02d:%02d" % (h, m)]
                if s != 0:
                    t += [":%02d.%01d" % (int(s), int(((10*s)-(10*int(s)))/10)) ]
            return ''.join(t)
        except:
            msg = "%sDate representation cannot be calculated\n" + \
                  "Try increasing the precision with `prec`."
            raise ValueError(msg % fln())

    def _units(self, val):
        '''Return a string representing val as a number of common time units.
        For example, if val == 1, this is 1 day, so '1 day' would be returned.
        '''
        y = mpf("365.25")
        w = mpf("7")
        mo = y/12
        s = mpf("86400")
        m = mpf("1440")
        h = mpf("24")
        Julian.fp.digits(3)
        f = Julian.fp.sig
        if val < 0:
            sign = -1
            val = abs(val)
        else:
            sign = 1
        if 0 <= val <= 1/m:
            return f(sign*s*val) + " seconds"
        elif 1/m < val <= 1/h:
            return f(sign*m*val) + " minutes"
        elif 1/h < val <= 1:
            return f(sign*h*val) + " hours"
        elif 1 < val <= w:
            return f(sign*val)   + " days"
        elif w < val <= mo:
            return f(sign*val/w) + " weeks"
        elif mo < val <= y:
            return f(sign*val/mo) + " months"
        else:
            return f(sign*val/y) + " years"

    def __str__(self):
        if self.value == inf: return "Julian(inf)"
        if self.value == -inf: return "Julian(-inf)"
        val = self.value - Julian.day_offset
        if isinstance(val, ctx_iv.ivmpf):
            if Julian.interval_representation == "a":
                a = self._st(val.mid)
                b = self._units(val.delta/mpf("2"))
                return a + " +-" + b
            elif Julian.interval_representation == "b":
                a = self._st(val.mid)
                p = mpf("100")*val.delta/(mpf("2")*val.mid)
                Julian.fp.digits(3)
                b = Julian.fp.sig(p).strip()
                return a + " (" + b + "%)"
            elif Julian.interval_representation == "c":
                s = " <<" + self._st(val.a) + ", " + \
                     self._st(val.b) + ">>"
                return s
            else:
                raise Exception("%sBad Julian.interval_representation" % fln())
        else:
            return " " + self._st(val)

    def __repr__(self):
        return "Julian(" + repr(self.value) + ")"

    def _convert_to_mpf_or_mpi(self, n):
        if isinstance(n, int) or isinstance(n, long):
            return mpf(n)
        if isinstance(n, Zn):
            return mpf(int(n))
        if isinstance(n, Rational):
            return n.mpf()
        if isinstance(n, mpc):
            return abs(n)
        if isinstance(n, mpf) or isinstance(n, ctx_iv.ivmpf):
            return n
        if isinstance(n, Julian):
            return n.value
        raise ValueError("%sBad type for operation with date/time" % fln())

    def __add__(self, other):
        return Julian(self.value + self._convert_to_mpf_or_mpi(other))
    def __sub__(self, other):
        return Julian(self.value - self._convert_to_mpf_or_mpi(other))
    def __mul__(self, other):
        return Julian(self.value * self._convert_to_mpf_or_mpi(other))
    def __div__(self, other):
        return Julian(self.value / self._convert_to_mpf_or_mpi(other))

    def __radd__(self, other):
        return Julian(self._convert_to_mpf_or_mpi(other) + self.value)
    def __rsub__(self, other):
        return Julian(self._convert_to_mpf_or_mpi(other) - self.value)
    def __rmul__(self, other):
        return Julian(self._convert_to_mpf_or_mpi(other) * self.value)
    def __rdiv__(self, other):
        raise Exception("%sMeaningless to divide by date/time" % fln())
    def __neg__(self):
        self.value = -self.value
        return self

    def __int__(self):
        if isinstance(self.value, mpf):
            return int(self.value)
        elif isinstance(self.value, ctx_iv.ivmpf):
            return int(self.value.mid)
        else:
            raise Exception("%sProgram bug:  unknown type" % fln())

    __long__ = __int__

    def to_mpf(self):
        if isinstance(self.value, mpf):
            return self.value
        elif isinstance(self.value, ctx_iv.ivmpf):
            return self.value.mid
        else:
            raise Exception("%sProgram bug:  unknown type" % fln())

if __name__ == "__main__":
    def TestNumericalInit():
        j = Julian(0)
        data = (
            # Test vectors from Meeus pg 62
            ( 1957, 10,           4, 19, 26, 24, mpf("2436116.31")),
            (  333,  1,          27, 12,  0,  0, mpf("1842713.0")),
            ( 2000,  1,  mpf("1.5"),  0,  0,  0, mpf("2451545.0")),
            ( 1999,  1,  mpf("1.0"),  0,  0,  0, mpf("2451179.5")),
            ( 1988,  6, mpf("19.5"),  0,  0,  0, mpf("2447332.0")),
            ( 1600,  1,  mpf("1.0"),  0,  0,  0, mpf("2305447.5")),
            ( 1600, 12, mpf("31.0"),  0,  0,  0, mpf("2305812.5")),
            (-123,  12, mpf("31.0"),  0,  0,  0, mpf("1676496.5")),
            (-122,   1,  mpf("1.0"),  0,  0,  0, mpf("1676497.5")),
            (-1000,  2, mpf("29.0"),  0,  0,  0, mpf("1355866.5")),
            (-4712,  1,  mpf("1.5"),  0,  0,  0, mpf("0.0")),
        )
        for item in data:
            y, M, d, h, m, s, jd = item
            assert j._to_julian(y, M, d, h, m, s) == jd
            y1, M1, d1, h1, m1, s1 = j._to_date(jd)
            if isinstance(d, mpf):
                D = d1 + (mpf(h1) + mpf(m1)/60 + mpf(s1)/3600)/mpf(24)
                assert (y, M, d) == (y1, M1, D)
            else:
                # Note this comparison is just to the nearest second
                assert (y, M, d, h, m, s) == (y1, M1, d1, h1, m1, int(s1))
    def TestStringInit():
        data = (
            # Test points from Meeus pg 62
            ("4Oct1957:19:26:24", mpf("2436116.31")),
            ("27Jan333:12",       mpf("1842713.0")),
            ("1.5jan2000",        mpf("2451545.0")),
            ("1.0jan1999",        mpf("2451179.5")),
            ("19.5JUN1988",       mpf("2447332.0")),
            ("1.0jan1600",        mpf("2305447.5")),
            ("31.0Dec1600",       mpf("2305812.5")),
            ("31.0Dec-123",       mpf("1676496.5")),
            ("1.0jan-122",        mpf("1676497.5")),
            ("29Feb-1000",        mpf("1355866.5")),
            ("1.5Jan-4712",       mpf("0.0")),
        )
        for s, jd in data:
            assert Julian(s).value == jd
        # Test the :h:m:s forms
        Julian.day_offset = mpf("0")
        fields = str(Julian(":2")).split(":")
        assert fields[1] == "02"; assert fields[2] == "00"
        fields = str(Julian(":22")).split(":")
        assert fields[1] == "22"; assert fields[2] == "00"
        fields = str(Julian(":22:45")).split(":")
        assert fields[1] == "22"; assert fields[2] == "45"
        fields = str(Julian(":22:45:12.3")).split(":")
        assert fields[1] == "22"; assert fields[2] == "45"
        assert fields[3] == "12.3"
    def TestStringRepresentations():
        # Test string representations (note leading spaces)
        Julian.day_offset = mpf("0")
        Julian.interval_representation = "c"
        j1, j2 = mpf("2451545.0"), mpf("2451545.5")
        assert str(Julian(j1)) == " 1Jan2000:12:00"
        assert str(Julian(mpi(j1, j2))) == " <<1Jan2000:12:00, 2Jan2000:00:00>>"
        Julian.day_offset = mpf("0.5")
        assert str(Julian(j1)) == " 1Jan2000:00:00"
        assert str(Julian(mpi(j1, j2))) == " <<1Jan2000:00:00, 1Jan2000:12:00>>"
    def TestArithmetic():
        Julian.day_offset = mpf("0")
        j = mpf("2451545.0")
        assert str(Julian(j)) == " 1Jan2000:12:00"
        assert str(Julian(j) + mpf("1")) == " 2Jan2000:12:00"
        assert str(Julian(j) - mpf("1")) == " 31Dec1999:12:00"
        assert str(Julian(j)*2) == " 8Feb8712:12:00"
        assert str(Julian(j)/2) == " 26Dec-1357:00:00"
        assert str(mpf("1") + Julian(j)) == " 2Jan2000:12:00"
        assert str(mpf("0") - Julian(-j)) == " 1Jan2000:12:00"
        assert str(2*Julian(j)) == " 8Feb8712:12:00"
    TestNumericalInit()
    TestStringInit()
    TestStringRepresentations()
    TestArithmetic()

# Number recognition regular expressions
integer = re.compile("^[+-]?\d+$")

cre=r'''
    %s                          # Match at beginning
    ([+-])%s                    # Optional leading sign
    %s                          # Placeholder for imaginary unit
    (\d+\.\d+|\d+\.?|\.\d+)     # Required digits and opt. decimal point
    (e[+-]?\d+)?                # Optional exponent
    %s                          # Match at end
'''
# Pure imaginary, xi or ix
I1 = cre % ("^", "?", "", "i$")
I2 = cre % ("^", "?", "i", "$")

# Reals
R = cre % ("^", "?", "", "$")

# True complex numbers:  x+iy or x+yi
C1 = (cre % ("^", "?", "", "")) + (cre % ("", "", "", "[ij]$"))
C2 = (cre % ("^", "?", "", "")) + (cre % ("", "", "[ij]", "$"))

# Degenerate complex numbers:  1-i, 3.7+i
num = r"([+-]?)(\d+\.\d+|\d+\.?|\.\d+)(e[+-]?\d+)?"
C3 = r"^%s([+-][ij])$" % num

# True complex numbers:  (x,y)
C4 = r"^\(%s(,)%s\)$" % (num, num)
del num

# Regular expressions
imag1 =    re.compile(I1, re.X | re.I)
imag2 =    re.compile(I2, re.X | re.I)
real =     re.compile(R,  re.X | re.I)
complex1 = re.compile(C1, re.X | re.I)
complex2 = re.compile(C2, re.X | re.I)
complex3 = re.compile(C3, re.X | re.I)
complex4 = re.compile(C4, re.X | re.I)

# Rationals:  "a/b", and "axb/c" forms are allowed where a and b are
# integers and x is one or more of the following characters: '+- '.
# b and c are always positive; a may be positive or negative.
Ra = r'''
    ^   # Mixed fraction
    ([-+])?             # Optional sign
    (\d+)               # Integer
    [-+ ]+              # Separation character
    (\d+)               # Numerator
    /                   # Fraction separator
    (\d+)               # Denominator
    $
    |                   # Or
    ^   # Canonical fraction
    ([-+]?\d+)          # Integer a with optional sign
    /                   # Fraction separator
    (\d+)               # Denominator
    $
'''
rational = re.compile(Ra, re.X | re.I)

ip = re.compile(r"^(\d{1,3})[.](\d{1,3})[.](\d{1,3})[.](\d{1,3})")
ip6 = re.compile(r"""^(([0-9a-f]{1,4}:){1,1}(:[0-9a-f]{1,4}){1,6})$|
^(([0-9a-f]{1,4}:){1,2}(:[0-9a-f]{1,4}){1,5})$|
^(([0-9a-f]{1,4}:){1,3}(:[0-9a-f]{1,4}){1,4})$|
^(([0-9a-f]{1,4}:){1,4}(:[0-9a-f]{1,4}){1,3})$|
^(([0-9a-f]{1,4}:){1,5}(:[0-9a-f]{1,4}){1,2})$|
^(([0-9a-f]{1,4}:){1,6}(:[0-9a-f]{1,4}){1,1})$""")
ip6 = re.compile(r"""
    ^((([0-9a-f]{1,4}:){1,6})(:[0-9a-f]{1,4}){1,1})$|
    ^((([0-9a-f]{1,4}:){1,5})(:[0-9a-f]{1,4}){1,2})$|
    ^((([0-9a-f]{1,4}:){1,4})(:[0-9a-f]{1,4}){1,3})$|
    ^((([0-9a-f]{1,4}:){1,3})(:[0-9a-f]{1,4}){1,4})$|
    ^((([0-9a-f]{1,4}:){1,2})(:[0-9a-f]{1,4}){1,5})$|
    ^((([0-9a-f]{1,4}:){1,1})(:[0-9a-f]{1,4}){1,6})$|
    ^((([0-9a-f]{1,4}:){7,7})([0-9a-f]{1,4}))$|
    ^((([0-9a-f]{1,4}:){1,7}|:):)$|
    ^(:(:[0-9a-f]{1,4}){1,7})$|
    ^(((([0-9a-f]{1,4}:){6})(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3}))$|
    ^((([0-9a-f]{1,4}:){5}[0-9a-f]{1,4}:(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3}))$|
    ^(([0-9a-f]{1,4}:){5}:[0-9a-f]{1,4}:(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3})$|
    ^(([0-9a-f]{1,4}:){1,1}(:[0-9a-f]{1,4}){1,4}:(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3})$|
    ^(([0-9a-f]{1,4}:){1,2}(:[0-9a-f]{1,4}){1,3}:(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3})$|
    ^(([0-9a-f]{1,4}:){1,3}(:[0-9a-f]{1,4}){1,2}:(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3})$|
    ^(([0-9a-f]{1,4}:){1,4}(:[0-9a-f]{1,4}){1,1}:(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3})$|
    ^((([0-9a-f]{1,4}:){1,5}|:):(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3})$|
    ^(:(:[0-9a-f]{1,4}){1,5}:(25[0-5]|2[0-4]\d|[0-1]?\d?\d)(\.(25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3})$
""", re.X | re.I)

vector = re.compile(r"\[\s*([.0])")

class Number(object):
    '''Used to generate a number object from a string.
    '''

    # If the following class variable is nonzero, then integers are
    # being restricted to a said number of bits.
    bits = 0

    # If signed is True, then integer arithmetic is signed; otherwise,
    # arithmetic is unsigned.
    signed = True

    def __init__(self):
        pass

    def __call__(self, s, tags=None):
        assert len(s) > 0
        suffix = 1
        if tags is not None:
            if 'ipaddr' in tags:
                return self.ip(s, tags)
        if s != "now" and s != "today":
            if len(s) > 1 and s[:2] != "0x":
                if s[-1] in suffixes_ln:
                    exponent = suffixes_ln[s[-1]]
                    if exponent >= 0:
                        suffix = Zn(10**exponent)
                    else:
                        suffix = mpf("1e" + str(exponent))
                    s = s[:-1]
        for func in (self.ip, self.j, self.i, self.q, self.v, self.r, self.c):
            x = func(s)
            if x != None:
                if suffix == 1:
                    return x
                if isint(x):
                    if isint(suffix): return Zn(suffix*x)
                    else:             return suffix*mpf(int(x))
                elif isinstance(x, Rational):
                    if isint(suffix): return Rational(suffix*x.n, x.d)
                    else:             return suffix*x.mpf()
                elif isinstance(x, mpf) or \
                     isinstance(x, mpc) or \
                     isinstance(x, mpi):
                    if isint(suffix): return mpf(int(suffix))*x
                    else:             return suffix*x
                else:
                    return None
        return None

    def j(self, s):
        '''Check to see if it's a Julian date/time form.  We only allow
        two forms:  'dS[y[:...]]' where S is a string for the month or
        '.+:.+' (regexp syntax) where it contains a colon and means a
        time today.  'now' and 'today' are also allowed.
        '''
        if s == "now" or s == "today" or ":" in s:
            return Julian(s)
        sl = s.lower()
        for month in Julian.month_names:
            if month.lower() in sl:
                return Julian(s)
        return None

    def i(self, s):
        # Handle special cases like 0x, 0o, and 0b
        try:
            value = 0
            match = False
            if len(s) > 2:
                if s[:2] == "0x":
                    value = int(s[2:], 16)
                    match = True
                elif s[:2] == "0o":
                    value = int(s[2:], 8)
                    match = True
                elif s[:2] == "0b":
                    value = int(s[2:], 2)
                    match = True
            if integer.match(s):
                return Zn(int(s))
            if match:
                return Zn(value)
        except ValueError:
            pass
        except Exception:
            raise
        return None

    def A(self, s):
        # is this is an array

        return None

    def V(self, s):
        # is this is a vector
        return None

    def ip(self, s, tags=None):
        def unpack(s):
            v = '0x'
            for i in range(len(s)):
                v += '%02x' % ord(s[i])
            return int(v,16)
        cidr = None
        if '/' in s:
            sparts = s.split('/')
            s = sparts[0]
            cidr = sparts[1]
        try:
            if tags is None:
                tags = ['ipv4', 'ipv6']
            if 'ipv4' in tags:
                mo = ip.match(s)
                if cidr is None:
                    cidr = 32
                if mo:
                    dquad = [ int(i) for i in mo.groups() if i ]
                    if max(dquad) > 255:
                        return None
                    ps = socket.inet_pton(socket.AF_INET, s)
                    return ipaddr(unpack(ps), cidr, 'ipv4')
            if 'ipv6' in tags:
                if cidr is None:
                    cidr = 128
                if ip6.match(s):
                    ps = socket.inet_pton(socket.AF_INET6, s)
                    return ipaddr(unpack(ps), cidr, 'ipv6')
        except Exception, e:
            print e
            pass
        return None

    def q(self, s):
        mo = rational.match(s)
        if mo:
            g = [i for i in mo.groups() if i]
            sign = ""
            if g[0] == "+":
                del g[0]
            elif g[0] == "-":
                sign = "-"
                del g[0]
            try:
                num = [int(i) for i in g]
            except:
                raise Exception("Bug:  rational match on non-integer")
            if not num:
                return None
            elif len(num) == 2:
                w, n, d = 0, num[0], num[1]
            elif len(num) == 3:
                w, n, d = num
            else:
                msg = "Program bug\nUnexpected number of matches on\n"
                msg += "'%s'" % s
                raise Exception(msg)
            n = int(sign + str(w*d + n))
            return Rational(n, d)
        else:
            return None

    def r(self, s):
        # Handle infinities
        if s == "inf": return inf
        if s == "-inf": return -inf
        # If the number begins with "E" or 'e', prepend a 1
        if s[0] == "E" or s[0] == "e":
            s = "1" + s
        mo = real.match(s)
        if mo:
            num = [i for i in mo.groups() if i]
            if not num:
                return None
            else:
                try:
                    return mpf(''.join(num))
                except:
                    return None
        else:
            return None

    def c(self, s):
        s = s.lower()
        s = s.replace("j", "i")
        if "i" not in s and "(" not in s:
            return None
        if s == "i" or s == "+i":
            return mpc(0, 1)
        if s == "-i":
            return mpc(0, -1)
        mo = complex3.match(s)
        if mo:
            n = ''.join([i for i in mo.groups() if i])
            try:
                ending = n[-2:]
                if ending == "+i":
                    return mpc(mpf(n[:-2]), 1)
                elif ending == "-i":
                    return mpc(mpf(n[:-2]), -1)
                else:
                    raise Exception("Program bug:  unexpected complex number")
            except:
                pass
        for expression in (imag1, imag2):
            mo = expression.match(s)
            if mo:
                num = [i for i in mo.groups() if i]
                if num:
                    try:
                        return mpc(0, mpf(''.join(num)))
                    except:
                        return None
                else:
                    return None
        for expression in (complex1, complex2):
            mo = expression.match(s)
            if mo:
                try:
                    g = mo.groups()
                    r = mpf(''.join([i for i in g[:3] if i]))
                    i = mpf(''.join([i for i in g[3:] if i]))
                    return mpc(r, i)
                except:
                    return None
        mo = complex4.match(s)
        if mo:
            try:
                s = "".join([i for i in mo.groups() if i])
                r, i = [mpf(i) for i in s.split(",")]
                return mpc(r, i)
            except:
                return None

    def v(self, s):
        '''Interval numbers:  allowed forms are
            1. 'a +- b'
            2. 'a (b%)'  % sign is optional
            3. '[a, b]'
        In 1, a is the midpoint of the interval and b is the half-width.
        In 2, a is the midpoint of the interval and b is the half-width.
        In 3, the interval is indicated directly.
        '''
        e = ValueError("Improperly formed interval number '%s'" %s)
        s = s.replace(" ", "")
        if "+-" in s:
            n = [mpf(strip(i)) for i in s.split("+-")]
            return mpi(n[0] - n[1], n[0] + n[1])
        elif "(" in s:
            if s[0] == "(":  # Don't confuse with a complex number (x,y)
                return None
            if ")" not in s:
                raise e
            s = s.replace(")", "")
            percent = False
            if "%" in s:
                if s[-1] != "%":
                    raise e
                percent = True
                s = s.replace("%", "")
            a, p = [mpf(strip(i)) for i in s.split("(")]
            d = p
            if percent:
                d = a*p/mpf(100)
            return mpi(a - d, a + d)
        elif "," in s:
            if "[" not in s: raise e
            if "]" not in s: raise e
            s = s.replace("[", "")
            s = s.replace("]", "")
            n = [mpf(strip(i)) for i in s.split(",")]
            return mpi(n[0], n[1])
        else:
            return None

if __name__ == "__main__":
    # Test cases
    nums = {
        # Integers
        Zn(0) : (
            "0", "+0", "-0", "000", "+000", "-000",
            ),
        Zn(1) : (
            "1", "+1", "01", "+01", "001", "+001",
            ),
        Zn(-1) : (
            "-1", "-01", "-001",
            ),
        Zn(124) : (
            "124", "+124", "0124", "+0124", "000124", "+000124",
            ),
        Zn(-123) : (
            "-123", "-000123",
            ),

        # Reals
        mpf(0) : (
            "0.0", "+0.0", "-0.0",
            "000.000", "+000.000", "-000.000",
            ),
        mpf(1) : (
            "1.", "+1.", "1.0", "+1.0",
            "1.0e0", "+1.0e0",
            "1.0E0", "+1.0E0",
            ),
        mpf(-2) : (
            "-2.", "-2.0", "-2.0e0", "-2.0000E000",
            ),
        mpf("-2.3") : (
            "-2.3", "-2.30", "-2.3000", "-2.3e0", "-2300e-3", "-0.0023e3",
            "-.23E1",
            ),
        mpf("2.345e-7") : (
            "2.345e-7", "2345e-10", "0.00000002345E+1", "0.0000002345",
            ),

        # Pure imaginaries
        mpc(0, 1) : (
            "i", "+i", "1i", "+1i", "+i1",
            "1.i", "+1.i", "+i1.",
            "1.0i", "+1.0i", "+i1.0",
            "1.00i", "+1.00i", "+i1.00",
            ),
        mpc(0, -1) : (
            "-i", "-i", "-1i", "-i1",
            "-1.i", "-1.i", "-i1.",
            "-1.0i", "-1.0i", "-i1.0",
            "-1.00i", "-1.00i", "-i1.00",
            ),
        mpc(0, 3) : (
            "3i", "+3i", "3.i", "+3.i", "3.0i", "+3.0i", "3.0e0i", "+3.0e0i",
            "I3", "+I3", "I3.", "+I3.", "I3.0", "+I3.0", "I3.0e0", "+I3.0e0",
            "3.000i", "i3.000", "3.000E0i", "i3.000E0",
            "3.000e-0J", "J3.000e-0", "3.000e+0J", "J3.000e+0",
            ),
        mpc(0, -8) : (
            "-8i", "-8.i", "-8.0i", "-8.0e0i",
            "-j8", "-j8.", "-j8.0", "-j8.0E0",
            ),
        mpc(0, mpf("-0.123")) : (
            "-.123i", "-.123j", "-0.123i", "-1.23e-1i",
            ),

        # Complex numbers
        mpc(1, -1) : (
            "1-1i", "1-1.i", "1.-1i", "1.-1.i",
            "1-j1", "1-j1.", "1.-j1", "1.-j1.",
            "1.00-1.00I", "1.00-I1.00", "1000e-3-100000e-5I",
            "1.00-J1.00", "1000E-3-J100000E-5",
            ),
        mpc(1, 1) : (
            "1+1i", "1+1.i", "1.+1i", "1.+1.i",
            "1+i", "1.0+i", "1.000+i",
            ),
        mpc("4.9", -1) : (
            "4.9-1i", "4.9-1.i", "49e-1-1i",
            "4.9-i", "49e-1-i",
            ),
        mpc("-7", -1) : (
            "-7-1i", "-7-1.i", "-70e-1-1i",
            "-7-i", "-70e-1-i",
            ),
        mpc(mpf("11.549e-59"), mpf("-8.31e89")) : (
            "11.549e-59-8.31e89I", "1.1549e-58-J831e87",
            ),
        mpc(1, 2): (
            "(1,2)", "(1.,2.)", "(1.0,2.0)", "(1.000,2.000)",
        ),

        # Rational numbers
        Rational(3, 8) : (
            "3/8", "6/16", "0 12/32", "0-15/40", "0+18/48",
            ),
        Rational(-3, 7) : (
            "-3/7", "-6/14", "-0 12/28", "-0-15/35", "-0+18/42",
            ),
        Rational(3, -7) : (
            "-3/7", "-6/14", "-0 12/28", "-0-15/35", "-0+18/42",
            ),
    }
    # Because of a bug in mpi == and != tests, we have to test them
    # differently.
    mpi_tests = {
        # Interval numbers
        mpi(1, 3) : (
            "[1, 3]", "[1.0, 3]", "[1, 3.0]", "[1.0, 3.0]",
            "[1,3]", "[1.0,3]", "[1,3.0]", "[1.0,3.0]",
            "[   1,     3]", "[   1.0,3]", "[     1,     3.0]",
            "1.5 +- 0.5", "1.5+-0.5", "1.5+-      0.5", "1.5     +-0.5",
            "1.5      +-      0.5", "15e-1 +- 500e-3",
            "1.5(33.33333333333333333333%)", "1.5  (33.33333333333333333333%)",
            "1.5    (     33.33333333333333333333%)",
            "1.5    (     33.33333333333333333333     % )",
            "1.5(     33.33333333333333333333    %  )",
            "1.5(33.33333333333333333333    %  )",
            ),
    }

    n = Number()
    status = 0
    for number in nums:
        for numstr in nums[number]:
            num = n(numstr)
            if num != number:
                print("Error for '%s'" % numstr)
                print("  Should be %s" % str(number))
                print("  Got       %s" % str(num))
                status += 1
    for number in mpi_tests:
        for numstr in mpi_tests[number]:
            num = n(numstr)
            if (num.a != number.a) and (num.b != number.b):
                print("Error for '%s'" % numstr)
                print("  Should be %s" % str(number))
                print("  Got       %s" % str(num))
                status += 1
    exit(status)


# Strings to identify types
INT = "z"
RAT = "q"
MPF = "r"
MPC = "c"
MPI = "i"
JUL = "t"

def Convert(x, arg_type, digits=0):
    '''Converts amongst the numerical types.  Some conversions lose
    information.  The digits argument controls the precision of a conversion
    of a real to a rational.
    '''
    e = SyntaxError("Unknown type")
    if arg_type == INT:
        if isint(x):                  return Zn(x)
        elif isinstance(x, Rational): return Zn(int(mpf(x.n)/mpf(x.d)))
        elif isinstance(x, mpf):      return Zn(int(x))
        elif isinstance(x, mpc):      return Zn(int(abs(x)))
        elif isinstance(x, ctx_iv.ivmpf):      return Zn(int(x.mid))
        elif isinstance(x, Julian):   return Zn(int(x))
        else: raise e
    elif arg_type == RAT:
        if isint(x):                  return Rational(int(x), 1)
        elif isinstance(x, Rational): return x
        elif isinstance(x, mpf):      return Rational().frac(x, digits)
        elif isinstance(x, mpc):      return Rational().frac(abs(x), digits)
        elif isinstance(x, ctx_iv.ivmpf):      return Rational(x.mid)
        elif isinstance(x, Julian):   return Rational().frac(x.to_mpf(), digits)
        else: raise e
    elif arg_type == MPF:
        if isint(x):                  return mpf(int(x))
        elif isinstance(x, Rational): return x.mpf()
        elif isinstance(x, mpf):      return x
        elif isinstance(x, mpc):      return abs(x)
        elif isinstance(x, ctx_iv.ivmpf):      return x.mid
        elif isinstance(x, Julian):   return x.to_mpf()
        else: raise e
    elif arg_type == MPC:
        if isint(x):                  return mpc(int(x))
        elif isinstance(x, Rational): return x.mpc()
        elif isinstance(x, mpf):      return mpc(x, 0)
        elif isinstance(x, mpc):      return x
        elif isinstance(x, ctx_iv.ivmpf):      return mpc(x.mid, 0)
        elif isinstance(x, Julian):   return mpc(x.to_mpf(), 0)
        else: raise e
    elif arg_type == MPI:
        if isint(x):                  return mpi(int(x))
        elif isinstance(x, Rational): return x.mpi()
        elif isinstance(x, mpf):      return mpi(x)
        elif isinstance(x, mpc):      return mpi(abs(x))
        elif isinstance(x, ctx_iv.ivmpf):      return x
        elif isinstance(x, Julian):
            if isinstance(x.value, mpf):  return mpi(x.value)
            else:                         return x.value
        else: raise e
    elif arg_type == JUL:
        if isint(x):                  return Julian(int(x))
        elif isinstance(x, Rational): return Julian(mpf(x.n)/mpf(x.d))
        elif isinstance(x, mpf):      return Julian(x)
        elif isinstance(x, mpc):      return Julian(abs(x))
        elif isinstance(x, ctx_iv.ivmpf):      return Julian(x)
        elif isinstance(x, Julian):   return x
        else: raise e
    else:
        raise SyntaxError("Unknown type")

if __name__ == "__main__":
    # Unit tests
    def TestConvert():
        n = 1
        number_types = (
            Zn(n),
            Rational(n, n),
            mpf(n),
            mpc(n, n),
            mpi(n),
            Julian(n),
        )
        results = (
            (INT, Zn),
            (RAT, Rational),
            (MPF, mpf),
            (MPC, mpc),
            (MPI, mpi),
            (JUL, Julian)
        )
        for number in number_types:
            for typename, type in results:
                assert isinstance(Convert(number, typename), type)
    TestConvert()
    exit(0)
