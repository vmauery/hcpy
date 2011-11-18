'''
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


from mpmath import mpf, mp, nan, inf, nstr
from mpmath.libmp.libmpf import to_digits_exp, fzero, finf, fninf, fnan
from si import suffixes_nl

try: from pdb import xx  # pdb.set_trace is xx; easy to find for debugging
except: pass

class mpFormatException(Exception): pass

class mpFormat(object):
    '''
    The format(number, format_type="sig") is a convenience function
    to let you specify the type of formatting you want with a string.

    Methods:
        digits(n)           Sets number of figures after the decimal
                            point for fix and significant figures for
                            the other modes.
        fix(number)         Fixed decimal point mode
        sig(number)         Significant figure mode
        sci(number)         Scientific notation
        eng(number)         Engineering notation
        engsi(number)       Engineering notation with an SI suffix

    For an example, run the following code:

        from mpformat import *
        from mpmath import *

        mp.dps = 30
        fp = mpFormat()
        mpFormat.comma_decorate = True
        x = mpf("12345.67890123456789")
        fp.digits(10)
        print fp.fix(x)
        print fp.sig(x)
        print fp.sci(x)
        print fp.eng(x)
        print fp.engsi(x)
        mpFormat.cuddle_si = True
        print fp.engsi(x)

    to get:

         12,345.67890 12346
         12,345.67890
         1.234567890e+4
         12.34567890e+3
         12.34567890 k
         12.34567890k

    Class variables:

        fix_low     If the number is less than this, the fix mode
                    underflows to the sci format.

        fix_high    If the number is greater than this, the fix mode
                    overflows to the sci format.

        sig_low     If the number is less than this, the sig mode
                    underflows to the sci format.

        sig_high    If the number is greater than this, the sig mode
                    overflows to the sci format.

    '''

    if 1:
        fix_low  = mpf("1e-11")
        fix_high = mpf("1e12")
        sig_low  = mpf("1e-11")
        sig_high = mpf("1e12")
    else:
        # For debugging
        fix_low  = mpf("0")
        fix_high = inf
        sig_low  = mpf("0")
        sig_high = inf

    # Since I am ignorant of how to deal with locale-specific issues,
    # you can set these things yourself.
    decimal_point      = "."
    exponent_character = "e"
    implicit_plus_sign = True   # " 4" instead of "4" if True
    explicit_plus_sign = False  # "+4" instead of "4" if True
    # Set the class variable comma_decorate to true if you want
    # comma decoration.
    comma_decorate = False
    left_comma_character = ","   # For decoration to the left of decimal point
    left_comma_spacing = 3       # Set to zero to disable
    right_comma_character = " "  # For decoration to the right of decimal point
    right_comma_spacing = 5      # Set to zero to disable

    # Set this as you wish to have your integer exponents formatted.
    # For example, you can include sign and leading zeros.  If the
    # exponent won't fit in this specification, then it is ignored.
    exponent_format    = "%+d"

    # Set this to True to cuddle the SI letter:  then "12.3 k" will be
    # "12.3k".
    cuddle_si = False

    # Set this if you want exponents of zero to be shown in sci and eng
    show_zero_exponent = False

    def __init__(self, num_digits=4):
        self.digits(num_digits)

    def digits(self, num_digits):
        self.num_digits = num_digits
        if num_digits < 0:
            self.num_digits = 0 # 0 works for fix; others will use 1
        if num_digits > mp.dps:
            self.num_digits = mp.dps

    def format(self, number, format_type="sig"):
        '''Convenience function to allow you to use a string to specify the
        formatting.  format_type must be one of sig, fix, sci, eng, or engsi.
        '''
        if format_type == "sig":
            return self.sig(number)
        elif format_type == "fix":
            return self.fix(number)
        elif format_type == "sci":
            return self.sci(number)
        elif format_type == "eng":
            return self.eng(number)
        elif format_type == "engsi":
            return self.engsi(number)
        else:
            return str(number)

    def _pathological(self, number):
        '''Return a string if the number is not 'normal'.'''
        if not isinstance(number, mpf):
            raise ValueError("mpFormat._pathological():  expected mpf")
        if str(number) == "nan":
            return "NaN"
        elif str(number) == "+inf":
            return "+inf"
        elif str(number) == "-inf":
            return "-inf"
        else:
            return None

    def _to_estr(self, number, dps):
        '''Hacked version of libmpf.to_str.  Returns (sign, mantissa,
        exponent) where sign is "" or "-", mantissa is a string of one
        of the forms "d." or "d." with a number of digits after the
        decimal point, and exponent is an integer.
        '''
        s = number._mpf_
        if dps < 0:
            raise ValueError("_to_estr:  dps must be >= 0")
        dps = max(1, min(dps, mp.dps))
        if not s[1]:
            if s == fzero:
                if mpFormat.implicit_plus_sign: sign = " "
                if mpFormat.explicit_plus_sign: sign = "+"
                return sign, '0', 0
            if s == finf or s == fninf or s == fnan:
                raise ValueError("_to_estr:  should have caught pathology")
        sign, digits, exponent = to_digits_exp(s, dps+3)
        if sign == "" and mpFormat.implicit_plus_sign: sign = " "
        if (sign == "" or sign == " ") and \
           mpFormat.explicit_plus_sign:
            sign = "+"
        if not dps:
            if digits[0] in '56789':
                exponent += 1
            digits = ".0"
        else:
            if len(digits) > dps and digits[dps] in '56789' and \
                (dps < 500 or digits[dps-4:dps] == '9999'):
                digits2 = str(int(digits[:dps]) + 1)
                if len(digits2) > dps:
                    digits2 = digits2[:dps]
                    exponent += 1
                digits = digits2
            else:
                digits = digits[:dps]
            digits = (digits[:1] + "." + digits[1:])
        return sign, digits, int(exponent)

    def fix(self, number):
        s = self._pathological(number)
        if s: return s
        if number != 0 and abs(number) <= mpFormat.fix_low or \
                           abs(number) >= mpFormat.fix_high:
            return self.sci(number)
        sign, mant, exp = self._to_estr(number, mp.dps)
        digits = self.num_digits
        dp = mpFormat.decimal_point
        mant = mant[0] + mant[2:]  # Remove the decimal point
        if digits < mp.dps:
            # Round the mantissa if needed
            if mant != "0":
                last_digit = exp + digits + 1
                m = mant[:last_digit + 1]   # We go one beyond to round
                n = 0
                if m: n = int(m) + 5        # Do the rounding
                new_m = str(n)
                if m and len(new_m) > len(m):     # Rounding increased size
                    new_m = new_m[:-1]
                    exp += 1
                mant = str(new_m)
        if exp >= 0:    # Move the decimal point to the right
            mant = mant[:exp+1] + dp + mant[exp+1:]
        else:           # Move the decimal point to the left
            if exp < -1: mant = ("0"*(-exp-1)) + mant
            mant = "0" + dp +  mant
        loc_dp = mant.find(dp)
        digits_to_right_of_dp = len(mant) - loc_dp
        digits = self.num_digits + 1
        if digits_to_right_of_dp < digits:
            # Append zeros
            mant += "0"*(digits - digits_to_right_of_dp)
        elif digits_to_right_of_dp > digits:
            # Truncate, as we have too many digits
            mant = mant[:-(digits_to_right_of_dp - digits)]
        if mpFormat.comma_decorate:
            mant = self.decorate_with_comma(mant)
        return sign + mant

    def sci(self, number):
        s = self._pathological(number)
        if s: return s
        sign, mant, exp = self._to_estr(number, self.num_digits)
        s = sign + mant + mpFormat.exponent_character + \
            (mpFormat.exponent_format % exp)
        if not mpFormat.show_zero_exponent and exp == 0:
            return sign + mant
        return s

    def eng(self, number):
        'Engineering format for a floating point number'
        s = self._pathological(number)
        if s: return s
        mant, exp = self._eng(number)
        s = mant + mpFormat.exponent_character + \
            (mpFormat.exponent_format % exp)
        if not mpFormat.show_zero_exponent and exp == 0:
            return mant
        return s

    def engsi(self, number):
        'Same as eng(), but decorate with SI suffix.'
        s = self._pathological(number)
        if s: return s
        mant, exp = self._eng(number)
        if exp in suffixes_nl:
            spc = " "
            if mpFormat.cuddle_si: spc = ""
            return mant + spc + suffixes_nl[exp]
        else:
            return self.eng(number)

    def _eng(self, number):
        'Return (mant, exp)'
        digits = self.num_digits
        sign, mant, exp = self._to_estr(number, digits)
        num3, dp = divmod(exp, 3)
        mant = mant[0] + mant[2:]
        if digits > 0:
            dp += 1
            if len(mant) < dp:
                mant += "0"*max(0, dp - len(mant))
            mant = mant[:dp] + mpFormat.decimal_point + mant[dp:]
        else:
            mant += "0"*dp
        return sign + mant, 3*num3

    def sig(self, number):
        s = self._pathological(number)
        if s: return s
        if number != 0 and abs(number) <= mpFormat.sig_low or \
                           abs(number) >= mpFormat.sig_high:
            return self.sci(number)
        digits = self.num_digits
        sign, mant, exp = self._to_estr(number, digits)
        dp = mpFormat.decimal_point
        mant = mant[0] + mant[2:]  # Remove the decimal point
        exp += 1  # Now implied decimal point is at left of mantissa
        if exp < 0:
            mant = "0" + dp + ("0"*abs(exp)) + mant
        elif exp == 0:
            if digits:
                if len(mant) < digits:
                    mant += "0"*max(0, digits - len(mant))
            mant = "0" + dp + mant
        else:
            mant += "0"*max(0, exp - len(mant))
            if digits:
                if len(mant) < digits:
                    mant += "0"*max(0, digits - len(mant))
                mant = mant[:exp] + dp + mant[exp:]
        if mpFormat.comma_decorate:
            mant = self.decorate_with_comma(mant)
        return sign + mant

    def decorate_with_comma(self, mant):
        dp = mant.find(mpFormat.decimal_point)
        # Separate into two strings on either side of decimal point
        fp = mant[dp+1:]  # Fractional part
        ip = mant[:dp]    # Integer part
        if mpFormat.left_comma_spacing > 0 and \
           dp >= mpFormat.left_comma_spacing:
            ip = list(ip)
            ip.reverse()
            s = []
            for i, digit in enumerate(ip):
                if i != 0 and i % mpFormat.left_comma_spacing == 0:
                    s.append(mpFormat.left_comma_character)
                s.append(digit)
            s.reverse()
            ip = "".join(s)
        if mpFormat.right_comma_spacing > 0:
            fp = list(fp)
            s = []
            for i, digit in enumerate(fp):
                if i != 0 and i % mpFormat.right_comma_spacing == 0:
                    s.append(mpFormat.right_comma_character)
                s.append(digit)
            fp = "".join(s)
        return ip + mpFormat.decimal_point + fp

if __name__ == "__main__":
    # Test the mpFormat class.
    mp.dps = 30
    # Set the class variables
    f = mpFormat()
    mpFormat.decimal_point      = "."
    mpFormat.exponent_character = "e"
    mpFormat.show_zero_exponent = True
    mpFormat.exponent_format    = "%+d"
    mpFormat.fix_low            = mpf("0")
    mpFormat.fix_high           = inf
    mpFormat.sig_low            = mpf("0")
    mpFormat.sig_high           = inf
    def Test_sig():
        f = mpFormat()
        x = mpf("1.2345678901234567890")
        f.digits( 0); assert(f.sig(x)  == " 1")
        f.digits( 0); assert(f.sig(-x) == "-1")
        f.digits( 1); assert(f.sig(x)  == " 1.")
        f.digits( 1); assert(f.sig(-x) == "-1.")
        f.digits( 2); assert(f.sig(x)  == " 1.2")
        f.digits( 3); assert(f.sig(x)  == " 1.23")
        f.digits( 4); assert(f.sig(x)  == " 1.235")
        f.digits( 5); assert(f.sig(x)  == " 1.2346")
        f.digits(20); assert(f.sig(x)  == " 1.2345678901234567890")
        f.digits(20); assert(f.sig(-x) == "-1.2345678901234567890")
        x /= mpf("10")
        f.digits( 0); assert(f.sig(x)  == " 0.1")
        f.digits( 0); assert(f.sig(-x) == "-0.1")
        f.digits( 1); assert(f.sig(x)  == " 0.1")
        f.digits( 1); assert(f.sig(-x) == "-0.1")
        f.digits( 2); assert(f.sig(x)  == " 0.12")
        f.digits( 3); assert(f.sig(x)  == " 0.123")
        f.digits( 4); assert(f.sig(x)  == " 0.1235")
        f.digits( 5); assert(f.sig(x)  == " 0.12346")
        f.digits(20); assert(f.sig(x)  == " 0.12345678901234567890")
        f.digits(20); assert(f.sig(-x) == "-0.12345678901234567890")
        mpFormat.high = mpf("2e6")
        x *= mpf("1e7")
        f.digits( 0); assert(f.sig(x)  == " 1000000")
        f.digits( 0); assert(f.sig(-x) == "-1000000")
        f.digits( 1); assert(f.sig(x)  == " 1000000.")
        f.digits( 2); assert(f.sig(x)  == " 1200000.")
        f.digits( 3); assert(f.sig(x)  == " 1230000.")
        f.digits( 7); assert(f.sig(x)  == " 1234568.")
        f.digits( 8); assert(f.sig(x)  == " 1234567.9")
        f.digits(20); assert(f.sig(x)  == " 1234567.8901234567890")
        f.digits(20); assert(f.sig(-x) == "-1234567.8901234567890")
        x /= mpf("1e12")
        mpFormat.low = mpf("1e-7")
        f.digits( 0); assert(f.sig(x)  == " 0.000001")
        f.digits( 0); assert(f.sig(-x) == "-0.000001")
        f.digits( 1); assert(f.sig(x)  == " 0.000001")
        f.digits( 2); assert(f.sig(x)  == " 0.0000012")
        f.digits( 6); assert(f.sig(x)  == " 0.00000123457")
        f.digits(20); assert(f.sig(x)  == " 0.0000012345678901234567890")
        f.digits(20); assert(f.sig(-x) == "-0.0000012345678901234567890")
    def Test_fix():
        f = mpFormat()
        x = mpf("1.2345678901234567890")
        f.digits( 0); assert(f.fix(x)  == " 1.")
        f.digits( 1); assert(f.fix(x)  == " 1.2")
        f.digits( 2); assert(f.fix(x)  == " 1.23")
        f.digits( 3); assert(f.fix(x)  == " 1.235")
        f.digits(20); assert(f.fix(x)  == " 1.23456789012345678900")
        f.digits(20); assert(f.fix(-x) == "-1.23456789012345678900")
        x *= mpf("1e9")
        f.digits( 0); assert(f.fix(x) == " 1234567890.")
        f.digits( 1); assert(f.fix(x) == " 1234567890.1")
        f.digits( 2); assert(f.fix(x) == " 1234567890.12")
        f.digits( 3); assert(f.fix(x) == " 1234567890.123")
        f.digits( 4); assert(f.fix(x) == " 1234567890.1235")
        f.digits(10); assert(f.fix(x) == " 1234567890.1234567890")
        x /= mpf("1e12")
        f.digits( 0); assert(f.fix(x)  == " 0.")
        f.digits( 0); assert(f.fix(-x) == "-0.")
        f.digits( 1); assert(f.fix(x)  == " 0.0")
        f.digits( 1); assert(f.fix(-x) == "-0.0")
        f.digits( 2); assert(f.fix(x)  == " 0.00")
        f.digits( 2); assert(f.fix(-x) == "-0.00")
        f.digits( 3); assert(f.fix(x)  == " 0.001")
        f.digits( 3); assert(f.fix(-x) == "-0.001")
        f.digits( 4); assert(f.fix(x)  == " 0.0012")
        f.digits( 5); assert(f.fix(x)  == " 0.00123")
        f.digits( 7); assert(f.fix(x)  == " 0.0012346")
        f.digits(25); assert(f.fix(x)  == " 0.0012345678901234567890000")
        f.digits(25); assert(f.fix(-x) == "-0.0012345678901234567890000")
    def Test_sci():
        f = mpFormat()
        mpFormat.show_zero_exponent = True
        x = mpf("1.2345678901234567890")
        f.digits( 0); assert(f.sci( x) == " 1.e+0")
        f.digits( 0); assert(f.sci(-x) == "-1.e+0")
        f.digits( 1); assert(f.sci(x)  == " 1.e+0")
        f.digits( 1); assert(f.sci(-x) == "-1.e+0")
        f.digits( 2); assert(f.sci(x)  == " 1.2e+0")
        f.digits( 3); assert(f.sci(x)  == " 1.23e+0")
        f.digits( 3); assert(f.sci(-x) == "-1.23e+0")
        f.digits(20); assert(f.sci(x)  == " 1.2345678901234567890e+0")
        mpFormat.show_zero_exponent = False
        f.digits( 0); assert(f.sci( x) == " 1.")
        f.digits( 0); assert(f.sci(-x) == "-1.")
        f.digits( 1); assert(f.sci(x)  == " 1.")
        f.digits( 1); assert(f.sci(-x) == "-1.")
        f.digits( 2); assert(f.sci(x)  == " 1.2")
        f.digits( 3); assert(f.sci(x)  == " 1.23")
        f.digits( 3); assert(f.sci(-x) == "-1.23")
        f.digits(20); assert(f.sci(x)  == " 1.2345678901234567890")
        x *= mpf("1e9")
        f.digits( 0); assert(f.sci(x)  == " 1.e+9")
        f.digits( 0); assert(f.sci(-x) == "-1.e+9")
        f.digits( 1); assert(f.sci(x)  == " 1.e+9")
        f.digits( 2); assert(f.sci(x)  == " 1.2e+9")
        f.digits( 2); assert(f.sci(-x) == "-1.2e+9")
        x /= mpf("1e18")
        f.digits( 0); assert(f.sci(x)  == " 1.e-9")
        f.digits( 1); assert(f.sci(x)  == " 1.e-9")
        f.digits( 2); assert(f.sci(x)  == " 1.2e-9")
        f.digits(20); assert(f.sci( x) == " 1.2345678901234567890e-9")
        f.digits(20); assert(f.sci(-x) == "-1.2345678901234567890e-9")
    def Test_eng():
        f = mpFormat()
        mpFormat.show_zero_exponent = True
        x = mpf("1.2345678901234567890")
        f.digits( 0); assert(f.eng( x) == " 1e+0")
        f.digits( 0); assert(f.eng(-x) == "-1e+0")
        f.digits( 1); assert(f.eng(x)  == " 1.e+0")
        f.digits( 1); assert(f.eng(-x) == "-1.e+0")
        f.digits( 2); assert(f.eng(x)  == " 1.2e+0")
        f.digits( 3); assert(f.eng(x)  == " 1.23e+0")
        f.digits( 3); assert(f.eng(-x) == "-1.23e+0")
        x = mpf("1.2345678901234567890e1")
        f.digits( 0); assert(f.eng( x) == " 10e+0")
        f.digits( 0); assert(f.eng(-x) == "-10e+0")
        f.digits( 1); assert(f.eng(x)  == " 10.e+0")
        f.digits( 1); assert(f.eng(-x) == "-10.e+0")
        f.digits( 2); assert(f.eng(x)  == " 12.e+0")
        f.digits( 3); assert(f.eng(x)  == " 12.3e+0")
        f.digits( 3); assert(f.eng(-x) == "-12.3e+0")
        x = mpf("1.2345678901234567890e-1")
        f.digits( 0); assert(f.eng( x) == " 100e-3")
        f.digits( 0); assert(f.eng(-x) == "-100e-3")
        f.digits( 1); assert(f.eng(x)  == " 100.e-3")
        f.digits( 1); assert(f.eng(-x) == "-100.e-3")
        f.digits( 2); assert(f.eng(x)  == " 120.e-3")
        f.digits( 2); assert(f.eng(-x) == "-120.e-3")
        f.digits( 3); assert(f.eng(x)  == " 123.e-3")
        f.digits( 3); assert(f.eng(-x) == "-123.e-3")
        x = mpf("1.2345678901234567890e2")
        f.digits( 0); assert(f.eng( x) == " 100e+0")
        f.digits( 0); assert(f.eng(-x) == "-100e+0")
        f.digits( 1); assert(f.eng(x)  == " 100.e+0")
        f.digits( 1); assert(f.eng(-x) == "-100.e+0")
        f.digits( 2); assert(f.eng(x)  == " 120.e+0")
        f.digits( 2); assert(f.eng(-x) == "-120.e+0")
        f.digits( 3); assert(f.eng(x)  == " 123.e+0")
        f.digits( 3); assert(f.eng(-x) == "-123.e+0")
        x = mpf("1.2345678901234567890e-2")
        f.digits( 0); assert(f.eng( x) == " 10e-3")
        f.digits( 0); assert(f.eng(-x) == "-10e-3")
        f.digits( 1); assert(f.eng(x)  == " 10.e-3")
        f.digits( 1); assert(f.eng(-x) == "-10.e-3")
        f.digits( 2); assert(f.eng(x)  == " 12.e-3")
        f.digits( 3); assert(f.eng(x)  == " 12.3e-3")
        f.digits( 3); assert(f.eng(-x) == "-12.3e-3")
        x = mpf("1.2345678901234567890e3")
        f.digits( 0); assert(f.eng( x) == " 1e+3")
        f.digits( 0); assert(f.eng(-x) == "-1e+3")
        f.digits( 1); assert(f.eng(x)  == " 1.e+3")
        f.digits( 1); assert(f.eng(-x) == "-1.e+3")
        f.digits( 2); assert(f.eng(x)  == " 1.2e+3")
        f.digits( 3); assert(f.eng(x)  == " 1.23e+3")
        f.digits( 3); assert(f.eng(-x) == "-1.23e+3")
        x = mpf("1.2345678901234567890e-3")
        f.digits( 0); assert(f.eng( x) == " 1e-3")
        f.digits( 0); assert(f.eng(-x) == "-1e-3")
        f.digits( 1); assert(f.eng(x)  == " 1.e-3")
        f.digits( 1); assert(f.eng(-x) == "-1.e-3")
        f.digits( 2); assert(f.eng(x)  == " 1.2e-3")
        f.digits( 3); assert(f.eng(x)  == " 1.23e-3")
        f.digits( 3); assert(f.eng(-x) == "-1.23e-3")
    def Test_engsi():
        f = mpFormat()
        f.digits(2)
        x = mpf(1.23e-27); assert(f.engsi(x) == " 1.2e-27")
        x = mpf(1.23e-26); assert(f.engsi(x) == " 12.e-27")
        x = mpf(1.23e-25); assert(f.engsi(x) == " 120.e-27")
        x = mpf(1.23e-24); assert(f.engsi(x) == " 1.2 y")
        x = mpf(1.23e-23); assert(f.engsi(x) == " 12. y")
        x = mpf(1.23e-22); assert(f.engsi(x) == " 120. y")
        x = mpf(1.23e-21); assert(f.engsi(x) == " 1.2 z")
        x = mpf(1.23e-20); assert(f.engsi(x) == " 12. z")
        x = mpf(1.23e-19); assert(f.engsi(x) == " 120. z")
        x = mpf(1.23e-18); assert(f.engsi(x) == " 1.2 a")
        x = mpf(1.23e-17); assert(f.engsi(x) == " 12. a")
        x = mpf(1.23e-16); assert(f.engsi(x) == " 120. a")
        x = mpf(1.23e-15); assert(f.engsi(x) == " 1.2 f")
        x = mpf(1.23e-14); assert(f.engsi(x) == " 12. f")
        x = mpf(1.23e-13); assert(f.engsi(x) == " 120. f")
        x = mpf(1.23e-12); assert(f.engsi(x) == " 1.2 p")
        x = mpf(1.23e-11); assert(f.engsi(x) == " 12. p")
        x = mpf(1.23e-10); assert(f.engsi(x) == " 120. p")
        x = mpf(1.23e-09); assert(f.engsi(x) == " 1.2 n")
        x = mpf(1.23e-08); assert(f.engsi(x) == " 12. n")
        x = mpf(1.23e-07); assert(f.engsi(x) == " 120. n")
        x = mpf(1.23e-06); assert(f.engsi(x) == " 1.2 u")
        x = mpf(1.23e-05); assert(f.engsi(x) == " 12. u")
        x = mpf(1.23e-04); assert(f.engsi(x) == " 120. u")
        x = mpf(1.23e-03); assert(f.engsi(x) == " 1.2 m")
        x = mpf(1.23e-02); assert(f.engsi(x) == " 12. m")
        x = mpf(1.23e-01); assert(f.engsi(x) == " 120. m")
        x = mpf(1.23e+00); assert(f.engsi(x) == " 1.2 ")
        x = mpf(1.23e+01); assert(f.engsi(x) == " 12. ")
        x = mpf(1.23e+02); assert(f.engsi(x) == " 120. ")
        x = mpf(1.23e+03); assert(f.engsi(x) == " 1.2 k")
        x = mpf(1.23e+04); assert(f.engsi(x) == " 12. k")
        x = mpf(1.23e+05); assert(f.engsi(x) == " 120. k")
        x = mpf(1.23e+06); assert(f.engsi(x) == " 1.2 M")
        x = mpf(1.23e+07); assert(f.engsi(x) == " 12. M")
        x = mpf(1.23e+08); assert(f.engsi(x) == " 120. M")
        x = mpf(1.23e+09); assert(f.engsi(x) == " 1.2 G")
        x = mpf(1.23e+10); assert(f.engsi(x) == " 12. G")
        x = mpf(1.23e+11); assert(f.engsi(x) == " 120. G")
        x = mpf(1.23e+12); assert(f.engsi(x) == " 1.2 T")
        x = mpf(1.23e+13); assert(f.engsi(x) == " 12. T")
        x = mpf(1.23e+14); assert(f.engsi(x) == " 120. T")
        x = mpf(1.23e+15); assert(f.engsi(x) == " 1.2 P")
        x = mpf(1.23e+16); assert(f.engsi(x) == " 12. P")
        x = mpf(1.23e+17); assert(f.engsi(x) == " 120. P")
        x = mpf(1.23e+18); assert(f.engsi(x) == " 1.2 E")
        x = mpf(1.23e+19); assert(f.engsi(x) == " 12. E")
        x = mpf(1.23e+20); assert(f.engsi(x) == " 120. E")
        x = mpf(1.23e+21); assert(f.engsi(x) == " 1.2 Z")
        x = mpf(1.23e+22); assert(f.engsi(x) == " 12. Z")
        x = mpf(1.23e+23); assert(f.engsi(x) == " 120. Z")
        x = mpf(1.23e+24); assert(f.engsi(x) == " 1.2 Y")
        x = mpf(1.23e+25); assert(f.engsi(x) == " 12. Y")
        x = mpf(1.23e+26); assert(f.engsi(x) == " 120. Y")
        x = mpf(1.23e+27); assert(f.engsi(x) == " 1.2e+27")
    def Test_comma_decoration():
        mp.dps = 100
        f = mpFormat()
        f.digits(11)
        mpFormat.fix_high = inf
        mpFormat.sig_high = inf
        mpFormat.decimal_point = "."
        mpFormat.left_comma_character = ","
        mpFormat.right_comma_character = " "
        mpFormat.left_comma_spacing = 3
        mpFormat.right_comma_spacing = 5
        from mpmath import pi
        x = pi*1e7
        mpFormat.comma_decorate = True
        assert f.fix(x) == " 31,415,926.53589 79323 8"
        f.digits(15)
        assert f.sig(x) == " 31,415,926.53589 79"
    Test_sig()
    Test_fix()
    Test_sci()
    Test_eng()
    Test_engsi()
    Test_comma_decoration()
