'''
$Id: si.py 1.4 2009/02/09 17:04:00 donp Exp $

Provides mappings from an integer exponent to an SI unit prefix (but
here we call them suffixes).

Also provides the inverse mapping.
'''

# SI suffixes, number to letter
suffixes_nl = { -24:"y", -21:"z", -18:"a", -15:"f", -12:"p", -9:"n",
                 -6:"u",  -3:"m",   0:"",    3:"k",   6:"M",  9:"G",
                 12:"T",  15:"P",  18:"E",  21:"Z",  24:"Y"}

# SI suffixes, letter to number
suffixes_ln = dict([[val, key] for key, val in suffixes_nl.items()])

if __name__ == "__main__":
    from sys import argv
    if len(argv) > 1:
        print "Number to letter"
        for i in suffixes_nl.items(): print "  ", i
        print "Letter to number"
        for i in suffixes_ln.items(): print "  ", i
    assert len(suffixes_nl) == len(suffixes_ln)
