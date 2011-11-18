#!/usr/bin/env python

import sys

try:
    import lhc.hc as hc
except ImportError as e:
    print "ERROR:",e 
    print "Your python path includes:"
    print sys.path
    print "Ensure the lhc package is included in the path and try again."
    sys.exit(1)

if __name__ == "__main__":
	hc.main(sys.argv)
