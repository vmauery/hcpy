#!/usr/bin/env python3

import sys

try:
    import lhc.hc as hc
except ImportError as e:
    print("ERROR:",e) 
    print("Your python path includes:")
    print(sys.path)
    print("Ensure the lhc package is included in the path and try again.")
    sys.exit(1)

if __name__ == "__main__":
    if '-g' in sys.argv:
        hc.main(sys.argv)
    else:
        try:
            hc.main(sys.argv)
        except:
            pass

