#!/usr/bin/env python
import sys
import imp

def main(argv=None):
    argv = sys.argv
    tool = argv[1]
    (file, pathname, description) = imp.find_module(tool)
    module = imp.load_module(tool, file, pathname, description)
    del sys.argv[0]
    module.main(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
