EnsurePythonVersion(2,7)

import os

PROJ_NAME = 'flatradix'
PROJ_SRC = 'src'

AddOption('--out',
    nargs=1,
    dest='outdir',
    action='store',
    default='bin',
    help='directory to use as the variant directory'
)
AddOption('--clang',
    action='store_true',
    help='use the clang/LLVM toolset rather than the GCC',
)
AddOption('--dbg',
    action='store_true',
    help='compile binary/libraries in debug mode (-O0, non-stripped)'
)
AddOption('--ugly',
    action='store_true',
    help='print the entire file-compilation string rather than pretty-printing'
)
AddOption('--valgrind',
    action='store_true',
    help='run the tests under valgrind'
)

env = Environment()

variant_dir = os.path.abspath(GetOption('outdir'))
SConscript('SConscript', variant_dir=variant_dir, duplicate=False)
Clean('.', variant_dir)
