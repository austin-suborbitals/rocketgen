EnsurePythonVersion(2,7)

import os
import subprocess

PROJ_NAME = 'rocketgen'
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


# initialize submodules
get_subs_cmd = [
    ['git', 'config', '--file', os.path.join(Dir('#').abspath, '.gitmodules'), '--name-only', '--get-regexp', 'path'],
    ['sed', "s/^submodule\.//"],
    ['sed', "s/\.path$//"],
]
get_subs = subprocess.Popen(get_subs_cmd[0], stdout=subprocess.PIPE)
filter_sub = subprocess.Popen(get_subs_cmd[1], stdin=get_subs.stdout, stdout=subprocess.PIPE)
get_subs.wait()
filter_sub.wait()
submodules = filter(None, subprocess.check_output(get_subs_cmd[2], stdin=filter_sub.stdout).split('\n'))
for mod in submodules:
    if not len([f for f in os.listdir(mod) if os.path.isfile(os.path.join(mod, f))]):
        subprocess.check_call(['git', 'submodule', 'init'])
        subprocess.check_call(['git', 'submodule', 'update', '--depth', '1'])
        break

env = Environment()

variant_dir = os.path.abspath(GetOption('outdir'))
SConscript('SConscript', variant_dir=variant_dir, duplicate=False)
Clean('.', variant_dir)
