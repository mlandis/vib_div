import os
import pandas as pd

# IO
fp = '/Users/mlandis/projects/gh_vib_div/'
in_fp = fp + 'output_raw/'
out_fp = fp + 'output/'

# get files
files = []
for r,d,f in os.walk(in_fp):
    for fn in f:
        files.append(fn)
files = [ x for x in files if 'out' in x ]

# process settings
n_it = 2001
it = [ x*50 for x in range(0,n_it) ]
#print(it)

# process files
#files = [ files[0] ]
#files = [ x for x in files if 'model' in x ]
for fn in files:
    in_fn = in_fp + fn
    df = pd.read_csv(in_fn, sep='\t', dtype=str)
    print(in_fn)
    #n_row = df.shape[0]
    print(df)
    df = df[-n_it:]
    df = df.assign(Iteration=it)
    print(df)
    print(df)
    df.to_csv(out_fp + fn, sep='\t', index=False)

