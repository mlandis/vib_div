in_fp = "/Users/mlandis/projects/gh_vib_div/output_raw/"
out_fp = "/Users/mlandis/projects/gh_vib_div/output/"
files = list.files(path=fp)

n_keep = 2001

#files = files[1]
for (fn in files) {
    #print(fn)
    x = read.table(paste0(in_fp,fn), sep="\t", header=T, stringsAsFactors=F)
    it = unique(x$Iteration)
    n_it = length(it)
    x = x[(n_it-n_keep+1):n_it, ]
    x$Iteration = seq(0,100000,by=50)
    #print(x)
    write.table( x, paste0(out_fp,fn), sep="\t", quote=F)
}