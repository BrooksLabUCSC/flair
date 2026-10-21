# Reading and writing the TSVs the flair R scripts share, with the same suffix
# driven compression the Python side uses: a path ending .gz or .bz2 is
# compressed, anything else is plain.
#
# Reading needs nothing added.  R's file() and data.table::fread both detect
# compression from the file itself, so read.table and fread already accept a
# compressed path.
#
# Writing does need it.  write.table given a file name opens a plain file
# whatever the name ends in, so a table written to x.tsv.gz would be an
# uncompressed file wearing a .gz name.  The connection has to be opened here
# instead.

# what --compress appends; gzip is the one flair writes, though the writer below
# honours whatever a path it is handed already ends in
tsvCompressSuffix <- ".gz"

tsvCompressedPath <- function(path, compress) {
    # the output name a run writes, given whether --compress was asked for
    if (isTRUE(compress)) {
        paste0(path, tsvCompressSuffix)
    } else {
        path
    }
}

tsvWriteConnection <- function(path) {
    if (grepl("\\.gz$", path)) {
        gzfile(path, "w")
    } else if (grepl("\\.bz2$", path)) {
        bzfile(path, "w")
    } else if (grepl("\\.bgz$", path)) {
        # gzfile would write ordinary gzip under a name that promises the block
        # format, which is worse than refusing
        stop(sprintf("cannot write bgzf: %s", path))
    } else {
        file(path, "w")
    }
}

tsvWriteTable <- function(x, path, ...) {
    # write.table with the tab separated, unquoted conventions flair uses, through
    # a connection that compresses when the name says to
    con <- tsvWriteConnection(path)
    on.exit(close(con))
    write.table(x, file = con, sep = "\t", quote = FALSE, ...)
}

tsvFread <- function(path, ...) {
    # data.table's fread needs the R.utils package to open a compressed file
    # itself, so a compressed one is decompressed through a pipe instead, the way
    # the Python side does it
    if (grepl("\\.(gz|bgz)$", path)) {
        data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(path)), ...)
    } else if (grepl("\\.bz2$", path)) {
        data.table::fread(cmd = sprintf("bzip2 -dc %s", shQuote(path)), ...)
    } else {
        data.table::fread(path, ...)
    }
}
