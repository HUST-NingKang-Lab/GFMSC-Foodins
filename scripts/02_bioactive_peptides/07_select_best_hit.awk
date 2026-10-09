BEGIN {
    FS = OFS = "\t"
}

!($1 in best_evalue) || ($15 + 0) < best_evalue[$1] {
    best_evalue[$1] = $15 + 0
    best_line[$1] = $0
}

END {
    for (query in best_line) {
        print best_line[query]
    }
}
