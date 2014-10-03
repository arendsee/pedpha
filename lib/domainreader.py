def _read_data(self, data, delimiter):
    out = collections.defaultdict(list)
    for line in data:
        row = line.split(delimiter)
        try:
            seqid, domid = row[0:2]
            start, stop = (int(s) for s in row[2:4])
            bounds = to_dna_interval(sorted([start, stop]))
            value = (domid, bounds)
            out[seqid].append(value)
            if start < 1 or stop < 1:
                raise ValueError
        except IndexError:
            sys.exit("Each interval line must have 4 columns")
        except ValueError:
            sys.exit("Interval coordinants must be integers greater than 0")
    return(out)

