import math


def correlation(calced, experimental):
    """Calculates correlation between calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    m = [0.0, 0.0, 0.0]
    d = [0.0, 0.0]
    match_count = 0

    for i in experimental:
        exp = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        m[0] += calc
        m[1] += exp
        m[2] += calc * exp

        match_count += 1

    m[0] /= match_count
    m[1] /= match_count
    m[2] /= match_count

    for i in experimental:
        exp = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        d[0] += (calc - m[0]) ** 2
        d[1] += (exp - m[1]) ** 2

    d[0] /= match_count
    d[0] = math.sqrt(d[0])
    d[1] /= match_count
    d[1] = math.sqrt(d[1])

    if d[0] * d[1] == 0:
        return -2
    else:
        return (m[2] - (m[0] * m[1])) / (d[0] * d[1])


def q_value(calced, experimental):
    """Calculates Q-value for calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    d2, e2 = 0, 0

    for i in experimental:
        exp = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        d2 += (calc - exp) ** 2
        e2 += exp ** 2

    q = 100 * math.sqrt(d2) / math.sqrt(e2)
    return round(q, 6)


def rmsd(calced, experimental):
    """Calculates RMSD for calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    d2 = 0
    match_count = 0

    for i in experimental:
        exp = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        d2 += (calc - exp) ** 2

        match_count += 1

    rmsd = math.sqrt(d2 / match_count)
    return round(rmsd, 6)
