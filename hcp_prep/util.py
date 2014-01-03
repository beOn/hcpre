def numberfy(s):
    n = s
    try:
        n = float(n)
        return n
    except Exception:
        return s

def float_or_none(s):
    n = s
    try:
        n = float(n)
        return n
    except Exception:
        return None

def int_or_none(s):
    n = s
    try:
        n = int(n)
        return n
    except ValueError:
        return None