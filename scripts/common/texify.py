import re
from functools import singledispatch

@singledispatch
def texify(x):
    return x


lookup = {
    # Variables
    "QGPV":         r"$q$",
    "<A>cos(φ)":    r"$\langle A \rangle \cos(\phi)$",
    "F1+F2+F3":     r"$\langle F_\lambda \rangle$",
    # Units
    "K":            r"$\mathrm{K}$",
    "km":           r"$\mathrm{km}$",
    "1/s":          r"$\mathrm{s}^{-1}$",
    "m/s":          r"$\mathrm{m} \; \mathrm{s}^{-1}$",
    "m²/s²":        r"$\mathrm{m}^2 \; \mathrm{s}^{-2}$",
    "K/m²/s/kg":    r"$\mathrm{K} \; \mathrm{m}^{-2} \; \mathrm{s}^{-1} \; \mathrm{kg}^{-1}$",
}

valunit = re.compile("^(.+)\s(km|K)$")

@texify.register
def _(x: str):
    if x in lookup:
        return lookup[x]
    # Good enough for flux terms
    if "F1" in x or "F2" in x or "F3" in x:
        return "$" + x.replace("F1", "F_1").replace("F2", "F_2").replace("F3", "F_3") + "$"
    # Vertical levels are already combined, try to sort it out with a regex
    m = valunit.match(x)
    if m is not None:
        val, unit = m.groups()
        return r"${} \; {}$".format(val, texify(unit).strip("$"))
    # Nothing works, return what was given
    return x

@texify.register
def _(x: dict):
    return { k: texify(v) for k, v in x.items() }

