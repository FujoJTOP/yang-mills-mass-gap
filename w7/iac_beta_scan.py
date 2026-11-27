import io, sys, math
p = r"research\v3\species\_ht7\faizal.txt"
lines = io.open(p, encoding="utf-8", errors="replace").read().split("\n")
print("== 23280..23300 ==")
for i in range(23280, 23301):
    print(i + 1, lines[i])
print("== beta-star hits ==")
for i, l in enumerate(lines):
    if ("\u03b2 \u22c6" in l) or ("\u03b2\u22c6" in l) or ("\u03b2 \u2217" in l):
        print(i + 1, l.strip()[:170])
print("== key arithmetic ==")
e = math.e
print("1/(2e^2*12)      %.10f" % (1 / (2 * e * e * 12)))
print("1/(2e^2*20)      %.10f" % (1 / (2 * e * e * 20)))
print("e^-2/2.5         %.10f" % ((1 / (e * e)) / 2.5))
print("1/(2e^2*12)/2.5  %.10f" % ((1 / (2 * e * e * 12)) / 2.5))
print("1/(2e^2*20)/2.5  %.10f" % ((1 / (2 * e * e * 20)) / 2.5))
print("2.5*e^2          %.10f" % (2.5 * e * e))
print("0.625*e^2        %.10f" % (0.625 * e * e))
print("(1/e^2)/(1/(2e^2*12)) %.6f" % ((1 / (e * e)) / (1 / (2 * e * e * 12))))
print("2e^65.2388       %.3e" % (2 * math.exp(2 * e * 12) * e * 12))
