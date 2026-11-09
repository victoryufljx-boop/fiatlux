"""Convert qlanth (zia-lab/qlanth @14563fd) symbolic LS tables into numeric JSON for n<=NMAX."""
import re, json, gzip, sys, math
from fractions import Fraction
src = sys.argv[1]; out = sys.argv[2]; NMAX = 4
P = ["F0","F2","F4","F6","zeta","M0","M2","M4","P2","P4","P6","T2","T2p","T3","T4","T6","T7","T8","T11","T11p","T12","T14","T15","T16","T17","T18","T19"]
KEY = re.compile(r'\{(\d+), (?:3, )?"(\w+)", "(\w+)"(?:, (\d+(?:/\d+)?))?\} -> ')
class V(dict):
    def __add__(s,o):
        r=V(s)
        if isinstance(o,V):
            for k,v in o.items(): r[k]=r.get(k,0)+v
            return r
        if o==0: return r
        raise ValueError(o)
    __radd__=__add__
    def __neg__(s): return V({k:-v for k,v in s.items()})
    def __sub__(s,o): return s+(-o)
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if isinstance(o,V): raise ValueError
        return V({k:v*o for k,v in s.items()})
    __rmul__=__mul__
    def __truediv__(s,o): return V({k:v/o for k,v in s.items()})
    def __pos__(s): return s
def parse(fname, t2=None, nkey=None):
    s = open(fname).read()
    s = s[s.find('<|')+2:]
    ms = list(KEY.finditer(s))
    env = {p: V({p:1.0}) for p in P}
    env.update(sqrt=math.sqrt, zeta=V({"zeta":1.0}), t2Switch=t2)
    res = {}
    for a, b in zip(ms, ms[1:]+[None]):
        n = int(a.group(1))
        if n > NMAX: continue
        e = s[a.end(): (b.start() if b else len(s))].strip().rstrip('>').rstrip('|').strip().rstrip(',')
        e = e.replace('\\[Zeta]','zeta').replace('Sqrt[','sqrt(').replace(']',')').replace('\n',' ').replace('^','**')
        val = eval(e, {"__builtins__":{}}, env)
        if isinstance(val, V): val = {k:v for k,v in val.items() if abs(v)>1e-14}
        else: val = float(val)
        if (isinstance(val,dict) and not val) or (not isinstance(val,dict) and val==0): continue
        J = a.group(4)
        res[f"{n}|{a.group(2)}|{a.group(3)}" + (f"|{float(Fraction(J))}" if J else "")] = val
    return res
d = {}
d["electrostatic"] = parse(f"{src}/ElectrostaticTable.m")
d["spin_orbit"] = parse(f"{src}/SpinOrbitTable.m")
d["soo_ecso"] = parse(f"{src}/SOOandECSOTable.m")
d["spin_spin"] = parse(f"{src}/SpinSpinTable.m")
d["three_body_t2switch1"] = parse(f"{src}/ThreeBodyTable.m", t2=1)
d["three_body_t2switch0"] = parse(f"{src}/ThreeBodyTable.m", t2=0)
# reduced U^k, key n|SL|SL' |k (k stored in J slot)
d["Uk"] = parse(f"{src}/ReducedUkTable.m")
terms = re.findall(r'\{"(\w+)", (\d+), "(\d+)", "(\d+)"\}', open(f"{src}/fnTerms.m").read())
d["terms_NK"] = terms
d["src"] = "zia-lab/qlanth commit 14563fd6c91b164dab47a5b3ff6f3686ea33a44b, data/*.m (Electrostatic, SpinOrbit, SOOandECSO (P prescaled), SpinSpin, ThreeBody, ReducedUk, fnTerms)"
for k,v in d.items(): print(k, len(v))
json.dump(d, gzip.open(out,"wt"))
