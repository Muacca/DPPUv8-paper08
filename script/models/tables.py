"""Readable Markdown tables generated from flat.json and nil.json (data/flat_thresholds.md, data/nil_thresholds.md)."""
from flint import arb

FLAT_HEAD = (
    "# 全28リフトの基準葉における最初の非零電荷\n\n"
    "各セルは正・負 $`Q_\\Sigma`$ の順。表記 $`v[d,O]`$ は $`v=(|\\lambda|/(2\\pi))^2`$ 、複素核次元 $`d`$ 、起源 $`O`$ "
    "（A は軸、H は非零横運動量）。したがって総電荷の絶対値は $`8\\pi V\\sqrt v`$ 。表の臨界値は二乗で、開いた分裂域は "
    "$`L^2>L_c^2`$ 。記号 — は全 $`L`$ で正負が等しいことを表す。スピンの順序と向きは付録Eの定義に従う。"
    "この表は [flat.json](flat.json) から `script/paper08.py recompute` が生成する。\n\n"
    "| 商 | spin | η | LC核 | L=1/2 | L=1 | L=2 | 分裂の $`L_c^2`$ | 両側が軸に達する $`L^2`$ |\n"
    "|---|---|---:|---:|---|---|---|---:|---:|\n"
)
NIL_HEAD = (
    "# 静的72断面のη・電荷閾値\n\n"
    "スピン欄は $`2δ_1,2δ_2,2δ_3`$ 。閾値は $`|Q|/(6V_Σ)`$ 。各欄の括弧に複素核次元と由来を記す。Pは対Hermiteの正枝、"
    "Uは非対系列、添字は $`k=r|m+δ_3|`$ 。Bの添字は $`4|l|^2`$ 、末尾は水平二重項の枝（B0の二枝は同一固有値）。"
    "丸め前のArb包含と総電荷は [nil.json](nil.json)。この表は `script/paper08.py recompute` が生成する。\n\n"
    "|r|2δ|t|η|非零中心 Q正|非零中心 Q負|全体 Q正|全体 Q負|全体 零電荷核|\n"
    "|---:|---|---:|---:|---|---|---|---|---:|\n"
)


def _mid(s):
    return float(arb(s).mid())


def flat_table(spins):
    rows = []
    for r in spins:
        cells = []
        for L in ("1/2", "1", "2"):
            pos, neg = r["thresholds"][L]
            cells.append(" ／ ".join(f"{c['lambda_squared_in_2pi_units']}[{c['dim']},{c['origin']}]" for c in (pos, neg)))
        crit = r["critical"]
        rows.append(f"| {r['group']} | ({', '.join(map(str, r['spin']))}) | {r['eta']} | {r['lc_kernel']} | "
                    + " | ".join(cells)
                    + f" | {crit.get('split_L_squared', '—')} | {crit.get('both_axial_L_squared', '—')} |")
    return FLAT_HEAD + "\n".join(rows) + "\n"


def nil_table(leaves):
    def cell(c):
        return f"{_mid(c['abs_u']):.7f} ({c['dim']}, {c['mode']})"

    rows = []
    for r in leaves:
        z, f = r["nonzero"], r["full"]
        rows.append(f"|{r['r']}|{''.join(map(str, r['spin2']))}|{r['t']}|{_mid(r['eta']):.8f}|"
                    f"{cell(z['Q_positive'])}|{cell(z['Q_negative'])}|{cell(f['Q_positive'])}|{cell(f['Q_negative'])}|"
                    f"{f['zero_kernel']}|")
    return NIL_HEAD + "\n".join(rows) + "\n"
