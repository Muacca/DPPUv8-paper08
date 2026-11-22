# 第2節　電荷、スピンの規約、スペクトル量

各コンパクトな空間切片 $`Y`$ に向き $`e^{123}>0`$ を与え、被覆変換群のスピンリフトを一つ指定する。四次元計量の符号は $`(-+++)`$ とし、 $`\epsilon_{0123}=+1`$ 、 $`T^a=De^a`$ を用いる。電荷と共役運動量を

```math
Q_\Sigma=\int_Y e_a\wedge T^a,\qquad p_\theta=-Q_\Sigma
```

で定義する。作用におけるねじれとの結合は $`\int\theta\,d(e_a\wedge T^a)`$ であり、部分積分すると、接続に対するバルクの源は $`-d\theta\wedge e_a\wedge T^a`$ となる。誘導余標構を固定したとき、場を固定する変分問題から電荷を固定する変分問題へ移るには、端点項 $`-[p_\theta\theta]_i^f=[\theta Q_\Sigma]_i^f`$ を加える。実定数の Holst 係数 $`c_H`$ を含める場合は $`h=2\kappa^2c_H`$ 、 $`d=1+h^2`$ とおく。接続を消去すると、 $`\phi_h=\sqrt{6}\,\kappa\theta/\sqrt{d}`$ に対する有効運動項 $`-\tfrac12(\partial\phi_h)^2`$ を得る。最小結合模型では $`h=0`$ である。斉次な軸性コントーションを $`K_{ij}=w\,\epsilon_{ijk}e^k`$ と書くと、固有時 $`\tau`$ において次式が成り立つ。

<a id="eq-shift"></a>

```math
\begin{aligned}
w&=\frac{\kappa^2}{d}\frac{d\theta}{d\tau},\qquad Q_\Sigma=-6\mathcal V w,
\\ D(Q_\Sigma)&=D_{\mathrm{LC}}+\tfrac32 w
           =D_{\mathrm{LC}}-\frac{Q_\Sigma}{4\mathcal V}.
\end{aligned}
\qquad(1)
```

ここで $`\mathcal V=\int_Y e^{123}`$ は物理体積である。Clifford 作用には $`c(e^j)=-i\sigma_j`$ 、 $`c_1c_2c_3=-1`$ を用いる。軸性ねじれによるシフトと固定電荷における諸量の対応は[付録A](DPPUv8-paper08_appA.md)に示す。実数 $`h`$ に対するこれらの式は、指定した斉次背景についてのものである。 $`\mathcal V`$ が与えられれば、電荷を固定するだけで $`w`$ が定まる。一方、再スケールする前のスカラー場の初速度を固定すると、 $`w`$ は $`h`$ に応じて変わる。

空間 Dirac 作用素の閉還元部分空間を $`\mathcal H`$ とし、次の二つの正の数を定義する。

<a id="eq-gaps"></a>

```math
\begin{aligned}
q_+&=4\mathcal V\min\bigl(\mathrm{spec}_{\mathcal H}D_{\mathrm{LC}}\cap(0,\infty)\bigr),
\\ q_-&=4\mathcal V\min\bigl(-\mathrm{spec}_{\mathcal H}D_{\mathrm{LC}}\cap(0,\infty)\bigr).
\end{aligned}
\qquad(2)
```

最初の接触、すなわち当該の準位が零に達するのは、 $`Q_\Sigma=q_+`$ および $`Q_\Sigma=-q_-`$ のときである。本稿の例ではどれもスペクトルが正負両方の固有値をもつ。 $`Q_\Sigma=0`$ ですでに存在する零モードは別に記録する。二つの閾値が異なり、一方の閾値に電荷の絶対値を合わせ、そのセクターを測定できるとき、この閾値によって符号を判別できる。絶対値が未知の任意の電荷について符号を同定できるとは限らない。

本稿で用いる η 不変量は、重複度を含めた和

```math
\eta_D(z)=\sum_{\lambda\ne0}\mathrm{sgn}(\lambda)\lvert\lambda\rvert^{-z}
```

を $`z=0`$ に解析接続した値である。和から零モードを除外し、核の次元の半分を加える補正も行わない。特記しない限り、η は切片全体におけるシフト前の作用素の値を指す。Sol の定理で測定する閾値は非零ファイバー運動量セクター、Nil の主要な閾値は非零中心運動量セクターに対するものである。

[英語原文](LaTeX/sections/sec02.tex)

[← 第1節](DPPUv8-paper08_sec01.md)　[目次](README_ja.md)　[第3節 →](DPPUv8-paper08_sec03.md)
