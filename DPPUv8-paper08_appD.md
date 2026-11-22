# 付録D　Nil のスペクトルと η 不変量の補正

水平格子の周期を $`P=2\pi`$ 、中心方向の周期を $`P^2/r`$ とする。体積は $`V_0W=P^4a^2c/r`$ である。物理的なスカラーシフトを一般に $`\mu`$ とすると、スペクトルの中心は $`b_0=\mu-t/4,\ t=c/a^2`$ である。中心方向の運動量が零のセクターは $`\delta_3=0`$ の場合に限り存在する。このセクターでは $`b_0`$ から測った準位が $`\pm\lvert\boldsymbol l\rvert/a,\ \boldsymbol l\in\mathbb Z^2+(\delta_1,\delta_2)`$ となり、 $`\boldsymbol l=0`$ のときも重複度は2である。非零の中心周波数 $`k=r\lvert j+\delta_3\rvert,\ \nu=k/(2\pi)`$ に対して、 $`b_0`$ から測った残りの準位は

$$
f=-\nu/c,\qquad
f_n^\pm=\pm\sqrt{(\nu/c)^2+2\nu n/a^2},\quad n\ge1
$$

となる。各準位の多重度は、中心運動量の両符号を含めて $`2k`$ である。実際、中心方向に Fourier 展開すると、水平方向の Clifford 微分は振動子の上昇・下降作用素の対となる。隣接する Hermite 次数が張る二次元ブロックの対角成分は $`\pm\nu/c`$ 、非対角成分は $`\sqrt{2\nu n}/a`$ である。対をもたずに残るベクトルの固有値は $`-\nu/c`$ となる。水平一周期の並進は第2の Fourier モード番号を $`k`$ だけずらすので、中心運動量の符号ごとに独立な剰余類が $`k`$ 個残る。対をなす塔の正の準位と、対をもたない塔の負の準位でそれぞれ原点に最も近いものを取ると、ギャップに関する仮定が成り立つ範囲で式(7) が従う。両閾値が一致する点は $`R_*=\nu_*/t+t/2`$ の両辺を二乗して得られる。[第5節・式(7)](DPPUv8-paper08_sec05.md#eq-nilmin)

中心運動量が非零のセクターでスピン接続の面積項 $`-t/4`$ を与えるのは Koszul の接続公式である。η 不変量を Gornet–Richardson の定理12と比較するには、同文献の幾何パラメータを $`A_{\rm GR}=1/t,\ r_{\rm GR}=P^2c/r,\ m_v=r`$ と置く。シフトを $`\mu=t/4`$ としてスペクトルの中心を原点に合わせると、対をなす塔の寄与が相殺する。対をもたない塔の寄与 $`-2r\zeta(-1,\alpha)=C_r`$ は、 $`\alpha=1`$ （周期）のとき $`r/6`$ 、 $`\alpha=1/2`$ （反周期）のとき $`-r/12`$ である。局所熱核係数は滑らかな部分の微分 $`\partial_\mu\eta_{\rm sm}=\mathcal V(R_3/24-\mu^2)/\pi^2,\ R_3=-t^2/2`$ を与える。[文献5・定理12](REFERENCES.md#ref5)

固有値が零を横切る寄与を向きに応じた符号で加えると、

$$
\begin{aligned}
\eta(\mu)&=C_r+\frac{\mathcal V}{\pi^2}
\left(\frac{(t/4)^3-\mu^3}{3}+\frac{(\mu-t/4)R_3}{24}\right)+J,\\
J&=\sum_f d_f\bigl\lbrace \mathrm{sgn}(f+\mu-t/4)-\mathrm{sgn}\thinspace f\bigr\rbrace ,
\qquad \mathrm{sgn}(0)=0 .
\end{aligned}
$$

この有限和に寄与するのは $`\lvert f\rvert\le\lvert\mu-t/4\rvert`$ の準位だけである。 $`\mu=0`$ では滑らかな項は $`\mathcal V t^3/(96\pi^2)`$ となる。 $`J`$ に水平モードが寄与するのは中心方向が周期的なスピン構造の場合に限る。 $`a=c=1,\ r=2,\ (\delta_1,\delta_2)=(0,1/2)`$ では、通常の η 不変量は中心方向が周期のとき約 $`1.15580037`$ 、反周期のとき約 $`0.65580037`$ である。 $`C_r`$ の符号が変わっても、両方の値は正である。

[英語原文](LaTeX/appendices/appD.tex)・[閾値表](data/nil_thresholds.md)

[← 付録C](DPPUv8-paper08_appC.md)　[目次](README_ja.md)　[付録E →](DPPUv8-paper08_appE.md)
