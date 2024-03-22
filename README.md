# Derivation_of_Kubo_Bastin
# 目的
Kubo-Bastin公式の導出を温度グリーン関数から行う。
# 内容
## Kubo-Bastin公式の導出
(静的)電気伝導度は電流自己相関関数
```math
\Phi_{j^{\mu}_e j^{\nu}_e}(\tau) = \frac{1}{V}\braket{T_\tau[j^{\mu}_e(\tau) j^{\nu}_e(0)]}
```
の松原表示を実周波数($`\hbar \omega + i\delta`$)に解析接続したのち、その周波数に関して差分をとることで
```math
\sigma_{\mu \nu} = \lim_{\omega\rightarrow 0} \frac{\Phi^{R}_{j^{\mu}_e j^{\nu}_e}(\omega) - \Phi^{R}_{j^{\mu}_e j^{\nu}_e}(0)}{i\omega}  
```
と得られる。(ちなみに光学伝導度は$`\omega`$を有限にすることで得られる。)ここで多体効果(電子間相互作用、電子-格子相互作用)から生じるバーテックス補正を無視する近似(独立電子近似)を行うと電流相関関数は
```math
\Phi_{j^{\mu}_e j^{\nu}_e}(\tau) = \frac{1}{V}\braket{T_\tau[j^{\mu}_e(\tau) j^{\nu}_e(0)]} = -\frac{e^2}{V} \sum_{\vec{k}}{\mathrm{tr}\vec{v}_\vec{k} G(\vec{k},\tau) \vec{v}_\vec{k} G(\vec{k},-\tau)} 
```
となる。$`G(\tau),G(-\tau)`$の２つの符号が出るのは電流演算子が生成、消滅演算子のペアからなり、消滅演算子が$`\tau`$生成演算子が$`0`$のペアを端点として持つグリーン関数が$`G(\tau)`$, 消滅演算子が$`0`$生成演算子が$`\tau`$のペアを端点としてもつグリーン関数が$`G(-\tau)`$となるため。また、負符号はグリーン関数のループから生じる。虚時間の式を松原表示すると
```math
\Phi_{j^{\mu}_e j^{\nu}_e}(i\omega_\nu) = -\frac{e^2}{\beta V} \sum_{\vec{k},n}{\mathrm{tr}v^{\mu}_\vec{k} G(\vec{k},i\omega_n ) v^{\nu}_\vec{k} G(\vec{k}, i\omega_n - i\omega_\nu)} 
```
次に、松原和を留数定理を使い、フェルミディラック積分に置き換えることを考える。フェルミディラック分布関数を$`\mu=0`$とした式は
```math
f(\varepsilon) = \frac{1}{e^{\beta(\varepsilon)} + 1}
```
であるが、これは$`i\omega_n = n \pi k_\mathrm{B} T (n:\mathrm{integer})`$に一位の極をもつ。これに松原周波数で正則な関数$`g(z)`$をかけた関数の松原周波数における留数は
$`-\frac{1}{\beta} g(i \omega_n)`$である。これより
```math
\Phi_{j^{\mu}_e j^{\nu}_e}(i\omega_\nu) = \frac{e^2}{2\pi i V} \int_\mathrm{C}{dz f_\mathrm{FD}(z)\sum_{\vec{k},n}{\mathrm{tr}v^{\mu}_\vec{k} G(\vec{k}, z) v^{\nu}_\vec{k} G(\vec{k}, z - i\omega_\nu)}} 
```
となる。ただし、エネルギーをシフトし、$f_{\mathrm{FD}}(\varepsilon) = \frac{1}{e^{\beta(\varepsilon - \mu} + 1}}$となるようにした。 
積分路Cは分枝線をよけ、外側では円状に取る。(詳しく知りたい場合はAGDなどをみる。)外側の円の部分の寄与は無限遠で0に近づくので残るのは分枝線の上側の実軸方向に負の無限大から
正の無限大への積分と、分枝線の下側の逆向きの積分の寄与だけとなる。
分枝線は$`\mathrm{Im}z = \omega_\nu,0`$に存在するので積分は
```math
\int_{-\infty}^{\infty}{d\varepsilon f_\mathrm{FD}(\varepsilon) \frac{e^2}{2\pi i V}\sum_{\vec{k}}{\mathrm{tr}\left[v^{\mu}_\vec{k} G(\vec{k}, \varepsilon + i\omega_\nu) v^{\nu}_\vec{k} (G(\vec{k},\varepsilon+ i\delta) - G(\vec{k},\varepsilon- i\delta )) + v^{\mu}_\vec{k} (G(\vec{k},\varepsilon+ i\delta) - G(\vec{k},\varepsilon- i\delta ))v^{\nu}_\vec{k}G(\vec{k}, \varepsilon - i\omega_\nu)\right]}}
```
となる。静的伝導度を計算するときは周波数を$`\hbar \omega + i\delta`$に解析接続した後に上の式の$`\omega`$の１次の項をとればよく
```math
\sigma_{\mu \nu} = \int_{-\infty}^{\infty}{d\varepsilon f_\mathrm{FD}(\varepsilon) \frac{\hbar e^2}{i \pi V}\sum_{\vec{k}}{\mathrm{tr}\left[v^{\mu}_\vec{k} \frac{\partial}{\partial \varepsilon}G(\vec{k}, \varepsilon + i\delta) v^{\nu}_\vec{k} \mathrm{Im} G(\vec{k},\varepsilon+ i\delta)  - v^{\mu}_\vec{k} \mathrm{Im} G(\vec{k},\varepsilon+ i\delta) v^{\nu}_\vec{k} \frac{\partial}{\partial \varepsilon} G(\vec{k}, \varepsilon - i\delta)\right]}}
```
これがKubo-Bastin公式である。
## Kubo-Streda公式の導出

## Kubo-Greenwood公式の導出
$`\mu = \nu`$つまり縦の電気伝導度の場合はトレースの巡回則を使い、フェルミディラック以外のすべての部分全体を微分した形に書ける。これを部分積分することで
```math
\sigma_{\mu \mu} = \int_{-\infty}^{\infty}{d\varepsilon \left(-\frac{df_\mathrm{FD}(\varepsilon)}{d\varepsilon}\right) \frac{\hbar e^2}{\pi V}\sum_{\vec{k}}{\mathrm{tr}\left[v^{\mu}_\vec{k} \mathrm{Im} G(\vec{k}, \varepsilon + i\delta) v^{\mu}_\vec{k} \mathrm{Im} G(\vec{k},\varepsilon+ i\delta) \right]}}
```
となる。これを$`T= 0\ \mathrm{K}`$で評価するとフェルミディラック分布の微分がデルタ関数になり、
```math
\sigma_{\mu \mu} = \frac{\hbar e^2}{\pi V}\sum_{\vec{k}}{\mathrm{tr}\left[v^{\mu}_\vec{k} \mathrm{Im} G(\vec{k}, \varepsilon + i\delta) v^{\mu}_\vec{k} \mathrm{Im} G(\vec{k},\varepsilon+ i\delta)\right]}
```

