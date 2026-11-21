import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
import csv
# ---- 收单员：命令行点菜 ----
ap = argparse.ArgumentParser()
ap.add_argument("--Kp", type=float, default=4.69)
ap.add_argument("--Ki", type=float, default=12.84)
ap.add_argument("--Kd", type=float, default=0.429)
ap.add_argument("--out", default="sim")
args = ap.parse_args()

TAUS = [0.2, 0.2, 0.2]
DT, T = 0.001, 12.0

def run(Kp, Ki=0.0, Kd=0.0):
    x = [0.0, 0.0, 0.0]
    integ, e_prev = 0.0, None
    ts, ys, es, us = [], [], [], []
    for i in range(int(T / DT)):
        e = 1.0 - x[2]
        es.append(e)
        integ += e * DT
        deriv = 0.0 if e_prev is None else (e - e_prev) / DT
        e_prev = e
        u = Kp * e + Ki * integ + Kd * deriv
        us.append(u)
        inp = [u, x[0], x[1]]
        nx = [x[k] + DT * ((-x[k] + inp[k]) / TAUS[k]) for k in range(3)]
        x = nx
        ts.append(i * DT)
        ys.append(x[2])
    return np.array(ts), np.array(ys), np.array(es), np.array(us)

def amp(ys, a, b):
    seg = ys[a:b]
    return seg.max() - seg.min()

def ratio(ys):
    w = int(2 / DT)
    return amp(ys, len(ys) - w, len(ys)) / amp(ys, len(ys) - 2 * w, len(ys) - w)

def peaks(ts, ys):
    idx = [i for i in range(1, len(ys) - 1) if ys[i] > ys[i-1] and ys[i] > ys[i+1]]
    return [(ts[i], ys[i]) for i in idx]

# ---- 标定段：粗扫 + 二分找 Ku ----
print("Kp      末2s/前2s    判定")
for kp in range(4, 12):
    ts, ys, _, _ = run(kp)
    r = ratio(ys)
    tag = "衰减" if r < 0.95 else ("发散" if r > 1.05 else "等幅")
    print(f"{kp:<6}{r:<13.3f}{tag}")

lo, hi = 7.0, 8.0
for _ in range(40):
    mid = (lo + hi) / 2
    ts, ys, _, _ = run(mid)
    if ratio(ys) > 1:
        hi = mid
    else:
        lo = mid
Ku = (lo + hi) / 2

# ---- 量 Tu：只用纯 P 的 run(Ku)，这几行里不许出现 args ----
t_ku, y_ku, _, _ = run(Ku)
p = peaks(t_ku, y_ku)
Tu = p[-1][0] - p[-2][0]
print(f"\nKu = {Ku:.2f}   Tu = {Tu:.3f} s   最后两峰: {p[-2][0]:.3f}s / {p[-1][0]:.3f}s")

Kp_zn = 0.60 * Ku
Ti_zn = 0.50 * Tu
Td_zn = 0.125 * Tu
Ki_zn = Kp_zn / Ti_zn
Kd_zn = Kp_zn * Td_zn
print(f"Kp={Kp_zn:.2f}  Ti={Ti_zn:.3f}s  Td={Td_zn:.3f}s  Ki={Ki_zn:.2f}  Kd={Kd_zn:.3f}")

# ---- 使用段：参赛名单（改动②）----
groups = [
    ("P only", Kp_zn, 0.0,   0.0),
    ("PI",     Kp_zn, 6.0,   0.0),
    ("Z-N",    Kp_zn, Ki_zn, Kd_zn),
    ("menu",   args.Kp, args.Ki, args.Kd),
]

results = []
for name, kp, ki, kd in groups:
    ts, ys, es, us = run(kp, ki, kd)
    results.append((name, ts, ys, es, us))

# ---- 画图：四个机位（改动③）----
def label_of(name):
    if name == "menu":
        return f"menu Kp={args.Kp} Ki={args.Ki} Kd={args.Kd}"
    return name

# 机位 1：四条 y
plt.figure(figsize=(9, 5))
for name, ts, ys, es, us in results:
    plt.plot(ts, ys, label=label_of(name))
plt.axhline(1, ls="--", c="gray", lw=1)
plt.legend()
plt.tight_layout()
plt.savefig(args.out + "_1.png", dpi=150)
print("saved", args.out + "_1.png")

# 机位 2：四条 e
plt.figure(figsize=(9, 5))
for name, ts, ys, es, us in results:
    plt.plot(ts, es, label=label_of(name))
plt.axhline(0, ls="--", c="gray", lw=1)
plt.legend()
plt.tight_layout()
plt.savefig(args.out + "_2.png", dpi=150)
print("saved", args.out + "_2.png")

# 机位 3：四条 u
plt.figure(figsize=(9, 5))
for name, ts, ys, es, us in results:
    plt.plot(ts, us, label=label_of(name))
plt.axhline(0, ls="--", c="gray", lw=1)
plt.legend()
plt.tight_layout()
plt.savefig(args.out + "_3.png", dpi=150)
print("saved", args.out + "_3.png")

# 机位 4：y 的最后 4 秒放大（看 ess）
plt.figure(figsize=(9, 5))
for name, ts, ys, es, us in results:
    m = ts > T - 4
    plt.plot(ts[m], ys[m], label=label_of(name))
plt.axhline(1, ls="--", c="gray", lw=1)
plt.legend()
plt.tight_layout()
plt.savefig(args.out + "_4.png", dpi=150)
print("saved", args.out + "_4.png")


rows = []
for name, ts, ys, es, us in results:
    peak = ys.max()                                   # 峰值
    ts_last = 0.0                                     # 最后一次离开 ±2% 带的时刻
    for i in range(len(ys)):
        if abs(ys[i] - 1) > 0.02:
            ts_last = ts[i]
    ess = es[ts > T - 0.5].mean()                     # 最后 0.5 秒的平均误差 = ess
    rows.append([name, f"{peak:.3f}", f"{(peak-1)*100:.1f}", f"{ts_last:.2f}", f"{ess:.3f}"])

with open(args.out + ".csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "peak", "overshoot_pct", "ts", "ess"])   # 表头
    w.writerows(rows)
print("saved", args.out + ".csv")