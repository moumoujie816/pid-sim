# PID 仿真

三阶惯性对象（τ=0.2×3）+ PID，argparse 可配置，一键四图 + CSV。
含 Ziegler-Nichols 临界比例度法自动整定（粗扫 + 二分找 Ku）。

## 复现
pip install numpy matplotlib
python3 sim_all.py                          # 默认参数，出 4 图 + CSV
python3 sim_all.py --Kp 3.8 --Ki 3.0 --Kd 0.5 --out hand   # 换参数对比

## 输出
sim_1.png ~ sim_4.png（y / e / u / 稳态放大）、sim.csv