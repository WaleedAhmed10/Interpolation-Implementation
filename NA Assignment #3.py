"""
Lab 10 - Numerical Differentiation using Interpolation Formulae
Numerical Analysis

Dataset : Temperature of a cup of hot water cooling in a room (Temperature vs Time)
Methods : Forward, Backward and Central differences (first and second derivative)

Run the GUI :  python lab10_numerical_differentiation.py
Run CLI only:  python lab10_numerical_differentiation.py --cli
Requires    :  numpy, matplotlib (tkinter ships with Python)
"""
import sys
import numpy as np


TIME = np.arange(0, 30, 2, dtype=float)
TEMP = np.array([95.00, 88.34, 82.31, 76.86, 71.92, 67.46, 63.42, 59.76,
                 56.45, 53.46, 50.75, 48.30, 46.08, 44.08, 42.26])
T_LABEL, Y_LABEL = "Time (min)", "Temperature (°C)"
D1_UNIT, D2_UNIT = "°C/min", "°C/min²"


def ref_f(t):  return 25 + 70 * np.exp(-0.05 * t)
def ref_d1(t): return -3.5 * np.exp(-0.05 * t)
def ref_d2(t): return 0.175 * np.exp(-0.05 * t)


def first_derivative(y, h, method):
    n = len(y)
    d = np.full(n, np.nan)
    if method == "forward":
        d[:-1] = (y[1:] - y[:-1]) / h
    elif method == "backward":
        d[1:] = (y[1:] - y[:-1]) / h
    elif method == "central":
        d[1:-1] = (y[2:] - y[:-2]) / (2 * h)
    else:
        raise ValueError(method)
    return d


def second_derivative(y, h, method):
    n = len(y)
    d = np.full(n, np.nan)
    if method == "forward":
        d[:-2] = (y[:-2] - 2 * y[1:-1] + y[2:]) / h**2
    elif method == "backward":
        d[2:] = (y[2:] - 2 * y[1:-1] + y[:-2]) / h**2
    elif method == "central":
        d[1:-1] = (y[2:] - 2 * y[1:-1] + y[:-2]) / h**2
    else:
        raise ValueError(method)
    return d


FIRST_METHODS = ["forward", "backward", "central"]
SECOND_METHODS = ["forward", "backward", "central"]
ORDER = {"forward": "O(h)", "backward": "O(h)", "central": "O(h²)"}

JUSTIFICATION = """WHY EACH METHOD IS CHOSEN
-------------------------
Forward difference  (y[i+1]-y[i])/h
   Uses only FUTURE points, so it is the only choice at the FIRST reading
   (t = 0) where no earlier data exists. Accuracy O(h): least accurate.

Backward difference (y[i]-y[i-1])/h
   Uses only PAST points, so it is the only choice at the LAST reading
   (t = 28) where no later data exists. Also O(h). It is natural for
   real-time monitoring because it needs no future data.

Central difference  (y[i+1]-y[i-1])/(2h)
   Uses points on BOTH sides, so it is used at all INTERIOR points.
   Errors of the two sides cancel -> O(h^2), the most accurate of the three.

Second derivative (y[i+1]-2y[i]+y[i-1])/h^2
   Central form is O(h^2) and used at interior points. Forward/backward
   forms need two extra points on one side and are only O(h).
   Second derivatives amplify rounding error of the data (error ~ eps/h^2),
   so they are noisier than first derivatives.
"""


def compute_all(t, y):
    h = t[1] - t[0]
    if not np.allclose(np.diff(t), h):
        raise ValueError("Time values must be equally spaced.")
    first = {m: first_derivative(y, h, m) for m in FIRST_METHODS}
    second = {m: second_derivative(y, h, m) for m in SECOND_METHODS}
    return h, first, second


def fmt(v, w=9, p=4):
    return " " * (w - 1) + "-" if np.isnan(v) else f"{v:{w}.{p}f}"


def table_text(t, y, first, second, show=("forward", "backward", "central")):
    L = []
    L.append(f"FIRST DERIVATIVE  dT/dt  [{D1_UNIT}]")
    head = f"{'t(min)':>7}{'T(°C)':>9}" + "".join(f"{m:>10}" for m in show) + f"{'exact':>10}"
    L += [head, "-" * len(head)]
    for i in range(len(t)):
        L.append(f"{t[i]:7.1f}{y[i]:9.2f}" + "".join(fmt(first[m][i], 10) for m in show)
                 + fmt(ref_d1(t[i]), 10))
    L.append("")
    L.append(f"SECOND DERIVATIVE  d²T/dt²  [{D2_UNIT}]")
    head = f"{'t(min)':>7}{'T(°C)':>9}" + "".join(f"{m:>10}" for m in SECOND_METHODS) + f"{'exact':>10}"
    L += [head, "-" * len(head)]
    for i in range(len(t)):
        L.append(f"{t[i]:7.1f}{y[i]:9.2f}" + "".join(fmt(second[m][i], 10, 5) for m in SECOND_METHODS)
                 + fmt(ref_d2(t[i]), 10, 5))
    return "\n".join(L)


def error_text(t, first, second):
    L = ["ACCURACY (mean absolute error vs exact model, over points where defined)", ""]
    L.append(f"{'First derivative':<22}{'Order':>8}{'MAE ['+D1_UNIT+']':>16}{'Max error':>12}")
    for m in FIRST_METHODS:
        e = np.abs(first[m] - ref_d1(t)); e = e[~np.isnan(e)]
        L.append(f"{m:<22}{ORDER[m]:>8}{e.mean():16.5f}{e.max():12.5f}")
    L.append("")
    L.append(f"{'Second derivative':<22}{'Order':>8}{'MAE ['+D2_UNIT+']':>16}{'Max error':>12}")
    for m in SECOND_METHODS:
        e = np.abs(second[m] - ref_d2(t)); e = e[~np.isnan(e)]
        L.append(f"{m:<22}{ORDER[m]:>8}{e.mean():16.5f}{e.max():12.5f}")
    return "\n".join(L)


def conclusion_text(t, first, second):
    mae1 = {m: np.nanmean(np.abs(first[m] - ref_d1(t))) for m in FIRST_METHODS}
    mae2 = {m: np.nanmean(np.abs(second[m] - ref_d2(t))) for m in SECOND_METHODS}
    b1, b2 = min(mae1, key=mae1.get), min(mae2, key=mae2.get)
    return (f"CONCLUSION\n"
            f"- Temperature falls from {TEMP[0]:.0f} °C to {TEMP[-1]:.0f} °C; the first derivative is negative "
            f"(cooling) and its magnitude shrinks from about {abs(ref_d1(t[0])):.2f} to {abs(ref_d1(t[-1])):.2f} °C/min, "
            f"while the second derivative is small and positive (the cooling rate is slowing down).\n"
            f"- Most accurate first-derivative method here: {b1} (MAE {mae1[b1]:.5f}).\n"
            f"- Most accurate second-derivative method here: {b2} (MAE {mae2[b2]:.5f}).\n"
            f"- Central difference beats forward/backward because its error is O(h²) versus O(h); "
            f"forward/backward differences systematically over/under-estimate on a curved (exponential) signal.\n"
            f"- Second derivatives are noticeably noisier because the 2-decimal rounding of the data is divided by h².")


STYLE = {"forward": ("tab:red", "s--"), "backward": ("tab:green", "^--"),
         "central": ("tab:blue", "o-")}


def draw_figure(fig, t, y, first, second, m1, m2, show_exact=True):
    fig.clear()
    ax = fig.subplots(3, 1, sharex=True)
    ax[0].plot(t, y, "ko-", ms=4, label="Measured data")
    ax[0].set_ylabel(Y_LABEL); ax[0].set_title("Temperature vs Time"); ax[0].grid(alpha=.3); ax[0].legend()
    tt = np.linspace(t[0], t[-1], 300)
    if show_exact:
        ax[1].plot(tt, ref_d1(tt), "k-", lw=1, alpha=.6, label="Exact (model)")
        ax[2].plot(tt, ref_d2(tt), "k-", lw=1, alpha=.6, label="Exact (model)")
    for m in m1:
        c, s = STYLE[m]; ax[1].plot(t, first[m], s, color=c, ms=4, label=m)
    for m in m2:
        c, s = STYLE[m]; ax[2].plot(t, second[m], s, color=c, ms=4, label=m)
    ax[1].set_ylabel(f"dT/dt ({D1_UNIT})"); ax[1].set_title("First derivative"); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8)
    ax[2].set_ylabel(f"d²T/dt² ({D2_UNIT})"); ax[2].set_title("Second derivative"); ax[2].grid(alpha=.3); ax[2].legend(fontsize=8)
    ax[2].set_xlabel(T_LABEL)
    fig.tight_layout()


def draw_error_figure(fig, t, first, second):
    fig.clear()
    a, b = fig.subplots(1, 2)
    for m in FIRST_METHODS:
        c, s = STYLE[m]; a.semilogy(t, np.abs(first[m] - ref_d1(t)) + 1e-12, s, color=c, ms=4, label=m)
    for m in SECOND_METHODS:
        c, s = STYLE[m]; b.semilogy(t, np.abs(second[m] - ref_d2(t)) + 1e-12, s, color=c, ms=4, label=m)
    a.set_title("|error| first derivative"); b.set_title("|error| second derivative")
    for x in (a, b): x.set_xlabel(T_LABEL); x.grid(alpha=.3, which="both"); x.legend(fontsize=8)
    a.set_ylabel(f"|error| ({D1_UNIT})"); b.set_ylabel(f"|error| ({D2_UNIT})")
    fig.tight_layout()


def run_cli(save=True):
    h, first, second = compute_all(TIME, TEMP)
    print(f"Data points: {len(TIME)}   step h = {h:g} min\n")
    print(JUSTIFICATION)
    print(table_text(TIME, TEMP, first, second, show=FIRST_METHODS)); print()
    print(error_text(TIME, first, second)); print()
    print(conclusion_text(TIME, first, second))
    if save:
        import matplotlib; matplotlib.use("Agg")
        from matplotlib.figure import Figure
        f = Figure(figsize=(9, 10)); draw_figure(f, TIME, TEMP, first, second, FIRST_METHODS[:3], SECOND_METHODS)
        f.savefig("lab10_derivatives.png", dpi=150)
        g = Figure(figsize=(11, 4.5)); draw_error_figure(g, TIME, first, second)
        g.savefig("lab10_errors.png", dpi=150)
        print("\nSaved lab10_derivatives.png and lab10_errors.png")


def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    root = tk.Tk()
    root.title("Lab 10 - Numerical Differentiation (Temperature vs Time)")
    root.geometry("1250x780")
    state = {"t": TIME.copy(), "y": TEMP.copy()}


    left = ttk.Frame(root, padding=8); left.pack(side="left", fill="y")
    ttk.Label(left, text="Data (editable)", font=("", 11, "bold")).pack(anchor="w")
    data_tree = ttk.Treeview(left, columns=("t", "y"), show="headings", height=16)
    data_tree.heading("t", text="Time (min)"); data_tree.heading("y", text="Temp (°C)")
    data_tree.column("t", width=90, anchor="center"); data_tree.column("y", width=90, anchor="center")
    data_tree.pack()

    def load_data_tree():
        data_tree.delete(*data_tree.get_children())
        for a, b in zip(state["t"], state["y"]):
            data_tree.insert("", "end", values=(f"{a:g}", f"{b:g}"))
    load_data_tree()

    def edit_cell(ev):
        item = data_tree.identify_row(ev.y); col = data_tree.identify_column(ev.x)
        if not item: return
        ci = int(col[1:]) - 1
        x, y, w, h_ = data_tree.bbox(item, col)
        e = ttk.Entry(data_tree); e.place(x=x, y=y, width=w, height=h_)
        e.insert(0, data_tree.item(item, "values")[ci]); e.focus()
        def done(_=None):
            vals = list(data_tree.item(item, "values")); vals[ci] = e.get()
            data_tree.item(item, values=vals); e.destroy()
        e.bind("<Return>", done); e.bind("<FocusOut>", done)
    data_tree.bind("<Double-1>", edit_cell)

    def read_tree():
        try:
            rows = [tuple(float(v) for v in data_tree.item(i, "values")) for i in data_tree.get_children()]
            t = np.array([r[0] for r in rows]); y = np.array([r[1] for r in rows])
            if len(t) < 5: raise ValueError("Need at least 5 points.")
            return t, y
        except Exception as ex:
            messagebox.showerror("Data error", str(ex)); return None

    def add_row():
        data_tree.insert("", "end", values=("0", "0"))
    def del_row():
        for i in data_tree.selection(): data_tree.delete(i)
    def load_csv():
        p = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if not p: return
        try:
            arr = np.genfromtxt(p, delimiter=",", skip_header=0)
            if np.isnan(arr[0]).any(): arr = np.genfromtxt(p, delimiter=",", skip_header=1)
            state["t"], state["y"] = arr[:, 0], arr[:, 1]; load_data_tree(); update()
        except Exception as ex:
            messagebox.showerror("CSV error", f"Expected two columns (time, value).\n{ex}")
    def reset():
        state["t"], state["y"] = TIME.copy(), TEMP.copy(); load_data_tree(); update()

    bf = ttk.Frame(left); bf.pack(fill="x", pady=4)
    for txt, cmd in [("Add row", add_row), ("Delete row", del_row), ("Load CSV", load_csv), ("Reset", reset)]:
        ttk.Button(bf, text=txt, command=cmd).pack(side="left", padx=1)

    ttk.Separator(left).pack(fill="x", pady=6)
    ttk.Label(left, text="First-derivative methods", font=("", 10, "bold")).pack(anchor="w")
    v1 = {m: tk.BooleanVar(value=m in ("forward", "backward", "central")) for m in FIRST_METHODS}
    for m in FIRST_METHODS:
        ttk.Checkbutton(left, text=f"{m}  {ORDER[m]}", variable=v1[m], command=lambda: update()).pack(anchor="w")
    ttk.Label(left, text="Second-derivative methods", font=("", 10, "bold")).pack(anchor="w", pady=(8, 0))
    v2 = {m: tk.BooleanVar(value=True) for m in SECOND_METHODS}
    for m in SECOND_METHODS:
        ttk.Checkbutton(left, text=f"{m}  {ORDER[m]}", variable=v2[m], command=lambda: update()).pack(anchor="w")
    exact_v = tk.BooleanVar(value=True)
    ttk.Checkbutton(left, text="Show exact model curve", variable=exact_v, command=lambda: update()).pack(anchor="w", pady=(8, 0))
    ttk.Button(left, text="Compute / Refresh", command=lambda: update()).pack(fill="x", pady=8)
    ttk.Button(left, text="Save plot (PNG)", command=lambda: save_plot()).pack(fill="x")
    ttk.Button(left, text="Export table (CSV)", command=lambda: export_csv()).pack(fill="x", pady=4)


    nb = ttk.Notebook(root); nb.pack(side="right", fill="both", expand=True, padx=6, pady=6)
    tab_plot, tab_tab, tab_err, tab_why = (ttk.Frame(nb) for _ in range(4))
    nb.add(tab_plot, text="Graphs"); nb.add(tab_tab, text="Tables")
    nb.add(tab_err, text="Accuracy & Conclusion"); nb.add(tab_why, text="Method justification")

    fig = Figure(figsize=(8, 7)); canvas = FigureCanvasTkAgg(fig, master=tab_plot)
    NavigationToolbar2Tk(canvas, tab_plot).update()
    canvas.get_tk_widget().pack(fill="both", expand=True)

    def mk_text(parent):
        tx = tk.Text(parent, font=("Courier New", 10), wrap="none"); tx.pack(fill="both", expand=True); return tx
    txt_tab, txt_err, txt_why = mk_text(tab_tab), mk_text(tab_err), mk_text(tab_why)
    txt_why.insert("1.0", JUSTIFICATION); txt_why.config(state="disabled")

    cache = {}
    def update():
        d = read_tree()
        if d is None: return
        t, y = d
        try:
            h, first, second = compute_all(t, y)
        except ValueError as ex:
            messagebox.showerror("Data error", str(ex)); return
        state["t"], state["y"] = t, y
        cache.update(t=t, y=y, first=first, second=second)
        m1 = [m for m in FIRST_METHODS if v1[m].get()]
        m2 = [m for m in SECOND_METHODS if v2[m].get()]
        draw_figure(fig, t, y, first, second, m1, m2, exact_v.get()); canvas.draw()
        for tx, s in ((txt_tab, f"h = {h:g} min,  n = {len(t)} points\n\n" + table_text(t, y, first, second, FIRST_METHODS)),
                      (txt_err, error_text(t, first, second) + "\n\n" + conclusion_text(t, first, second))):
            tx.config(state="normal"); tx.delete("1.0", "end"); tx.insert("1.0", s); tx.config(state="disabled")

    def save_plot():
        p = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png")])
        if p: fig.savefig(p, dpi=150)
    def export_csv():
        if not cache: return
        p = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not p: return
        cols = [cache["t"], cache["y"]] + [cache["first"][m] for m in FIRST_METHODS] + [cache["second"][m] for m in SECOND_METHODS]
        hdr = "t,T," + ",".join(f"d1_{m}" for m in FIRST_METHODS) + "," + ",".join(f"d2_{m}" for m in SECOND_METHODS)
        np.savetxt(p, np.column_stack(cols), delimiter=",", header=hdr, comments="", fmt="%.6f")

    update()
    root.mainloop()


if __name__ == "__main__":
    if "--cli" in sys.argv:
        run_cli()
    else:
        try:
            run_gui()
        except Exception as ex:
            print("GUI could not start:", ex, "\nFalling back to CLI mode.\n")
            run_cli()
