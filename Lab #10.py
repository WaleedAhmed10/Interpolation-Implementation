import math
import sys

H = 0.1
N = 11
KMAX = 4

COEF_LABEL = {
    1: "1",
    2: "(2p-1)/2",
    3: "(3p^2-6p+2)/6",
    4: "(4p^3-18p^2+22p-6)/24",
}
SUP = ["", "Δ", "Δ²", "Δ³", "Δ⁴", "Δ⁵", "Δ⁶", "Δ⁷", "Δ⁸", "Δ⁹", "Δ¹⁰"]


def f(x):
    return x ** 3 + 2 * x ** 2 + x + 1


def df_exact(x):
    return 3 * x ** 2 + 4 * x + 1


def grid():
    return [round(i * H, 10) for i in range(N)]


def clean(v):
    v = round(v, 10)
    return 0.0 if abs(v) < 1e-9 else v


def build_table(xs):
    t = [[round(f(x), 4) for x in xs]]
    for _ in range(1, len(xs)):
        p = t[-1]
        t.append([clean(p[i + 1] - p[i]) for i in range(len(p) - 1)])
    return t


def constant_order(t):
    for k in range(1, len(t)):
        row = t[k]
        if len(row) >= 2 and max(row) - min(row) < 1e-9:
            return k
    return None


def binom(p, k):
    v = 1.0
    for j in range(k):
        v *= (p - j)
    return v / math.factorial(k)


def dbinom(p, k):
    if k == 0:
        return 0.0
    s = 0.0
    for m in range(k):
        v = 1.0
        for j in range(k):
            if j != m:
                v *= (p - j)
        s += v
    return s / math.factorial(k)


def evaluate(t, xs, idx, xe):
    x0 = xs[idx]
    p = round((xe - x0) / H, 10)
    exact = df_exact(xe)
    rows, cum = [], 0.0
    for k in range(1, KMAX + 1):
        d = t[k][idx]
        c = dbinom(p, k)
        term = c * d / H
        cum += term
        err = abs(cum - exact) / abs(exact) * 100 if exact else float("nan")
        rows.append(dict(k=k, label=COEF_LABEL[k], coef=c, delta=d,
                         term=term, cum=cum, err=err))
    interp = sum(binom(p, k) * t[k][idx] for k in range(KMAX + 1))
    return dict(x=xe, x0=x0, p=p, exact=exact, rows=rows,
                interp=interp, fexact=f(xe))


def fmt(v):
    v = 0.0 if abs(v) < 5e-13 else v
    return f"{v:.4f}"


def verdict(A, B, cord):
    ea, eb = A["rows"][cord - 1]["err"], B["rows"][cord - 1]["err"]
    e1a, e1b = A["rows"][0]["err"], B["rows"][0]["err"]
    same = abs(ea - eb) < 5e-5
    s = (f"Keeping all terms up to Δ^{cord}: error at x={A['x']:.4g} is {ea:.4f}% "
         f"and at x={B['x']:.4g} is {eb:.4f}% -> ")
    s += ("same accuracy (both exact to 4 d.p. because f is a cubic and Δ^4 = 0). "
          if same else "the two positions do NOT have the same accuracy. ")
    s += (f"They differ strongly with early truncation: keeping only Δ gives "
          f"{e1a:.4f}% at p={A['p']:.4g} versus {e1b:.4f}% at p={B['p']:.4g}; "
          f"the point further from the origin (|p| larger) is the extrapolated, "
          f"less reliable one.")
    return s


def paragraph(A, B, cord):
    far = B if abs(B["p"]) > abs(A["p"]) else A
    return (f"Origin used: x0 = {A['x0']:.1f} (so p = (x - x0)/h with h = {H}); "
            f"f'({A['x']:.4g}) sits exactly on the origin (p = {A['p']:.4g}) and "
            f"f'({B['x']:.4g}) is reached by stepping back from it (p = {B['p']:.4g}). "
            f"Terms kept: {cord} (Δ to Δ^{cord}), because Δ^{cord} is constant and "
            f"Δ^{cord + 1} and higher are zero, so the series is exact for this cubic "
            f"and truncation error is zero. The estimate I trust least is "
            f"f'({far['x']:.4g}): it has the largest |p| = {abs(far['p']):.4g}, i.e. "
            f"it is extrapolated away from the origin, and with fewer terms its "
            f"error grows fastest ({far['rows'][0]['err']:.4f}% with Δ alone).")


def build_report(idx=1, xa=0.1, xb=0.0):
    xs = grid()
    t = build_table(xs)
    cord = constant_order(t)
    A, B = evaluate(t, xs, idx, xa), evaluate(t, xs, idx, xb)
    L = []
    add = L.append
    bar = "=" * 100
    add(bar)
    add("LAB 09 - NUMERICAL DIFFERENTIATION USING NEWTON'S FORWARD INTERPOLATION")
    add("f(x) = x^3 + 2x^2 + x + 1,   0 <= x <= 1,   h = 0.1")
    add(bar)

    add("\nTASK 1 - FUNCTION VALUES (4 d.p.)")
    add("x    : " + " ".join(f"{x:7.1f}" for x in xs))
    add("f(x) : " + " ".join(f"{f(x):7.4f}" for x in xs))
    add(f"The grid is equally spaced because x_i = x_0 + i*h with the same step "
        f"h = {H} for i = 0..10 (1.0/10). This uniform spacing is what makes the "
        f"forward differences and the single variable p = (x - x0)/h valid.")

    add("\nTASK 2 - FORWARD DIFFERENCE TABLE")
    hdr = ["x", "f(x)"] + [f"Delta^{k}" for k in range(1, N)]
    add(" ".join(f"{h:>9}" for h in hdr))
    for i in range(N):
        row = [f"{xs[i]:9.4f}"] + [f"{fmt(t[k][i]):>9}" for k in range(N - i)]
        add(" ".join(row))
    add(f"\nDelta^{cord} is constant = {fmt(t[cord][0])} (= 3! * a * h^3 = 6*1*0.001); "
        f"Delta^{cord + 1} and higher are 0  ->  degree of polynomial = {cord}.")

    def deriv_block(R, task):
        add(f"\n{task}")
        add(f"Origin x0 = {R['x0']:.1f},  x = {R['x']:.4g},  p = (x - x0)/h = {R['p']:.4g}")
        add("f'(x0 + p h) = (1/h) [ D f0 + (2p-1)/2 D^2 f0 + (3p^2-6p+2)/6 D^3 f0 "
            "+ (4p^3-18p^2+22p-6)/24 D^4 f0 + ... ]")
        add("Term-by-term:")
        for r in R["rows"]:
            add(f"  Term {r['k']}: (1/{H}) * [{r['label']} = {r['coef']:+.4f}] * "
                f"({fmt(r['delta'])}) = {r['term']:+.4f}")
        add(f"Truncated after order {cord}: Delta^{cord} is constant, so Delta^{cord + 1}.. = 0 "
            f"and the remaining terms vanish (zero truncation error).")
        add(f"{'Truncation':<14}{'Numerical':>12}{'Exact':>12}{'Error %':>14}")
        for r in R["rows"]:
            e = "n/a" if math.isnan(r["err"]) else f"{r['err']:.4f}%"
            add(f"{'up to D^' + str(r['k']):<14}{fmt(r['cum']):>12}{fmt(R['exact']):>12}{e:>14}")
        add(f"Result (Delta^{cord}): f'({R['x']:.4g}) = {fmt(R['rows'][cord - 1]['cum'])}")
        add(f"Same table, interpolation: f({R['x']:.4g}) = {fmt(R['interp'])}  "
            f"(exact {fmt(R['fexact'])})")

    deriv_block(A, f"TASK 3 - NUMERICAL f'({xa:.4g}) FROM ORIGIN x0 = {A['x0']:.1f}")
    add(f"\nTASK 4 - EXACT DERIVATIVE: f'(x) = 3x^2 + 4x + 1  ->  f'({xa:.4g}) = "
        f"{fmt(A['exact'])};  % error (D^{cord}) = {A['rows'][cord - 1]['err']:.4f}%")
    deriv_block(B, f"TASK 5 - f'({xb:.4g}) FROM THE SAME TABLE (p = {B['p']:.4g})")
    add("\n" + verdict(A, B, cord))
    add("\nSUMMARY PARAGRAPH (draft - rewrite in your own words)")
    add(paragraph(A, B, cord))
    return "\n".join(L)


BG, PANEL, PANEL2 = "#0b1220", "#141e33", "#1b2842"
FG, MUTED = "#e6edf7", "#8da2c0"
ACC, GOOD, WARN, PINK, BAD = "#38bdf8", "#34d399", "#fbbf24", "#f472b6", "#f87171"
MONO = ("Consolas", 10)


class Plot:

    def __init__(s, cv, xr, yr, xt, yt, title, w=560, h=400):
        s.cv, s.xr, s.yr, s.xt, s.yt, s.title = cv, xr, yr, xt, yt, title
        s.w, s.h = w, h
        s.l, s.t, s.r, s.b = 52, 40, 16, 40
        cv.delete("all")

    def X(s, x):
        return s.l + (x - s.xr[0]) / (s.xr[1] - s.xr[0]) * (s.w - s.l - s.r)

    def Y(s, y):
        return s.h - s.b - (y - s.yr[0]) / (s.yr[1] - s.yr[0]) * (s.h - s.t - s.b)

    def axes(s):
        c = s.cv
        c.create_text(s.l, 16, text=s.title, anchor="w", fill=FG, font=("Segoe UI", 11, "bold"))
        for x in s.xt:
            c.create_line(s.X(x), s.t, s.X(x), s.h - s.b, fill="#22304d")
            c.create_text(s.X(x), s.h - s.b + 14, text=f"{x:g}", fill=MUTED, font=("Segoe UI", 8))
        for y in s.yt:
            c.create_line(s.l, s.Y(y), s.w - s.r, s.Y(y), fill="#22304d")
            c.create_text(s.l - 8, s.Y(y), text=f"{y:g}", fill=MUTED, anchor="e", font=("Segoe UI", 8))
        c.create_line(s.l, s.h - s.b, s.w - s.r, s.h - s.b, fill=MUTED)
        c.create_line(s.l, s.t, s.l, s.h - s.b, fill=MUTED)

    def poly(s, pts, **kw):
        s.cv.create_line(*[v for p in pts for v in (s.X(p[0]), s.Y(p[1]))], **kw)

    def dot(s, x, y, r=5, **kw):
        s.cv.create_oval(s.X(x) - r, s.Y(y) - r, s.X(x) + r, s.Y(y) + r, **kw)


class App:
    def __init__(self):
        import tkinter as tk
        from tkinter import ttk, messagebox, filedialog
        self.tk, self.ttk, self.mb, self.fd = tk, ttk, messagebox, filedialog
        self.root = tk.Tk()
        r = self.root
        r.title("Lab 09 · Newton's Forward Differentiation")
        r.geometry("1220x820")
        r.minsize(1100, 720)
        r.configure(bg=BG)

        self.xs = grid()
        self.t = build_table(self.xs)
        self.cord = constant_order(self.t)
        self.var_o = tk.StringVar(value="0.1")
        self.var_a = tk.StringVar(value="0.1")
        self.var_b = tk.StringVar(value="0.0")
        self.var_k = tk.IntVar(value=self.cord)

        self._style()
        self._header()
        self._controls()
        self._tabs()
        self.refresh()

    def _style(self):
        ttk = self.ttk
        s = ttk.Style()
        s.theme_use("clam")
        s.configure(".", background=BG, foreground=FG, fieldbackground=PANEL2)
        s.configure("TNotebook", background=BG, borderwidth=0)
        s.configure("TNotebook.Tab", background=PANEL, foreground=MUTED,
                    padding=(16, 8), font=("Segoe UI", 10, "bold"))
        s.map("TNotebook.Tab", background=[("selected", PANEL2)],
              foreground=[("selected", ACC)])
        s.configure("Treeview", background=PANEL, fieldbackground=PANEL, foreground=FG,
                    rowheight=27, borderwidth=0, font=("Consolas", 10))
        s.configure("Treeview.Heading", background=PANEL2, foreground=ACC,
                    font=("Segoe UI", 9, "bold"), relief="flat")
        s.map("Treeview", background=[("selected", "#274064")])
        s.configure("TCombobox", fieldbackground=PANEL2, background=PANEL2, foreground=FG,
                    selectbackground=PANEL2, selectforeground=FG, arrowcolor=FG)
        s.map("TCombobox", fieldbackground=[("readonly", PANEL2)],
              foreground=[("readonly", FG)], selectbackground=[("readonly", PANEL2)],
              selectforeground=[("readonly", FG)])
        self.root.option_add("*TCombobox*Listbox.background", PANEL2)
        self.root.option_add("*TCombobox*Listbox.foreground", FG)
        s.configure("Accent.TButton", background=ACC, foreground="#06101f",
                    font=("Segoe UI", 10, "bold"), padding=(14, 6))
        s.map("Accent.TButton", background=[("active", "#7dd3fc")])
        s.configure("Ghost.TButton", background=PANEL2, foreground=FG,
                    font=("Segoe UI", 10), padding=(12, 6))
        s.configure("TLabelframe", background=BG, bordercolor="#22304d")
        s.configure("TLabelframe.Label", background=BG, foreground=WARN,
                    font=("Segoe UI", 10, "bold"))
        s.configure("TFrame", background=BG)
        s.configure("Horizontal.TScale", background=BG)

    def _header(self):
        tk = self.tk
        hd = tk.Frame(self.root, bg=BG)
        hd.pack(fill="x", padx=20, pady=(14, 4))
        tk.Label(hd, text="Newton's Forward Differentiation", bg=BG, fg=ACC,
                 font=("Segoe UI", 20, "bold")).pack(side="left")
        tk.Label(hd, text="   f(x) = x³ + 2x² + x + 1   ·   0 ≤ x ≤ 1   ·   h = 0.1",
                 bg=BG, fg=MUTED, font=("Segoe UI", 11)).pack(side="left", pady=(8, 0))
        self.chip = tk.Label(hd, text="", bg=PANEL2, fg=GOOD, padx=12, pady=4,
                             font=("Segoe UI", 10, "bold"))
        self.chip.pack(side="right")

    def _controls(self):
        tk, ttk = self.tk, self.ttk
        bar = tk.Frame(self.root, bg=PANEL, padx=14, pady=10)
        bar.pack(fill="x", padx=20, pady=8)

        def lab(txt):
            tk.Label(bar, text=txt, bg=PANEL, fg=MUTED, font=("Segoe UI", 10)).pack(side="left", padx=(8, 4))

        lab("Table origin x₀")
        cb = ttk.Combobox(bar, textvariable=self.var_o, width=6, state="readonly",
                          values=[f"{x:.1f}" for x in self.xs[:N - KMAX]])
        cb.pack(side="left")
        cb.bind("<<ComboboxSelected>>", lambda e: self.refresh())
        lab("Point A")
        ea = tk.Entry(bar, textvariable=self.var_a, width=7, bg=PANEL2, fg=WARN,
                      insertbackground=FG, relief="flat", font=("Consolas", 11))
        ea.pack(side="left", ipady=3)
        lab("Point B")
        eb = tk.Entry(bar, textvariable=self.var_b, width=7, bg=PANEL2, fg=PINK,
                      insertbackground=FG, relief="flat", font=("Consolas", 11))
        eb.pack(side="left", ipady=3)
        for e in (ea, eb):
            e.bind("<Return>", lambda ev: self.refresh())
        ttk.Button(bar, text="▶  Compute", style="Accent.TButton",
                   command=self.refresh).pack(side="left", padx=14)
        ttk.Button(bar, text="Export report…", style="Ghost.TButton",
                   command=self.export).pack(side="right")
        tk.Label(bar, text="(PDF: A = 0.1 with origin 0.1;  B = 0.0 for Task 5)",
                 bg=PANEL, fg=MUTED, font=("Segoe UI", 9)).pack(side="right", padx=12)

    def _tabs(self):
        tk, ttk = self.tk, self.ttk
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill="both", expand=True, padx=20, pady=(0, 16))
        names = ["1 · Function values", "2 · Difference table", "3-5 · Derivatives",
                 "Graphs", "Report"]
        self.fr = []
        for n in names:
            fr = tk.Frame(self.nb, bg=BG)
            self.nb.add(fr, text=n)
            self.fr.append(fr)

        g = self.fr[3]
        top = tk.Frame(g, bg=BG)
        top.pack(fill="x", padx=10, pady=10)
        tk.Label(top, text="Truncate the series after order:", bg=BG, fg=FG,
                 font=("Segoe UI", 11)).pack(side="left")
        for k in range(1, KMAX + 1):
            tk.Radiobutton(top, text=f"Δ{k if k > 1 else ''}"
                           .replace("Δ2", "Δ²").replace("Δ3", "Δ³").replace("Δ4", "Δ⁴"),
                           variable=self.var_k, value=k, command=self.draw_graphs,
                           bg=BG, fg=FG, selectcolor=PANEL2, activebackground=BG,
                           activeforeground=ACC, font=("Segoe UI", 11, "bold"),
                           indicatoron=False, padx=14, pady=4, bd=0).pack(side="left", padx=4)
        row = tk.Frame(g, bg=BG)
        row.pack()
        self.c1 = tk.Canvas(row, width=560, height=400, bg=PANEL, highlightthickness=0)
        self.c2 = tk.Canvas(row, width=560, height=400, bg=PANEL, highlightthickness=0)
        self.c1.pack(side="left", padx=8)
        self.c2.pack(side="left", padx=8)
        self.ginfo = tk.Label(g, text="", bg=BG, fg=FG, font=("Consolas", 11), justify="left")
        self.ginfo.pack(pady=10)

        self.txt = tk.Text(self.fr[4], bg=PANEL, fg=FG, insertbackground=FG,
                           font=("Consolas", 10), wrap="none", relief="flat", padx=12, pady=10)
        sy = ttk.Scrollbar(self.fr[4], command=self.txt.yview)
        sx = ttk.Scrollbar(self.fr[4], orient="horizontal", command=self.txt.xview)
        self.txt.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.pack(side="right", fill="y")
        sx.pack(side="bottom", fill="x")
        self.txt.pack(fill="both", expand=True)

    def refresh(self):
        try:
            xa, xb = float(self.var_a.get()), float(self.var_b.get())
            if not (0 <= xa <= 1 and 0 <= xb <= 1):
                raise ValueError
            self.idx = round(float(self.var_o.get()) / H)
        except ValueError:
            self.mb.showerror("Invalid input", "Points A and B must be numbers in [0, 1].")
            return
        self.A = evaluate(self.t, self.xs, self.idx, xa)
        self.B = evaluate(self.t, self.xs, self.idx, xb)
        self.chip.config(text=f"Δ{SUP[self.cord][1:]} constant = {fmt(self.t[self.cord][0])}"
                              f"   →   degree {self.cord}")
        self.fill_vals()
        self.fill_diff()
        self.fill_der()
        self.draw_graphs()
        self.txt.delete("1.0", "end")
        self.txt.insert("1.0", build_report(self.idx, xa, xb))

    def _clear(self, fr):
        for w in fr.winfo_children():
            w.destroy()

    def fill_vals(self):
        tk, ttk = self.tk, self.ttk
        fr = self.fr[0]
        self._clear(fr)
        tv = ttk.Treeview(fr, columns=("i", "x", "fx"), show="headings", height=11)
        for c, t, w in (("i", "i", 80), ("x", "xᵢ", 160), ("fx", "f(xᵢ)  (4 d.p.)", 200)):
            tv.heading(c, text=t)
            tv.column(c, width=w, anchor="center")
        for i, x in enumerate(self.xs):
            tv.insert("", "end", values=(i, f"{x:.1f}", f"{f(x):.4f}"))
        tv.pack(side="left", padx=20, pady=20, anchor="n")
        info = tk.Frame(fr, bg=PANEL, padx=20, pady=18)
        info.pack(side="left", padx=10, pady=20, anchor="n", fill="x", expand=True)
        msg = ("Why equally spaced?\n\n"
               f"xᵢ = x₀ + i·h with the same step h = {H} for every i = 0…10\n"
               "(the interval [0,1] is cut into 10 equal pieces).\n\n"
               "Uniform spacing is exactly what lets us use forward differences\n"
               "Δf(xᵢ) = f(xᵢ₊₁) − f(xᵢ) and the single variable p = (x − x₀)/h.\n\n"
               "Remember: ONE difference table gives both the interpolated\n"
               "value and the derivative.")
        tk.Label(info, text=msg, bg=PANEL, fg=FG, justify="left",
                 font=("Segoe UI", 11)).pack(anchor="w")

    def fill_diff(self):
        tk = self.tk
        fr = self.fr[1]
        self._clear(fr)
        grid_fr = tk.Frame(fr, bg=BG)
        grid_fr.pack(padx=16, pady=16, anchor="nw")
        heads = ["x", "f(x)"] + [SUP[k] for k in range(1, N)]
        for c, h in enumerate(heads):
            tk.Label(grid_fr, text=h, width=8, bg=PANEL2, fg=ACC,
                     font=("Segoe UI", 10, "bold")).grid(row=0, column=c, padx=1, pady=1)
        for i in range(N):
            is_o = (i == self.idx)
            tk.Label(grid_fr, text=f"{self.xs[i]:.1f}", width=8, font=MONO,
                     bg=WARN if is_o else PANEL2, fg="#1a1300" if is_o else FG
                     ).grid(row=i + 1, column=0, padx=1, pady=1)
            for k in range(N - i):
                v = self.t[k][i]
                bg, fg = PANEL, FG
                if k == self.cord:
                    bg, fg = "#0e3b2e", GOOD
                elif k > self.cord:
                    bg, fg = BG, "#3a4a68"
                if is_o:
                    bg = "#3b3210" if k <= self.cord else bg
                tk.Label(grid_fr, text=fmt(v), width=8, font=MONO, bg=bg, fg=fg
                         ).grid(row=i + 1, column=k + 1, padx=1, pady=1)
        leg = tk.Frame(fr, bg=BG)
        leg.pack(anchor="w", padx=16, pady=6)
        for col, txt in ((GOOD, f"{SUP[self.cord]} column constant → polynomial degree {self.cord}"),
                         ("#3a4a68", f"orders above {SUP[self.cord]} are all 0"),
                         (WARN, f"row x₀ = {self.xs[self.idx]:.1f}: the origin whose differences feed the formulas")):
            tk.Label(leg, text="■ ", bg=BG, fg=col, font=("Segoe UI", 12)).pack(side="left")
            tk.Label(leg, text=txt + "      ", bg=BG, fg=MUTED, font=("Segoe UI", 10)).pack(side="left")

    def fill_der(self):
        tk, ttk = self.tk, self.ttk
        fr = self.fr[2]
        self._clear(fr)
        cols = ("term", "coef", "cval", "delta", "contrib", "cum", "exact", "err")
        heads = ("Term", "Coefficient C′(p,k)", "Value", "Δᵏf(x₀)", "Contribution",
                 "Cumulative", "Exact f′", "Error %")
        widths = (80, 220, 90, 100, 120, 120, 110, 110)
        for R, name, col in ((self.A, "Point A", WARN), (self.B, "Point B", PINK)):
            lf = ttk.Labelframe(fr, text=f" {name}:  f′({R['x']:.4g})   ·   x₀ = {R['x0']:.1f}   ·   "
                                         f"p = {R['p']:.4g}   ·   exact = {fmt(R['exact'])} ")
            lf.pack(fill="x", padx=14, pady=(10, 2))
            tv = ttk.Treeview(lf, columns=cols, show="headings", height=KMAX)
            for c, h, w in zip(cols, heads, widths):
                tv.heading(c, text=h)
                tv.column(c, width=w, anchor="center")
            tv.tag_configure("keep", foreground=GOOD)
            tv.tag_configure("zero", foreground="#5b6d8c")
            for r in R["rows"]:
                e = "n/a" if math.isnan(r["err"]) else f"{r['err']:.4f}%"
                tv.insert("", "end", tags=("keep" if r["k"] <= self.cord else "zero",),
                          values=(f"Δ{r['k']}" if r["k"] == 1 else SUP[r["k"]],
                                  r["label"], f"{r['coef']:+.4f}", fmt(r["delta"]),
                                  f"{r['term']:+.4f}", fmt(r["cum"]), fmt(R["exact"]), e))
            tv.pack(fill="x", padx=6, pady=6)
            res = R["rows"][self.cord - 1]
            tk.Label(lf, bg=BG, fg=col, font=("Segoe UI", 10, "bold"), anchor="w",
                     text=f"Truncated after {SUP[self.cord]} (constant → higher terms vanish):  "
                          f"f′({R['x']:.4g}) ≈ {fmt(res['cum'])}   (error {res['err']:.4f}%)"
                          f"        Interpolated f({R['x']:.4g}) = {fmt(R['interp'])}"
                          f" (exact {fmt(R['fexact'])})").pack(fill="x", padx=8, pady=(0, 6))
        tk.Label(fr, text=verdict(self.A, self.B, self.cord), bg=PANEL, fg=FG, justify="left",
                 wraplength=1120, font=("Segoe UI", 10), padx=14, pady=10
                 ).pack(fill="x", padx=14, pady=10)

    def draw_graphs(self):
        k = self.var_k.get()
        A, B = self.A, self.B
        p1 = Plot(self.c1, (-0.08, 1.08), (0.4, 5.6), [0, .2, .4, .6, .8, 1.0], [1, 2, 3, 4, 5],
                  "f(x) with tangents (slope from the series)")
        p1.axes()
        p1.poly([(i / 100, f(i / 100)) for i in range(101)], fill=ACC, width=2)
        for x in self.xs:
            p1.dot(x, f(x), 3, fill=MUTED, outline="")
        for R, col in ((A, WARN), (B, PINK)):
            x, y = R["x"], f(R["x"])
            x1, x2 = max(-0.05, x - 0.2), min(1.05, x + 0.2)
            for slope, dash, c in ((R["exact"], (4, 3), "#cbd5e1"),
                                   (R["rows"][k - 1]["cum"], None, col)):
                p1.poly([(x1, y + slope * (x1 - x)), (x2, y + slope * (x2 - x))],
                        fill=c, width=2, dash=dash)
            p1.dot(x, y, 6, fill=col, outline="white")
        self.c1.create_text(p1.w - 14, 20, anchor="e", fill=MUTED, font=("Segoe UI", 8),
                            text="dashed = exact tangent · solid = numerical")

        p2 = Plot(self.c2, (-0.08, 1.08), (0, 8.8), [0, .2, .4, .6, .8, 1.0], [0, 2, 4, 6, 8],
                  "f′(x): exact curve vs numerical points")
        p2.axes()
        p2.poly([(i / 100, df_exact(i / 100)) for i in range(101)], fill=GOOD, width=2)
        for R, col in ((A, WARN), (B, PINK)):
            p2.dot(R["x"], R["exact"], 9, outline="white", fill="")
            p2.dot(R["x"], R["rows"][k - 1]["cum"], 5, fill=col, outline="")
        self.c2.create_text(p2.w - 14, 20, anchor="e", fill=MUTED, font=("Segoe UI", 8),
                            text="ring = exact · dot = numerical")
        ra, rb = A["rows"][k - 1], B["rows"][k - 1]
        self.ginfo.config(
            text=f"Series up to {SUP[k]}:\n"
                 f"  A  f′({A['x']:.4g}) ≈ {fmt(ra['cum'])}   exact {fmt(A['exact'])}   error {ra['err']:.4f}%\n"
                 f"  B  f′({B['x']:.4g}) ≈ {fmt(rb['cum'])}   exact {fmt(B['exact'])}   error {rb['err']:.4f}%")

    def export(self):
        path = self.fd.asksaveasfilename(defaultextension=".txt", initialfile="lab09_report.txt",
                                         filetypes=[("Text", "*.txt")])
        if path:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(self.txt.get("1.0", "end"))
            self.mb.showinfo("Saved", f"Report written to\n{path}")

    def run(self):
        self.root.mainloop()


def main():
    if "--cli" in sys.argv:
        print(build_report())
        return
    try:
        App().run()
    except Exception as e:
        if e.__class__.__name__ == "TclError" or isinstance(e, ImportError):
            print(f"(GUI unavailable: {e}) - falling back to terminal report\n")
            print(build_report())
        else:
            raise


if __name__ == "__main__":
    main()