import numpy as np

# import matplotlib as mp

# -----------------------Constantes----------------------------
resist = 100
cap = 1e-6
rl = 1
ra = 1000
ind = 2e-3
kc1 = 1
kc2 = 0

CARD = 5
NQ = 21
INCR = 1
NPG = 200
# NA = NPG / INCR
NA = int(200)

# -------------------------------------------------------------
var = 5
sinze = 2
XNmin = 0
YNmin = 0
XNmax = 1
YNmax = 1
ndv = 0
ndh = 0

matrix = np.zeros([sinze, sinze], dtype=float)
vector = np.zeros([sinze], dtype=float)


def list_generator(x: []):
    if not isinstance(x, list):
        raise TypeError
    for i in range(1):
        x.append(i)


def cdx(W, Wmin, Wmax, nd):
    result = (((XNmax - XNmin) / (Wmax - Wmin) * (W - Wmin) + XNmin) * (nd - 1) + 0.5)
    result = int(result)
    return result


def cdy(W, Wmin, Wmax, nd):
    result = (nd - 1 - ((YNmax - YNmin) / (Wmax - Wmin) * (W - Wmin) + YNmin) * (nd - 1) + 0.5)
    result = int(result)
    return result


# def imp_escala(Xinicio, Xincr, Xfim, Yinicio,
#                Yincr, Yfim, XWmin, YWmin,
#                XWmax, YWmax, Xdesloc, Ydesloc):

def imp_escala(Xinicio, Xincr, Xfim, Yinicio,
               Yincr, Yfim, XWmin, YWmin,
               XWmax, YWmax):
    ne = int(((Xfim - Xinicio) / Xincr) + 1)

    for i in range(ne):
        # if Xinicio != 0:
            # x1 = cdx(Xinicio, XWmin, XWmax, ndh)
            # y1 = cdy(0, YWmin, YWmax, ndv)
            # line(x1,y1,x1,y1+6);
            # outtextxy(x1-Xdesloc,y1+9,gcvt(Xinicio,3," "))
        Xinicio += Xincr

    ne = int(((Yfim - Yinicio) / Yincr) + 1)

    for i in range(ne):
        # if Yinicio != 0:
            # x1 = cdx(0, XWmin, XWmax, ndh)
            # y1 = cdy(Yinicio, YWmin, YWmax, ndv)
            # line(x1,y1,x1-5,y1);
            # outtextxy(x1-Ydesloc,y1-3,gcvt(Yinicio,3," ")
        Yinicio += Yincr


# -------------------------------------------------------------------------------------

# ------------------------SUBROTINA QUE MULTIPLICA DUAS MATRIZES---------------------
def multmm(a1, a2, a):
    for i in range(sinze):
        for j in range(sinze):

            a[i, j] = float(0)
            for k in range(sinze):
                a[i, j] = a[i, j] + (a1[i, j] * a2[i, j])


# -----------------------------------------------------------------------------------

# ------------SUBROTINA QUE MULTIPLICA UMA MATRIZ POR UM VETOR-----------------------
def multmv(w1: [], v1: [], v: []):

    w1 = np.zeros([sinze, sinze], dtype=float)
    v0 = np.zeros([sinze], dtype=float)
    v1 = np.zeros([sinze], dtype=float)
    for i in range(sinze):
        v[i] = float(0)
        for j in range(sinze):
            v[i] = v0[i] + (w1[i, j] * v1[j])


# -----------------------------------------------------------------------------------

def rungeKutta(A, B, H, U, X):
    w2 = np.zeros([sinze, sinze], dtype=float)
    w3 = np.zeros([sinze, sinze], dtype=float)
    w4 = np.zeros([sinze, sinze], dtype=float)

    k1 = np.zeros([sinze, sinze], dtype=float)
    k2 = np.zeros([sinze, sinze], dtype=float)
    k3 = np.zeros([sinze, sinze], dtype=float)
    k4 = np.zeros([sinze, sinze], dtype=float)

    w5 = np.zeros([sinze], dtype=float)
    r1 = np.zeros([sinze], dtype=float)
    r2 = np.zeros([sinze], dtype=float)

    identity = np.identity(sinze, dtype=float)

    multmm(A, A, w2)
    multmm(A, w2, w3)
    multmm(w2, w2, w4)
    multmv(B, U, w5)

    print(A, w2, w3)

    for i in range(sinze):
        for j in range(sinze):
            k1[i, j] = identity[i, j] + (H * A[i, j]) + (0.5 * (H ** 2) * w2[i, j]) + \
                       (0.167 * (H ** 3) * w3[i, j]) + (0.0417 * (H ** 4) * w4[i, j])

            k2[i, j] = identity[i, j] + (H * A[i, j]) + (0.5 * (H ** 2) * w2[i, j]) + \
                       (0.25 * (H ** 3) * w3[i, j])

            k3[i, j] = identity[i, j] + (0.5 * H * A[i, j]) + (0.125 * (H ** 2) * w2[i][j])

            k4[i, j] = 0.167 * H * identity[i, j] + (0.167 * H * k2[i, j]) + (0.667 * H * k3[i, j])

    multmv(k1, X, r1)
    multmv(k4, w5, r2)
    for i in range(sinze):
        X[i] = r1[i] + r2[i]


def fuzificacao(vnf: float, vqv, pos):
    # difant = difatual = float(0)

    pos = np.zeros([1])

    difant = np.abs(vnf - vqv[0])  # difant = fabs(vnf - vqv[0]);
    i = 0
    while pos[0] == -1 & i < NQ:

        difatual = np.abs(vnf - vqv[i + 1])  # difatual = fabs(vnf - vqv[i+1]);

        if difant <= difatual:
            pos[0] = i

        difant = difatual
        i = i + 1

        if i == NQ:
            pos[0] = i - 1


def obtem_regras_inferencia(cFuzzy, ymi, xmi, nri, ri):
    minimo = xmi
    if ymi < minimo:
        minimo = ymi
    for i in range(NQ):
        num = int(cFuzzy)
        if minimo < num:
            cFuzzy = minimo
        ri[nri, i] = cFuzzy


def uniao_das_regras(nri: int, ri: [], r: [NQ]):
    for j in range(NQ):
        maximo = ri[0][j]
        for i in range(nri):
            if ri[i][j] > maximo:
                maximo = ri[i, j]
        r[j] = maximo


def defuzificacao(r: [NQ], vdelta: [NQ], dd):
    soma = float(0.00)
    produto = float(0)

    dd = np.zeros([1])

    for i in range(NQ):
        soma += r[i]
        produto += r[i] * vdelta[i]
    dd[0] = produto / soma


def controlador_fuzzy(e, dfi, dd):
    v_delta = np.array([-30, -27, -24, -21, -18, -15, -12, -9, -6,
                        -3, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30], dtype=int)

    v_e = np.array([-20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2,
                    4, 6, 8, 10, 12, 14, 16, 18, 20], dtype=int)

    v_dfi = np.array([-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
                      0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], dtype=float)

    cFuzzy = np.array([[1, 1, 0.8, 0.5, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0.3, 0.5, 0.8, 1, 1, 1, 0.8, 0.5, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.8, 1, 1, 1, 0.8, 0.5, 0.3, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.8, 1, 1, 1, 0.8, 0.5, 0.3, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.8, 1, 1]], dtype=float)

    struct = {'mi': 20,
              'nc': 0}
    x = []
    y = []
    for i in range(CARD):
        x.append(struct)
        y.append(struct)

    br = np.array([[0, 2, 4, 4, 4],
                   [0, 1, 3, 4, 4],
                   [0, 1, 2, 3, 4],
                   [0, 0, 1, 3, 4],
                   [0, 0, 0, 2, 4]])

    ie = np.zeros([1])
    idfi = np.zeros([1])

    npe = 0
    ri = np.zeros([4, NQ], dtype=float)
    r = np.zeros([NQ], dtype=float)

    fuzificacao(e, v_e, ie)
    fuzificacao(dfi, v_dfi, idfi)

    for i in range(CARD):
        ie_0 = int(ie[0])
        if cFuzzy[i, ie_0] != 0:
            x[npe]['nc'] = i
            x[npe]['mi'] = cFuzzy[i, ie_0]
            npe += 1
    npdfi = 0
    for i in range(CARD):
        idfi_0 = int(idfi[0])
        if cFuzzy[i][idfi_0] != 0.0:
            y[npdfi]['nc'] = i
            y[npdfi]['mi'] = cFuzzy[i, idfi_0]
            npdfi += 1
    nri = 0
    for i in range(CARD):
        for j in range(npdfi):
            ncasaida = br[y[j]['nc'], x[i]['nc']]
            obtem_regras_inferencia(cFuzzy[ncasaida, 0], y[j]['mi'], x[i]['mi'], nri, ri)
    uniao_das_regras(nri, ri, r)
    defuzificacao(r, v_delta, dd)
    for i in range(CARD):
        print(y[i]['mi'])


def dinamica_Sistema(h, dt, v: [], u: []):
    const = 1
    a = np.zeros([sinze, sinze], dtype=float)
    b = np.zeros([sinze, const], dtype=float)
    vu = np.array([0, 0], dtype=float)
    x = np.zeros([sinze], dtype=float)

    a[0, 0] = -1 / (resist * cap)
    a[0, 1] = kc1 / cap
    a[1, 0] = -kc1 / ind
    a[1, 1] = -((kc1 * rl) + (kc2 * ra)) / ind

    b[1, 0] = kc1 / ind
    vu[1] = float(0)
    vu[0] = 100 * np.sin(2500 * dt / 100000)

    print("Vector in position 0 - vu[0]: {}".format(vu[0]))

    rungeKutta(a, b, h, vu, x)
    u[0] = x[0]
    v[0] = x[1]


def sistema():
    dfi = np.zeros([1])
    dd = np.zeros([1])
    fi = np.zeros([1])

    h = 1e-6
    ve = np.zeros([NA], dtype=float)
    vdd = np.zeros([NA], dtype=float)
    vda = np.zeros([NA], dtype=float)
    vdfi = np.zeros([NA], dtype=float)
    vfi = np.zeros([NA], dtype=float)
    vfic = np.zeros([NA], dtype=float)

    for t in range(NA):
        # print("Apresentação numero:", t)
        fic = 20
        vfic[t] = -100 * np.sin(2500 * t / 100000)
        e = fic - fi[0]
        ve[t] = e

        dfi_0 = dfi[0]
        controlador_fuzzy(e, dfi_0, dd)

        vdd[t] = dd[0]
        dinamica_Sistema(h, t, dfi, fi)
        vdfi[t] = dfi
        vfi[t] = fi[0]
    # print("  ----- A G U A R D E -----")
    # funcao_grafica(vfic, ve, vdd, vda, vdfi, vfi)
