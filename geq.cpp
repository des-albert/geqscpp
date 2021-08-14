#include <cmath>
#include <iostream>
#include <fstream>

#define EXTERN

#include "geq.h"

using namespace std;

extern void bndmat();

extern void eqsil(double *);

extern void splnco(double *);

extern void startt();

extern void compar();

extern void xcur();

extern double gfl(double, double, double, double);

extern void condit(double *, double, double, int, double **, int, int);

extern void saddle();

extern void curf();

extern void plotit();

double **array2d(int, int);

void array2del(double **, int);

int main() {

    double Offset, El, Rxpsn, Zxpsn, tri, elxp, trixp, ang, al1, al2, rc1, rc2, anga, angb, ang1, ang2;
    double rac, zac, exc, Ra[Ng][Mc], Za[Ng][Mc], Ex[Ng][Mc], Rl[Ng][Mc];
    double Rc[Nmax], Zc[Nmax], Rcc[llmax], Zcc[llmax];
    double raxis, zaxis, zdes, alp, xt1, xt2, xt3, curr;
    int ic[Mc], icl, nof, jn, kk, k;

    cout << "Garching Tokamak Equilibrium" << endl;

    int Nexp = 6;
    Mr = 1 + (1 << Nexp);
    Nz = Mr;
    MN = Mr * Nz;
    Mm1 = Mr - 1;
    Nm1 = Nz - 1;
    llp = 2 * (Mr + Nz) - 8;
    pi = 3.141592653589793238462643383279502884197;

    ifstream fin;

    fin.open("solver.dat");

    fin >> Rmin >> Rmax >> Zmin >> Zmax >> error >> meshfg >> mprfg;
    fin >> Rmpl >> Offset >> Apl >> El >> tri >> Rxpsn >> Zxpsn;


    if (meshfg > 0) {
        Rmin = Rmpl - 32.0 * Apl / 20.0;
        Rmax = Rmpl + 32.0 * Apl / 20.0;
        Zmin = Offset - 32.0 * (Offset + fabs(Zxpsn)) / 20.0;
        Zmax = Offset + 32.0 * (Apl * El) / 20.0;
    }

    R = new double[Mr];
    Z = new double[Nz];
    ip = new int[llp];
    jp = new int[llp];
    ityp = new int[Nmax];
    pr = new double[Mr];
    bt2 = new double[Mr];
    cjt = new double[Mr];

    aux = array2d(llp, llp);


    dr = (Rmax - Rmin) / (double) Mm1;
    dz = (Zmax - Zmin) / (double) Nm1;

    for (int i = 0; i < Mr; i++) {
        R[i] = Rmin + dr * (double) i;
    }
    for (int j = 0; j < Nz; j++) {
        Z[j] = Zmin + dz * (double) j;
    }

    alpha = (Rmax - Rmin) / (Rmax + Rmin);
    sh = dz / dr;
    ss = sh * sh;

    cout << "Rmin  = " << Rmin << " Rmax = " << Rmax << " dr = " << dr << endl;
    cout << "Zmin  = " << Zmin << " zmax = " << Zmax << " dz = " << dz << endl;
    cout << "alpha = " << alpha << "   sh = " << sh << " ss = " << ss << endl;

    bndmat();

    /*
        Read in Poloidal Field Coil Data
    */

    k = 0;

    while (true) {
        ic[k] = 0;
        while (true) {
            fin >> rac >> zac >> exc;
            if (rac <= 0.0) goto L10;

            Ra[ic[k]][k] = rac;
            Za[ic[k]][k] = zac;
            Ex[ic[k]][k] = exc;
            Rl[ic[k]][k] = 1.0e-20 * rac;
            ic[k] = ic[k] + 1;

        }
        L10:
        if (ic[k] < 1) goto L20;
        k += 1;
    }

    L20:
    Mmax = k;

    if (mprfg != 0) {
        cout << "Conductor groups available for optimization" << endl;
        for (int j = 0; j < Mmax; j++) {
            cout << " Group : " << j + 1 << endl;
            for (int i = 0; i < ic[j]; i++) {
                printf(" %7.3f   %7.3f   %7.3f \n", Ra[i][j], Za[i][j], Ex[i][j]);
            }
        }
    }

    /*
        Conditions to be satisfied by resulting equilibrium
    */

    elxp = fabs(Offset - Zxpsn) / Apl;
    trixp = (Rmpl - Rxpsn) / Apl;

    for (int j = 0; j < 3; j++) {
        ityp[j] = j;
        Rc[j] = Rxpsn;
        Zc[j] = Zxpsn;
    }

    ityp[3] = 0;
    Rc[3] = Rmpl - Apl;
    Zc[3] = Offset;
    ityp[4] = 0;
    Rc[4] = Rmpl + Apl;
    Zc[4] = Offset;
    ityp[5] = 0;
    Rc[5] = Rmpl - Apl * tri;
    Zc[5] = Offset + Apl * El;

    for (int j = 0; j < 4; j++) {
        ang = (pi * double(j + 1)) / 10.0;
        Rcc[j] = Rmpl + Apl * cos(ang + tri * sin(ang));
        Zcc[j] = Offset + El * Apl * sin(ang);
        Rcc[j + 4] = Rmpl + Apl * cos(ang + pi / 2.0 + tri * sin(ang + pi / 2.0));
        Zcc[j + 4] = Offset + El * Apl * sin(ang + pi / 2.0);
    }

    for (int j = 0; j < 4; j++) {
        al1 = Apl * (((1.0 + trixp) * (1.0 + trixp)) + elxp * elxp) / (2.0 * (1.0 + trixp));
        al2 = Apl * (((1.0 - trixp) * (1.0 - trixp)) + elxp * elxp) / (2.0 * (1.0 - trixp));
        anga = atan(2.0 * elxp * (1.0 + trixp) / (elxp * elxp - (1.0 + trixp) * (1.0 + trixp)));
        angb = atan(2.0 * elxp * (1.0 - trixp) / (elxp * elxp - (1.0 - trixp) * (1.0 - trixp)));
        rc1 = Rmpl + Apl - al1;
        rc2 = Rmpl - Apl + al2;
        ang1 = anga * (j + 1) / 5.0;
        ang2 = angb * (j + 1) / 5.0;
        Rcc[j + 8] = rc1 + al1 * cos(ang1);
        Zcc[j + 8] = Offset - al1 * sin(ang1);
        Rcc[j + 12] = rc2 - al2 * cos(ang2);
        Zcc[j + 12] = Offset - al2 * sin(ang2);
    }

    if (mprfg >= 0) {

        cout << "  Single null case: boundary points " << endl;

        for (int j = 0; j < llmax; j++) {
            printf(" %7.3f   %7.3f\n", Rcc[j], Zcc[j]);
        }
    }

    auto *expsi = new double[MN];
    auto *fool = new double[MN];
    fk = new double[Mc + Nmax];
    double **psiext = array2d(MN, Mc);

    bb = array2d(Nmax, Mc);
    eb = array2d(llmax, Mc);
    cl = array2d(Mc + 1, Mc);


    for (kk = 0; kk < Mmax; kk++) {
        icl = ic[kk];
        for (int i = 0; i < Nz; i++) {
            nof = i * Mr;
            for (int j = 0; j < Mr; j++) {
                jn = nof + j;
                expsi[jn] = 0.;
            }
        }

        for (int i = 0; i < icl; i++) {
            if (!(((Zmax - Za[i][kk]) * (Zmin - Za[i][kk]) <= 0.) && ((Rmax - Ra[i][kk]) * (Rmin - Ra[i][kk]) <= 0.))) {
                for (int k = 0; k < Nz; k += Nm1) {
                    nof = k * Mr;
                    for (int j = 0; j < Mr; j++) {
                        jn = nof + j;
                        expsi[jn] += Ex[i][kk] * gfl(R[j], Ra[i][kk], Z[k] - Za[i][kk], 0.);

                    }
                }

                for (int k = 0; k < Nz; k++) {
                    nof = k * Mr;
                    for (int j = 0; j < Mr; j += Mm1) {
                        jn = nof + j;
                        expsi[jn] += Ex[i][kk] * gfl(R[j], Ra[i][kk], Z[k] - Za[i][kk], 0.);


                    }
                }
            }
        }

        eqsil(expsi);

        for (int i = 0; i < icl; i++) {
            if (!(((Zmax - Za[i][kk]) * (Zmin - Za[i][kk]) > 0.) || ((Rmax - Ra[i][kk]) * (Rmin - Ra[i][kk]) > 0.))) {
                for (int k = 0; k < Nz; k++) {
                    nof = k * Mr;
                    for (int j = 0; j < Mr; j++) {
                        jn = nof + j;
                        expsi[jn] += Ex[i][kk] * gfl(R[j], Ra[i][kk], Z[k] - Za[i][kk], 0.);
                    }
                }
            }
        }

        for (int j = 0; j < MN; j++) {
            psiext[j][kk] = expsi[j];
        }
        /*
            Computation of matrix elements for exact conditions
        */

        splnco(expsi);
        for (int j = 0; j < Nmax; j++) {
            condit(expsi, Rc[j], Zc[j], ityp[j], bb, j, kk);
        }
        for (int j = 0; j < llmax; j++) {
            condit(expsi, Rcc[j], Zcc[j], 0, eb, j, kk);
        }

        /*
            Computation of inductances
        */


        cl[kk][kk] = 0.;
        cl[Mmax][kk] = 0;

        for (int i = 0; i < icl; i++) {
            cl[Mmax][kk] += Ex[i][kk];
            cl[kk][kk] += Ex[i][kk] * Ex[i][kk] * 1.0e6 * (0.58 + log(Ra[i][kk] / Rl[i][kk])) / (2.0 * pi);
        }

        for (int i = 0; i < icl; i++) {
            int ii = i + 2;
            if (ii <= ic[kk]) {
                for (int j = ii - 1; j < icl; j++) {
                    cl[kk][kk] += 2. * Ex[i][kk] * Ex[j][kk] * gfl(Ra[j][kk], Ra[i][kk], Za[j][kk] - Za[i][kk], 0.);
                }
            }
        }
        int lp1 = kk + 1;
        if ((lp1 + 1) <= Mmax) {
            for (int l = lp1; l < Mmax; l++) {
                int icm = ic[l];
                cl[kk][l] = 0.;
                for (int i = 0; i < icl; i++) {
                    for (int j = 0; j < icm; j++) {
                        cl[kk][l] += Ex[i][kk] * Ex[j][l] * gfl(Ra[j][l], Ra[i][kk], Za[j][l] - Za[i][kk], 0.);
                    }
                }
            }
        }
    }

    fin >> icops >> value;
    mpnmax = Mmax + Nmax + 1;
    if (icops >= 1) {
        mpnmax += 1;
    }

    fin >> totcurr >> betapol >> alfac;
    fin >> raxis >> zaxis >> zdes >> alp;
    fin.close();

    jaxis = (int) floor(0.1 + (raxis - R[0]) / dr);
    raxis = R[jaxis];
    naxis = (int) floor(0.1 + (zaxis - Z[0]) / dz);
    zaxis = Z[naxis];
    ndes = (int) floor(0.1 + (fabs(zdes) - Z[0]) / dz);
    if (zdes > 0.) {
        zdes = Z[ndes];
    }

    printf(" Magnetic Axis r = %8.3f  z = %7.3f\n", raxis, zaxis);
    printf(" Rail Limiter  z = %8.3f  alp factor = %10.3e\n", zdes, alp);

    alph = 2. * alp * pi / (llmax * raxis);

    g = new double[Mr * Nz];

    startt();

    com = new double[Mr];
    icycle = 1;
    idecis = 0;

    while (icycle < 20) {
        printf(" ==== Cycle number %2i  ====\n",icycle);

        eqsil(g);

        compar();

        for (int j = 0; j < MN; j++) {
            fool[j] = 0.;
            expsi[j] = g[j];
        }

        splnco(expsi);

        for (int j = 0; j < Nmax; j++) {
            condit(expsi, Rc[j], Zc[j], ityp[j], bb, j, Mmax);
        }
        for (int j = 0; j < llmax; j++) {
            condit(expsi, Rcc[j], Zcc[j], 0, eb, j, Mmax);
        }

        xcur();

        xt1 = 0.;
        xt2 = 0.;
        xt3 = 0.;

        for (int i = 0; i < Mmax; i++) {
            for (int j = 0; j < ic[i]; j++) {
                curr = Ex[j][i] * fk[i];
                xt1 += curr * curr;
                xt2 += fabs(curr);
                xt3 += fabs(curr * Ra[j][i]);
            }
        }

        if (mprfg != 0) {
            printf(" SIG(I**2) = %12.4f  SIG{ABS(I)) = %12.4f  SIG(ABS(RI)) = %12.4f \n", xt1, xt2, xt3);
        }

        for (k = 0; k < Mmax; k++) {
            for (int j = 0; j < MN; j++) {
                fool[j] += fk[k] * psiext[j][k];
                g[j] += fk[k] * psiext[j][k];
            }
        }

        saddle();

        if (mprfg != 0) {
            printf(" Saddle point r = %12.5f  z = %12.5f \n", R[irsp], Z[izsp]);
        }

        if (idecis > 0)
            goto L30;

        curf();

        icycle += 1;
    }

    L30:

    plotit();


    delete[] expsi;
    delete[] fool;
    delete[] g;
    delete[] fk;
    delete[] pr;
    delete[] bt2;
    delete[] cjt;

    array2del(psiext, MN);
    array2del(bb, Nmax);
    array2del(eb, llmax);
    array2del(cl, Mc + 1);

}

double **array2d(int m, int n) {
    auto **ptr = new double *[m];
    for (int i = 0; i < m; i++) {
        ptr[i] = new double[n];
    }
    return ptr;
}

void array2del(double **ptr, int m) {
    for (int i = 0; i < m; i++) {
        delete[] ptr[i];
    }
    delete[] ptr;
}