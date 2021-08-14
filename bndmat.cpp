#include <cmath>
#include <algorithm>
#include "geq.h"

using namespace std;

double gfl(double, double, double, double);

void bndmat() {

    auto *rt = new double[llp];
    auto *zt = new double[llp];

    double ar = 2.0 / (double) Mm1;
    double r0 = 1.0 / alpha;
    double az = sh * ar;
/*
    Index vector interior neighbours to boundary, bottom-top
*/

    int j1 = Mr;
    int j2 = Mr * (Nz - 2);
    int kp = 0;

    for (int i = 2; i <= Mm1; i++) {
        *(ip + kp) = i + j1 - 1;
        *(ip + kp + 1) = i + j2 - 1;
        kp += 2;
    }

/*
    Left - Right
*/

    for (int j = Mr; j <= j2; j += Mr) {
        ip[kp] = j + 1;
        ip[kp + 1] = j + Nm1 - 1;
        kp += 2;
    }

/*
     r - coordinates of boundary points, index vector
*/

    double ra = r0 - 1.;
    double rb = r0 + 1.;
    double za = 0.;
    double zb = az * (double) Nm1;
    int m2 = Mr - 2;
    j2 = Mr * Nm1;

/*
    Bottom - Top
*/
    int nl = 0;
    for (int i = 1; i <= m2; i++) {
        rt[nl] = ra + ar * (double) i;
        rt[nl + 1] = rt[nl];
        zt[nl] = za;
        zt[nl + 1] = zb;
        jp[nl] = i;
        jp[nl + 1] = i + j2;
        nl += 2;
    }

/*
    Left - Right
*/

    int n2 = Nz - 2;
    for (int j = 1; j <= n2; j++) {
        rt[nl] = ra;
        rt[nl + 1] = rb;
        zt[nl] = az * (double) j;
        zt[nl + 1] = zt[nl];
        jp[nl] = j * Mr;
        jp[nl + 1] = (j + 1) * Mr - 1;
        nl += 2;
    }

/*
    Matrix elements

*/
    double arh = 0.5 * ar ;

    for (int i = 0; i < kp; i++) {

/*
    Bottom - Top
*/

        nl = 0;
        for (int j = 0; j < m2; j++) {
            aux[nl][i] = gfl(rt[i], rt[nl], zt[i] - zt[nl], ar) / (sh * rt[nl]);
            aux[nl + 1][i] = gfl(rt[i], rt[nl + 1], zt[i] - zt[nl + 1], ar)/ (sh * rt[nl + 1]);
            nl += 2;
        }
/*
    Left - Right
*/
        for (int j = 0; j < n2; j++) {
            aux[nl][i] = gfl(rt[i], rt[nl], zt[i] - zt[nl], az) * sh / (rt[nl] + arh);
            aux[nl + 1][i] = gfl(rt[i], rt[nl + 1], zt[i] - zt[nl + 1], az)* sh /(rt[nl + 1] - arh);
            nl += 2;
        }

    }

    delete[] rt;
    delete[] zt;

}

double gfl(double rv, double rst, double zv, double del)
{

    double xdl;
    double ak = 4.0 * rv * rst /((rv + rst)*(rv + rst) + zv*zv);
    double x = clamp(((rv - rst)*(rv - rst) + zv*zv)/((rv + rst)*(rv + rst) + zv*zv), 5.e-6, 1.0);
    double y = sqrt(1. - x);

    return sqrt(rv * rst/ak) * ((1. - 0.5 * ak) * comp_ellint_1(y) - comp_ellint_2(y)) / pi;


}

