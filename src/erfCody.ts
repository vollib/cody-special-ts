const A = [3.1611237438705656, 113.864154151050156, 377.485237685302021, 3209.37758913846947, 0.185777706184603153];
const B = [23.6012909523441209, 244.024637934444173, 1282.61652607737228, 2844.23683343917062];
const C = [
  0.564188496988670089, 8.88314979438837594, 66.1191906371416295, 298.635138197400131,
  881.95222124176909, 1712.04761263407058, 2051.07837782607147, 1230.33935479799725,
  2.15311535474403846e-8
];
const D = [
  15.7449261107098347, 117.693950891312499, 537.181101862009858, 1621.38957456669019,
  3290.79923573345963, 4362.61909014324716, 3439.36767414372164, 1230.33935480374942
];
const P = [0.305326634961232344, 0.360344899949804439, 0.125781726111229246, 0.0160837851487422766, 6.58749161529837803e-4, 0.0163153871373020978];
const Q = [2.56852019228982242, 1.87295284992346047, 0.527905102951428412, 0.0605183413124413191, 0.00233520497626869185];

const ZERO = 0.0;
const HALF = 0.5;
const ONE = 1.0;
const TWO = 2.0;
const FOUR = 4.0;
const SQRPI = 0.56418958354775628695;
const THRESH = 0.46875;
const SIXTEEN = 16.0;
const XINF = Number.MAX_VALUE;
const XNEG = -26.628;
const XSMALL = 1.11e-16;
const XBIG = 26.543;
const XHUGE = 6.71e7;
const XMAX = 2.53e307;

function dInt(x: number): number {
  return x > 0 ? Math.floor(x) : -Math.floor(-x);
}

function fixUpForNegativeArgument(jint: 0 | 1 | 2, result: number, x: number): number {
  if (jint === 0) {
    result = (HALF - result) + HALF;
    if (x < ZERO) result = -result;
  } else if (jint === 1) {
    if (x < ZERO) result = TWO - result;
  } else if (x < ZERO) {
    if (x < XNEG) {
      result = XINF;
    } else {
      const ysq = dInt(x * SIXTEEN) / SIXTEEN;
      const del = (x - ysq) * (x + ysq);
      const y = Math.exp(ysq * ysq) * Math.exp(del);
      result = y + y - result;
    }
  }
  return result;
}

export function calerf(x: number, jint: 0 | 1 | 2): number {
  const y = Math.abs(x);
  let ysq: number;
  let xnum: number;
  let xden: number;
  let result: number;

  if (y <= THRESH) {
    ysq = y > XSMALL ? y * y : ZERO;
    xnum = A[4] * ysq;
    xden = ysq;
    for (let i = 0; i < 3; i += 1) {
      xnum = (xnum + A[i]) * ysq;
      xden = (xden + B[i]) * ysq;
    }
    result = x * (xnum + A[3]) / (xden + B[3]);
    if (jint !== 0) result = ONE - result;
    if (jint === 2) result *= Math.exp(ysq);
    return result;
  }

  if (y <= FOUR) {
    xnum = C[8] * y;
    xden = y;
    for (let i = 0; i < 7; i += 1) {
      xnum = (xnum + C[i]) * y;
      xden = (xden + D[i]) * y;
    }
    result = (xnum + C[7]) / (xden + D[7]);
    if (jint !== 2) {
      ysq = dInt(y * SIXTEEN) / SIXTEEN;
      const del = (y - ysq) * (y + ysq);
      result *= Math.exp(-ysq * ysq) * Math.exp(-del);
    }
  } else {
    result = ZERO;
    if (y >= XBIG) {
      if (jint !== 2 || y >= XMAX) return fixUpForNegativeArgument(jint, result, x);
      if (y >= XHUGE) return fixUpForNegativeArgument(jint, SQRPI / y, x);
    }
    ysq = ONE / (y * y);
    xnum = P[5] * ysq;
    xden = ysq;
    for (let i = 0; i < 4; i += 1) {
      xnum = (xnum + P[i]) * ysq;
      xden = (xden + Q[i]) * ysq;
    }
    result = ysq * (xnum + P[4]) / (xden + Q[4]);
    result = (SQRPI - result) / y;
    if (jint !== 2) {
      ysq = dInt(y * SIXTEEN) / SIXTEEN;
      const del = (y - ysq) * (y + ysq);
      result *= Math.exp(-ysq * ysq) * Math.exp(-del);
    }
  }
  return fixUpForNegativeArgument(jint, result, x);
}

export function erfCody(x: number): number {
  return calerf(x, 0);
}

export function erfcCody(x: number): number {
  return calerf(x, 1);
}

export function erfcxCody(x: number): number {
  return calerf(x, 2);
}
