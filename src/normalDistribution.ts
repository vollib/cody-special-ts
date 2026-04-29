import { DBL_EPSILON, DBL_MAX, ONE_OVER_SQRT_TWO, ONE_OVER_SQRT_TWO_PI, SQRT_TWO } from "./constants.js";
import { erfcCody } from "./erfCody.js";

const normCdfAsymptoticExpansionFirstThreshold = -10.0;
const normCdfAsymptoticExpansionSecondThreshold = -1 / Math.sqrt(DBL_EPSILON);
const uMax = 0.3413447460685429;

export function normPdf(x: number): number {
  return ONE_OVER_SQRT_TWO_PI * Math.exp(-0.5 * x * x);
}

export function normCdf(z: number): number {
  if (z <= normCdfAsymptoticExpansionFirstThreshold) {
    let sum = 1;
    if (z >= normCdfAsymptoticExpansionSecondThreshold) {
      const zsqr = z * z;
      let i = 1;
      let g = 1;
      let a = DBL_MAX;
      let lasta: number;
      do {
        lasta = a;
        const x = (4 * i - 3) / zsqr;
        const y = x * ((4 * i - 1) / zsqr);
        a = g * (x - y);
        sum -= a;
        g *= y;
        i += 1;
        a = Math.abs(a);
      } while (lasta > a && a >= Math.abs(sum * DBL_EPSILON));
    }
    return -normPdf(z) * sum / z;
  }
  return 0.5 * erfcCody(-z * ONE_OVER_SQRT_TWO);
}

export function inverseNormCdfForLowProbabilities(p: number): number {
  const r = Math.sqrt(-Math.log(p));
  if (r < 6.7) {
    if (r < 3.41) {
      if (r < 2.05) {
        return (3.691562302945566191 + r * (4.7170590600740689449e1 + r * (6.5451292110261454609e1 + r * (-7.4594687726045926821e1 + r * (-8.3383894003636969722e1 - 1.3054072340494093704e1 * r))))) /
          (1 + r * (2.0837211328697753726e1 + r * (7.1813812182579255459e1 + r * (5.9270122556046077717e1 + r * (9.2216887978737432303 + 1.8295174852053530579e-4 * r)))));
      }
      return (3.2340179116317970288 + r * (1.449177828689122096e1 + r * (6.8397370256591532878e-1 + r * (-1.81254427791789183e1 + r * (-1.005916339568646151e1 - 1.2013147879435525574 * r))))) /
        (1 + r * (8.8820931773304337525 + r * (1.4656370665176799712e1 + r * (7.1369811056109768745 + r * (8.4884892199149255469e-1 + 1.0957576098829595323e-5 * r)))));
    }
    return (3.1252235780087584807 + r * (9.9483724317036560676 + r * (-5.1633929115525534628 + r * (-1.1070534689309368061e1 + r * (-2.8699061335882526744 - 1.5414319494013597492e-1 * r))))) /
      (1 + r * (7.076769154309171622 + r * (8.1086341122361532407 + r * (2.0307076064309043613 + r * (1.0897972234131828901e-1 + 1.3565983564441297634e-7 * r)))));
  }
  if (r < 12.9) {
    return (2.6161264950897283681 + r * (2.250881388987032271 + r * (-3.688196041019692267 + r * (-2.9644251353150605663 + r * (-4.7595169546783216436e-1 - 1.612303318390145052e-2 * r))))) /
      (1 + r * (3.2517455169035921495 + r * (2.1282030272153188194 + r * (3.3663746405626400164e-1 + r * (1.1400087282177594359e-2 + 3.0848093570966787291e-9 * r)))));
  }
  return (2.3226849047872302955 + r * (-4.2799650734502094297e-2 + r * (-2.5894451568465728432 + r * (-8.6385181219213758847e-1 + r * (-6.5127593753781672404e-2 - 1.0566357727202585402e-3 * r))))) /
    (1 + r * (1.9361316119254412206 + r * (6.1320841329197493341e-1 + r * (4.6054974512474443189e-2 + r * (7.471447992167225483e-4 + 2.3135343206304887818e-11 * r)))));
}

export function inverseNormCdfmHalfForMidrangeProbabilities(u: number): number {
  const s = uMax * uMax - u * u;
  return u * ((2.92958954698308805 + s * (5.0260572167303103e1 + s * (3.01870541922933937e2 + s * (7.4997781456657924e2 + s * (6.90489242061408612e2 + s * (1.34233243502653864e2 - 7.58939881401259242 * s)))))) /
    (1 + s * (1.8918538074574598e1 + s * (1.29404120448755281e2 + s * (3.86821208540417453e2 + s * (4.79123914509756757e2 + 1.79227008508102628e2 * s))))));
}

export function inverseNormCdf(p: number): number {
  const u = p - 0.5;
  if (Math.abs(u) < uMax) return inverseNormCdfmHalfForMidrangeProbabilities(u);
  return u > 0 ? -inverseNormCdfForLowProbabilities(1 - p) : inverseNormCdfForLowProbabilities(p);
}

export function erfInv(e: number): number {
  if (Math.abs(e) < 2 * uMax) {
    return inverseNormCdfmHalfForMidrangeProbabilities(0.5 * e) / SQRT_TWO;
  }
  return (e < 0
    ? inverseNormCdfForLowProbabilities(0.5 * e + 0.5)
    : -inverseNormCdfForLowProbabilities(-0.5 * e + 0.5)) / SQRT_TWO;
}
