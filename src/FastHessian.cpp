#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// package names: blobomics, blobscopeR

// SolveLinearSystem, DotProduct, findMaximum, interpFeature, and fitQuadrat functions adapted/borrowed from https://github.com/herbertbay/SURF

// Solve the square system of linear equations, Ax=b, where A is given
// in matrix "sq" and b in the vector "solution".  Result is given in
// solution.  Uses Gaussian elimination with pivoting
void SolveLinearSystem(double *solution, double sq[3][3], int size)
{
  int row, col, c, pivot = 0, i;
  double maxc, coef, temp, mult, val;

  // Triangularize the matrix
  for (col = 0; col < size - 1; col++)
  {
    // Pivot row with largest coefficient to top
    maxc = -1.0;
    for (row = col; row < size; row++)
    {
      coef = sq[row][col];
      coef = (coef < 0.0 ? -coef : coef);
      if (coef > maxc)
      {
        maxc = coef;
        pivot = row;
      }
    }
    if (pivot != col)
    {
      // Exchange "pivot" with "col" row (this is no less efficient
      // than having to perform all array accesses indirectly)
      for (i = 0; i < size; i++)
      {
        temp = sq[pivot][i];
        sq[pivot][i] = sq[col][i];
        sq[col][i] = temp;
      }
      temp = solution[pivot];
      solution[pivot] = solution[col];
      solution[col] = temp;
    }
    // Do reduction for this column
    for (row = col + 1; row < size; row++)
    {
      mult = sq[row][col] / sq[col][col];
      for (c = col; c < size; c++) // Could start with c=col+1
        sq[row][c] -= mult * sq[col][c];
      solution[row] -= mult * solution[col];
    }
  }

  // Do back substitution.  Pivoting does not affect solution order
  for (row = size - 1; row >= 0; row--)
  {
    val = solution[row];
    for (col = size - 1; col > row; col--)
      val -= solution[col] * sq[row][col];
    solution[row] = val / sq[row][row];
  }
}

// Return dot product of two vectors with given length
double DotProd(double *v1, double *v2, int len)
{
  int i;
  double sum = 0.0;

  for (i = 0; i < len; i++)
    sum += v1[i] * v2[i];
  return sum;
}

NumericMatrix image(NumericVector intensity, IntegerVector x, IntegerVector y, int xmx, int ymx)
{
  NumericMatrix im(xmx + 1, ymx + 1);
  IntegerVector::iterator it_x, it_y;
  NumericVector::iterator it_intensity;

  // create image matrix
  for (it_x = x.begin(), it_y = y.begin(), it_intensity = intensity.begin(); it_x != x.end(); ++it_x, ++it_y, ++it_intensity)
  {
    im(*it_x, *it_y) = *it_intensity;
  }

  return im;
}

NumericMatrix iimage(NumericMatrix im)
{
  // compute integral
  for (int x = 1; x < im.nrow(); x++)
  {
    for (int y = 1; y < im.ncol(); y++)
    {
      im(x, y) += im(x - 1, y) + im(x, y - 1) - im(x - 1, y - 1);
    }
  }

  return im;
}

IntegerMatrix iimage(IntegerMatrix im)
{
  // compute integral
  for (int x = 1; x < im.nrow(); x++)
  {
    for (int y = 1; y < im.ncol(); y++)
    {
      im(x, y) += im(x - 1, y) + im(x, y - 1) - im(x - 1, y - 1);
    }
  }

  return im;
}

IntegerVector getBoxSizes(int w)
{
  IntegerVector sizes(5);
  // compute boxes (uses integer division)
  sizes(0) = w;
  sizes(1) = w / 2;                     // overall filter offset
  sizes(2) = (int)round(2.0 * w / 9.0); // skip witdh for x & y filters
  sizes(2) = w / 2 - sizes(2);          // box offset for x, y filters - dim1
  sizes(3) = w / 6;                     // box offset for x, y filters - dim2
  sizes(4) = w / 3;                     // box heights for all filters

  return sizes;
}

IntegerVector getOvlap(int w)
{
  IntegerVector ovlap(3);
  IntegerVector sizes = getBoxSizes(w);

  ovlap(0) = (sizes(2) * 2 + 1) * w / 3; // minus x, y filters
  ovlap(1) = ovlap(0) * 2;               // plus x, y filters
  ovlap(2) = sizes(4) * sizes(4) * 2;    // plus and minus of xy filter

  return ovlap;
}

double getSum(NumericMatrix iim, int x1, int y1, int x2, int y2)
{
  return iim(x1 - 1, y1 - 1) + iim(x2, y2) - iim(x2, y1 - 1) - iim(x1 - 1, y2);
}

int getSum(IntegerMatrix iim, int x1, int y1, int x2, int y2)
{
  return iim(x1 - 1, y1 - 1) + iim(x2, y2) - iim(x2, y1 - 1) - iim(x1 - 1, y2);
}

LogicalMatrix getNMSBounds(IntegerMatrix imask)
{
  const int nrow = imask.nrow(), ncol = imask.ncol();
  LogicalMatrix out(nrow - 1, ncol - 1);

  // compute integral
  for (int x = 2; x < nrow - 1; x++)
  {
    for (int y = 2; y < ncol - 1; y++)
    {
      out(x - 1, y - 1) = getSum(imask, x - 1, y - 1, x + 1, y + 1) == 9;
    }
  }

  return out;
}

std::array<NumericMatrix, 6> maskWeights(IntegerMatrix imask, LogicalMatrix mask, IntegerVector sizes, IntegerVector ovlap)
{
  IntegerVector bx(16), by(16);
  const int nrow = imask.nrow(), ncol = imask.ncol();
  int cxxn, cxxp, cyyn, cyyp, cxyn, cxyp;
  NumericMatrix wtxxn(nrow - 1, ncol - 1);
  NumericMatrix wtxxp = clone(wtxxn);
  NumericMatrix wtyyn = clone(wtxxn);
  NumericMatrix wtyyp = clone(wtxxn);
  NumericMatrix wtxyn = clone(wtxxn);
  NumericMatrix wtxyp = clone(wtxxn);

  // compute weights at each loc
  for (int x = 1; x < nrow; x++)
  {
    // vector of box coords
    bx = {x - sizes(3), x + sizes(3), x - sizes(1), x + sizes(1), x - sizes(2), x + sizes(2), x - sizes(2), x + sizes(2), x - sizes(4), x - 1, x + 1, x + sizes(4), x - sizes(4), x - 1, x + 1, x + sizes(4)};
    // deal with out of bounds points
    if ((x - sizes(1)) <= 0 || (x + sizes(1)) >= nrow)
    {
      bx = pmax(pmin(bx, nrow - 1), 1);
    }

    for (int y = 1; y < ncol; y++)
    {
      if (!mask(x, y))
        continue;

      // vector of box coords
      by = {y - sizes(2), y + sizes(2), y - sizes(2), y + sizes(2), y - sizes(3), y + sizes(3), y - sizes(1), y + sizes(1), y + 1, y + sizes(4), y - sizes(4), y - 1, y - sizes(4), y - 1, y + 1, y + sizes(4)};
      // deal with out of bounds points
      if ((y - sizes(1)) <= 0 || (y + sizes(1)) >= ncol)
      {
        by = pmax(pmin(by, ncol - 1), 1);
      }

      // Dxx
      // count
      cxxn = getSum(imask, bx[0], by[0], bx[1], by[1]);
      if (cxxn == 0)
        continue;
      cxxp = getSum(imask, bx[2], by[2], bx[3], by[3]) - cxxn;
      if (cxxp == 0)
        continue;

      // Dyy
      // count
      cyyn = getSum(imask, bx[4], by[4], bx[5], by[5]);
      if (cyyn == 0)
        continue;
      cyyp = getSum(imask, bx[6], by[6], bx[7], by[7]) - cyyn;
      if (cyyp == 0)
        continue;

      // Dxy
      // count
      cxyn = getSum(imask, bx[8], by[8], bx[9], by[9]);
      cxyn += getSum(imask, bx[10], by[10], bx[11], by[11]);
      if (cyyn == 0)
        continue;
      cxyp = getSum(imask, bx[12], by[12], bx[13], by[13]);
      cxyp += getSum(imask, bx[14], by[14], bx[15], by[15]);
      if (cyyp == 0)
        continue;

      // store
      wtxxn(x - 1, y - 1) = -2 * ovlap(0) / (double)cxxn;
      wtxxp(x - 1, y - 1) = ovlap(1) / (double)cxxp;
      wtyyn(x - 1, y - 1) = -2 * ovlap(0) / (double)cyyn;
      wtyyp(x - 1, y - 1) = ovlap(1) / (double)cyyp;
      wtxyn(x - 1, y - 1) = -ovlap(2) / (double)cxyn;
      wtxyp(x - 1, y - 1) = ovlap(2) / (double)cxyp;
    }
  }

  return std::array<NumericMatrix, 6>{wtxxn, wtxxp, wtyyn, wtyyp, wtxyn, wtxyp};
}

NumericMatrix getHessian(NumericMatrix iim, LogicalMatrix mask, IntegerVector sizes, std::array<NumericMatrix, 6> wts)
{
  const int nrow = iim.nrow(), ncol = iim.ncol();
  IntegerVector bx(16), by(16);
  double pos, neg, norm;
  NumericMatrix out(nrow - 1, ncol - 1);
  NumericMatrix wtxxn = wts[0];
  NumericMatrix wtxxp = wts[1];
  NumericMatrix wtyyn = wts[2];
  NumericMatrix wtyyp = wts[3];
  NumericMatrix wtxyn = wts[4];
  NumericMatrix wtxyp = wts[5];

  norm = 81.0 / pow(sizes(0), 4);
  for (int x = 1; x < nrow; x++)
  {
    // vector of box coords
    bx = {x - sizes(3), x + sizes(3), x - sizes(1), x + sizes(1), x - sizes(2), x + sizes(2), x - sizes(2), x + sizes(2), x - sizes(4), x - 1, x + 1, x + sizes(4), x - sizes(4), x - 1, x + 1, x + sizes(4)};
    // deal with out of bounds points
    if ((x - sizes(1)) <= 0 || (x + sizes(1)) >= nrow)
    {
      bx = pmax(pmin(bx, nrow - 1), 1);
    }

    for (int y = 1; y < ncol; y++)
    {
      if (!mask(x, y))
        continue;

      // vector of box coords
      by = {y - sizes(2), y + sizes(2), y - sizes(2), y + sizes(2), y - sizes(3), y + sizes(3), y - sizes(1), y + sizes(1), y + 1, y + sizes(4), y - sizes(4), y - 1, y - sizes(4), y - 1, y + 1, y + sizes(4)};
      // deal with out of bounds points
      if ((y - sizes(1)) <= 0 || (y + sizes(1)) >= ncol)
      {
        by = pmax(pmin(by, ncol - 1), 1);
      }

      // Dxx
      neg = getSum(iim, bx[0], by[0], bx[1], by[1]);
      pos = getSum(iim, bx[2], by[2], bx[3], by[3]) - neg;
      out(x - 1, y - 1) = pos * wtxxp(x - 1, y - 1) + neg * wtxxn(x - 1, y - 1);

      // Dyy
      // weight
      neg = getSum(iim, bx[4], by[4], bx[5], by[5]);
      pos = getSum(iim, bx[6], by[6], bx[7], by[7]) - neg;
      out(x - 1, y - 1) *= pos * wtyyp(x - 1, y - 1) + neg * wtyyn(x - 1, y - 1);

      // Dxy
      // weight
      neg = getSum(iim, bx[8], by[8], bx[9], by[9]);
      neg += getSum(iim, bx[10], by[10], bx[11], by[11]);
      pos = getSum(iim, bx[12], by[12], bx[13], by[13]);
      pos += getSum(iim, bx[14], by[14], bx[15], by[15]);
      out(x - 1, y - 1) -= 0.6 * pow(pos * wtxyp(x - 1, y - 1) + neg * wtxyn(x - 1, y - 1), 2);
      out(x - 1, y - 1) *= norm;
    }
  }

  return out;
}

double fitQuadrat(double *offset, std::array<NumericMatrix, 4> responses, std::array<int, 4> scales, int x, int y, int k)
{
  double g[3], H[3][3];

  // gradients
  offset[0] = g[0] = (responses[k + 1](x, y) - responses[k - 1](x, y)) / 2.0;
  offset[1] = g[1] = (responses[k](x + 1, y) - responses[k](x + 1, y)) / 2.0;
  offset[2] = g[2] = (responses[k](x, y + 1) - responses[k](x, y + 1)) / 2.0;

  // Hessian
  H[0][0] = (responses[k + 1](x, y) - 2.0 * responses[k](x, y) + responses[k - 1](x, y));
  H[1][1] = (responses[k](x + 1, y) - 2.0 * responses[k](x, y) + responses[k](x - 1, y));
  H[2][2] = (responses[k](x, y + 1) - 2.0 * responses[k](x, y) + responses[k](x, y - 1));
  H[0][1] = H[1][0] = ((responses[k + 1](x + 1, y) - responses[k + 1](x - 1, y)) - (responses[k - 1](x + 1, y) - responses[k - 1](x - 1, y))) / 4.0;
  H[0][2] = H[2][0] = ((responses[k + 1](x, y + 1) - responses[k + 1](x, y - 1)) - (responses[k - 1](x, y + 1) - responses[k - 1](x, y - 1))) / 4.0;
  H[1][2] = H[2][1] = ((responses[k](x + 1, y + 1) - responses[k](x + 1, y - 1)) - (responses[k](x - 1, y + 1) - responses[k](x - 1, y - 1))) / 4.0;

  SolveLinearSystem(offset, H, 3);

  return (responses[k](x, y) + 0.5 * DotProd(offset, g, 3));
}

bool interpFeature(double *pt, std::array<NumericMatrix, 4> responses, std::array<int, 4> scales, LogicalMatrix bounds, const double threshold, int movesRemain)
{
  bool change = true;
  double offset[3];

  // Interpolate the detected maximum in order to
  // get a more accurate location
  const double strength = fitQuadrat(offset, responses, scales, pt[1], pt[2], pt[0]);

  if (offset[1] > 0.6 && bounds(pt[1] + 1, pt[2]))
  {
    pt[1]++;
    change = true;
  }
  if (offset[1] < -0.6 && bounds(pt[1] - 1, pt[2]))
  {
    pt[1]--;
    change = true;
  }
  if (offset[2] > 0.6 && bounds(pt[1], pt[2] + 1))
  {
    pt[2]++;
    change = true;
  }
  if (offset[2] < -0.6 && bounds(pt[1], pt[2] - 1))
  {
    pt[2]--;
    change = true;
  }

  if (movesRemain > 0 && change)
  {
    return interpFeature(pt, responses, scales, bounds, threshold, movesRemain - 1);
  }

  // Do not create a keypoint if interpolation still remains far
  // outside expected limits, or if magnitude of peak value is below
  // threshold (i.e., contrast is too low)
  if (isnan(offset[0]) || isnan(offset[1]) || isnan(offset[2]))
    return false;
  if (fabs(offset[0]) > 1.5 || fabs(offset[1]) > 1.5 ||
      fabs(offset[2]) > 1.5 || strength < threshold)
    return false;

  // change points
  pt[0] += offset[0];
  pt[1] += offset[1];
  pt[2] += offset[2];
  pt[3] = strength;

  return true;
}

void findMaximum(vector<list<double>> &df, std::array<NumericMatrix, 4> responses, std::array<int, 4> scales, LogicalMatrix bounds, const double threshold)
{
  const int nrow = responses[0].nrow(), ncol = responses[0].ncol();
  double best;
  double pt[4];
  int r, c, s, ss;
  int dr, dc, ds;
  int cas;

  // find maxima
  for (int k = 1; k < 3; k += 2)
  {
    for (int x = 1; x < nrow - 1; x += 2)
    {
      for (int y = 1; y < ncol - 1; y += 2)
      {
        // skip if not in bounds
        if (!bounds(x, y))
          continue;

        best = responses[k](x, y);
        cas = 0;
        if (responses[k](x + 1, y) > best)
        {
          best = responses[k](x + 1, y);
          cas = 1;
        }
        if (responses[k](x, y + 1) > best)
        {
          best = responses[k](x, y + 1);
          cas = 2;
        }
        if (responses[k](x + 1, y + 1) > best)
        {
          best = responses[k](x + 1, y + 1);
          cas = 3;
        }

        if (responses[k + 1](x, y) > best)
        {
          best = responses[k + 1](x, y);
          cas = 4;
        }
        if (responses[k + 1](x + 1, y) > best)
        {
          best = responses[k + 1](x + 1, y);
          cas = 5;
        }
        if (responses[k + 1](x, y + 1) > best)
        {
          best = responses[k + 1](x, y + 1);
          cas = 6;
        }
        if (responses[k + 1](x + 1, y + 1) > best)
        {
          best = responses[k + 1](x + 1, y + 1);
          cas = 7;
        }
        if (best < threshold)
          continue;
        if (k + 1 == 3 && cas > 3)
          continue; // first 3 = 4 (max scale) - 1

        c = x;
        r = y;
        s = k;
        dc = -1;
        dr = -1;
        ds = -1;
        if (cas != 0)
        {
          if (cas == 1)
          {
            c = x + 1;
            dc = 1;
          }
          else if (cas == 2)
          {
            r = y + 1;
            dr = 1;
          }
          else if (cas == 3)
          {
            c = x + 1;
            r = y + 1;
            dc = 1;
            dr = 1;
          }
          else
          {
            s++;
            ds = 1;
            if (cas == 5)
            {
              c = x + 1;
              dc = 1;
            }
            else if (cas == 6)
            {
              r = y + 1;
              dr = 1;
            }
            else if (cas == 7)
            {
              c = x + 1;
              r = y + 1;
              dc = 1;
              dr = 1;
            }
          }
        }

        ss = s + ds;
        if (best < responses[ss](c - 1, r - dr))
          continue;
        if (best < responses[ss](c, r - dr))
          continue;
        if (best < responses[ss](c + 1, r - dr))
          continue;
        if (best < responses[ss](c - 1, r))
          continue;
        if (best < responses[ss](c, r))
          continue;
        if (best < responses[ss](c + 1, r))
          continue;
        if (best < responses[ss](c - 1, r + dr))
          continue;
        if (best < responses[ss](c, r + dr))
          continue;
        if (best < responses[ss](c + 1, r + dr))
          continue;
        if (best < responses[s](c - 1, r + dr))
          continue;
        if (best < responses[s](c, r + dr))
          continue;
        if (best < responses[s](c + 1, r + dr))
          continue;
        if (best < responses[s](c + dc, r))
          continue;
        if (best < responses[s](c + dc, r - dr))
          continue;
        ss = s - ds;
        if (best < responses[ss](c - 1, r + dr))
          continue;
        if (best < responses[ss](c, r + dr))
          continue;
        if (best < responses[ss](c + 1, r + dr))
          continue;
        if (best < responses[ss](c + dc, r))
          continue;
        if (best < responses[ss](c + dc, r - dr))
          continue;

        // add to results
        pt[0] = s;
        pt[1] = c;
        pt[2] = r;
        pt[3] = responses[s](c, r);
        if (interpFeature(pt, responses, scales, bounds, threshold, 5))
        {
          df[0].push_back(pt[0] / 3.0 * (scales[3] - scales[0]));
          df[1].push_back(pt[1]);
          df[2].push_back(pt[2]);
          df[3].push_back(pt[3]);
        }
      }
    }
  }
}

IntegerVector getAllScales(int octaves)
{
  IntegerVector scales(octaves * 2 + 2);
  scales[0] = 9;
  scales[1] = 15;

  for (int o = 0; o < octaves; o++)
  {
    scales[2 * o + 2] = 3 + pow(2.0, (o)) * (2 + 1) * 6;
    scales[2 * o + 3] = 3 + pow(2.0, (o)) * (3 + 1) * 6;
  }

  return scales;
}

// [[Rcpp::export]]
DataFrame fastHessian(NumericMatrix emat, NumericVector maskvec, IntegerVector x, IntegerVector y, int octaves, int threshold)
{
  DataFrame res;
  vector<list<double>> df(4);
  int xmx = max(x), ymx = max(y);
  NumericMatrix iim;
  IntegerMatrix imask;
  LogicalMatrix mask, bounds;
  IntegerVector ov;
  IntegerVector scales;
  std::vector<std::array<NumericMatrix, 6>> maskwts(octaves * 2 + 2);
  std::vector<IntegerVector> boxSizes(octaves * 2 + 2);
  std::array<NumericMatrix, 4> responses;
  std::array<int, 4> sc;

  // precompute mask and its weights
  imask = (IntegerMatrix)image(maskvec, x, y, xmx, ymx);
  mask = (LogicalMatrix)clone(imask);
  imask = iimage(imask);
  bounds = getNMSBounds(imask); // NMS boundaries

  // prepare scales and weights
  scales = getAllScales(octaves);
  for (int k = 0; k < scales.size(); k++)
  {
    // prepare scales
    boxSizes[k] = getBoxSizes(scales[k]);
    ov = getOvlap(scales[k]);
    maskwts[k] = maskWeights(imask, mask, boxSizes[k], ov);
  }

  for (int g = 0; g < emat.ncol(); g++)
  {
    // create integral image
    iim = iimage(image(emat(_, g), x, y, xmx, ymx));

    // compute responses for first two scales

    responses[1] = getHessian(iim, mask, boxSizes[0], maskwts[0]);
    responses[3] = getHessian(iim, mask, boxSizes[1], maskwts[1]);
    sc[1] = scales[0];
    sc[3] = scales[1];

    for (int o = 0; o < octaves; o++)
    {
      // re-use computations from previous run
      responses[0] = responses[1];
      responses[1] = responses[3];
      responses[2] = getHessian(iim, mask, boxSizes[2 * o + 2], maskwts[2 * o + 2]);
      responses[3] = getHessian(iim, mask, boxSizes[2 * o + 3], maskwts[2 * o + 3]);

      // find maxima
      sc[0] = sc[1];
      sc[1] = sc[3];
      sc[2] = scales[2 * o + 2];
      sc[3] = scales[2 * o + 3];
      findMaximum(df, responses, sc, bounds, threshold);
    }
  }

  res = wrap(df);
  res.attr("names") = CharacterVector::create("scale", "x", "y", "response");
  return res;
}

// [[Rcpp::export]]
NumericMatrix smoothScales(DataFrame features, NumericVector maskvec, IntegerVector x, IntegerVector y)
{
  // precompute mask and its weights
  int xmx = max(x), ymx = max(y);
  LogicalMatrix mask = (LogicalMatrix)image(maskvec, x, y, xmx, ymx);
  IntegerMatrix count = IntegerMatrix(xmx, ymx);
  NumericMatrix scale = NumericMatrix(xmx, ymx);
  NumericVector vx = features["x"], vy = features["y"], vs = features["scale"];
  double cx, cy, cs, cinvdist;

  for (int i = 0; i < features.nrow(); i++)
  {
    // define circle
    cx = vx[i] - 1;
    cy = vy[i] - 1;
    cs = vs[i];

    for (int x = max((int)(cx - cs), 0); x <= min((int)(cx + cs), xmx - 1); x++)
    {
      for (int y = max((int)(cy - cs), 0); y <= min((int)(cy + cs), ymx - 1); y++)
      {
        // distance from edge of circle
        cinvdist = (cs - sqrt(pow(x - cx, 2.0) + pow(y - cy, 2.0))) / cs;
        if (mask(x + 1, y + 1) && cinvdist > 0)
        {
          scale(x, y) += pow(cinvdist, 1.0) * cs;
          count(x, y) += 1;
        }
      }
    }
  }

  // compute average
  for (int x = 0; x < xmx; x++)
  {
    for (int y = 0; y < ymx; y++)
    {
      if (count(x, y) > 0)
        scale(x, y) /= (double)count(x, y);
    }
  }

  return scale;
}
