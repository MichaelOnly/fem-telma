namespace fem_telma.SLAE;

public class SLAE
{
    public Matrix M;
    public double[] RHSVector;
    public double[] Result;

    public SLAE(Grid.Grid grid)
    {
        M = new Matrix(grid);
        RHSVector = new double[M.Size];
        Result = new double[M.Size];
        var g = new[,]
        {
            { 1.0, -1.0 },
            { -1.0, 1.0 }
        };
        var m = new[,]
        {
            { 1.0 / 3.0, 1.0 / 6.0 },
            { 1.0 / 6.0, 1.0 / 3.0 }
        };
        var locG = new double[4, 4];
        foreach (var element in grid.Elements)
        {
            double hx = grid.Points[element.ElementNumbers[1]].X -
                        grid.Points[element.ElementNumbers[0]].X;
            double hy = grid.Points[element.ElementNumbers[2]].Y -
                        grid.Points[element.ElementNumbers[0]].Y;
            var locLambda = 1.0 / element.MagneticPermeability;
            locG[0, 0] = locLambda * (hy * g[0, 0] * m[0, 0] / hx + hx * g[0, 0] * m[0, 0] / hy);
            locG[0, 1] = locLambda * (hy * g[0, 1] * m[0, 0] / hx + hx * g[0, 0] * m[0, 1] / hy);
            locG[0, 2] = locLambda * (hy * g[0, 0] * m[0, 1] / hx + hx * g[0, 1] * m[0, 0] / hy);
            locG[0, 3] = locLambda * (hy * g[0, 1] * m[0, 1] / hx + hx * g[0, 1] * m[0, 1] / hy);
            locG[1, 1] = locLambda * (hy * g[1, 1] * m[0, 0] / hx + hx * g[0, 0] * m[1, 1] / hy);
            locG[1, 2] = locLambda * (hy * g[1, 0] * m[0, 1] / hx + hx * g[0, 1] * m[1, 0] / hy);
            locG[1, 3] = locLambda * (hy * g[1, 1] * m[0, 1] / hx + hx * g[0, 1] * m[1, 1] / hy);
            locG[2, 2] = locLambda * (hy * g[0, 0] * m[1, 1] / hx + hx * g[1, 1] * m[0, 0] / hy);
            locG[2, 3] = locLambda * (hy * g[0, 1] * m[1, 1] / hx + hx * g[1, 1] * m[0, 1] / hy);
            locG[3, 3] = locLambda * (hy * g[1, 1] * m[1, 1] / hx + hx * g[1, 1] * m[1, 1] / hy);
            for (var i = 0; i < 4; i++)
            {
                for (var j = i + 1; j < 4; j++)
                {
                    locG[j, i] = locG[i, j];
                }
            }

            var currentDensity = element.CurrentDensity;
            for (var i = 0; i < 4; i++)
            {
                M.DiagonalElements[element.ElementNumbers[i]] += locG[i, i];
                RHSVector[element.ElementNumbers[i]] += hx * hy * currentDensity / 4.0;
            }
            
            for (var i = 1; i < 4; i++)
            {
                for (var j = 0; j < i; j++)
                {
                    var indexOfElement = GetIndex(element.ElementNumbers[i], element.ElementNumbers[j]);
                    M.LowTriangleElements[indexOfElement] += locG[i, j];
                    M.UpperTriangleElements[indexOfElement] += locG[j, i];
                }
            }
        }
        ApplyBoundary(grid);
    }

    private int GetIndex(int i, int j)
    {
        var index = M.Ia[i];
        while (M.Ja[index] != j)
        {
            index++;
        }

        return index;
    }

    private void ApplyBoundary(Grid.Grid grid)
    {
        var yNumberSegments = grid.YNumberOfSegments;
        var xNumberSegments = grid.XNumberOfSegments;
        var boundaryPoints = new List<int> { grid.Elements[0].ElementNumbers[0] };
        for (var i = 0; i < yNumberSegments; i++)
        {
            boundaryPoints.Add(grid.Elements[i * xNumberSegments].ElementNumbers[2]);
        }
        for (var i = 0; i < xNumberSegments; i++)
        {
            boundaryPoints.Add(grid.Elements[(yNumberSegments - 1)*xNumberSegments + i].ElementNumbers[3]);
        }

        for (int i = yNumberSegments - 1; i >= 0; i--)
        {
            boundaryPoints.Add(grid.Elements[i * xNumberSegments + xNumberSegments - 1].ElementNumbers[1]);
        }

        var bigDouble = 1e30;
        foreach (var boundaryPoint in boundaryPoints)
        {
            M.DiagonalElements[boundaryPoint] = bigDouble;
            RHSVector[boundaryPoint] = 0.0;
        }
    }

    private double[] LUForwardProp(Matrix mx, double[] rhsVector)
    {
        var result = new double[rhsVector.Length];
        rhsVector.AsSpan().CopyTo(result);
        for (var i = 0; i < result.Length; i++)
        {
            var sum = 0.0;
            for (var j = mx.Ia[i]; j < mx.Ia[i + 1]; j++)
                sum += mx.LowTriangleElements[j] * result[mx.Ja[j]];
            result[i] -= sum;
            result[i] /= mx.DiagonalElements[i];
        }

        return result;
    }

    private double[] LUBackwardProp(Matrix mx, double[] rhsVector)
    {
        var result = new double[rhsVector.Length];
        rhsVector.AsSpan().CopyTo(result);
        for (var i = mx.Size - 1; i >= 0; i--)
        {
            result[i] /= mx.DiagonalElements[i];
            for (var j = mx.Ia[i]; j < mx.Ia[i + 1]; j++)
                result[mx.Ja[j]] -= mx.UpperTriangleElements[j] * result[i];
        }

        return result;
    }

    private double ScalarProd(double[] x, double[] y)
    {
        var result = 0.0;
        for (var i = 0; i < x.Length; i++)
            result += x[i] * y[i];

        return result;
    }

    public void SolveWithLOSPrecond(double eps, int maxIter)
    {
        var factorizedM = M.LUsqFactorization();
        var k = 0;
        var buf = M.MultMatrixOnVector(Result);
        for (var i = 0; i < buf.Length; i++)
            buf[i] = RHSVector[i] - buf[i];
        var r = LUForwardProp(factorizedM, buf);
        var error = ScalarProd(r, r);
        var error1 = 3*error / 2;
        var z = LUBackwardProp(factorizedM, r);
        buf = M.MultMatrixOnVector(z);
        var p = LUForwardProp(factorizedM, buf);
        while (error > eps && k < maxIter && Math.Abs(error - error1) > eps)
        {
            var pp = ScalarProd(p, p);
            var pr = ScalarProd(p, r);
            var alpha = pr / pp;
            error1 = error;
            error -= alpha * alpha * pp;
            error = Math.Abs(error);
            for (var i = 0; i < M.Size; i++)
            {
                Result[i] += alpha * z[i];
                r[i] -= alpha * p[i];
            }

            var Ur = LUBackwardProp(factorizedM, r);
            buf = M.MultMatrixOnVector(Ur);
            buf = LUForwardProp(factorizedM, buf);
            var betta = -(ScalarProd(p, buf) / pp);
            for (var i = 0; i < M.Size; i++)
            {
                z[i] = Ur[i] + betta * z[i];
                p[i] = buf[i] + betta * p[i];
            }
            k++;
        }
        Console.WriteLine("iter: " + k.ToString()+" Error: "+ error.ToString());
    }
}