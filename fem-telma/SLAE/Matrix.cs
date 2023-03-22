namespace fem_telma.SLAE;

public class Matrix
{
    public int[] Ia;
    public int[] Ja;
    public double[] LowTriangleElements;
    public double[] UpperTriangleElements;
    public double[] DiagonalElements;
    public int Size;

    public Matrix(Grid.Grid grid)
    {
        Size = grid.Points.Length;
        DiagonalElements = new double[Size];
        Ia = new int[Size + 1];
        Ia[0] = 0;
        Ia[1] = 0;
        for (var i = 1; i < Size + 1; i++)
        {
            Ia[i + 1] = Ia[i] + grid.PointsBounds[i].Count;
        }

        Ja = new int[Ia.Last()];
        for (var i = 1; i < Size + 1; i++)
        {
            var stringSize = Ia[i + 1] - Ia[i];
            for (var j = 0; j < stringSize; j++)
            {
                Ja[Ia[i] + j] = grid.PointsBounds[i][j];
            }
        }
        LowTriangleElements = new double[Ia.Last()];
        UpperTriangleElements = new double[Ia.Last()];
    }

    public Matrix()
    {
        Size = 1;
        DiagonalElements = new double[Size];
        Ia = new int[Size + 1];
        Ja = new int[Size];
        LowTriangleElements = new double[Size];
        UpperTriangleElements = new double[Size];
    }

    public double[] MultMatrixOnVector(double[] vec)
    {
        var result = new double[vec.Length];
        for (var i = 0; i < vec.Length; i++)
        {
            result[i] = DiagonalElements[i] * vec[i];
            for (var j = Ia[i]; j < Ia[i+1]; j++)
            {
                result[i] += LowTriangleElements[j] * vec[Ja[j]];
                result[Ja[j]] += UpperTriangleElements[j] * vec[i];
            }
        }

        return result;
    }

    private Matrix CopyMatrix()
    {
        var mxCopy = new Matrix();
        Ia.AsSpan().CopyTo(mxCopy.Ia);
        mxCopy.Size = Size;
        Ja.AsSpan().CopyTo(mxCopy.Ja);
        DiagonalElements.AsSpan().CopyTo(mxCopy.DiagonalElements);
        LowTriangleElements.AsSpan().CopyTo(mxCopy.LowTriangleElements);
        UpperTriangleElements.AsSpan().CopyTo(mxCopy.UpperTriangleElements);
        return mxCopy;
    }
    public Matrix LUsqFactorization()
    {
        var factorizedMx = CopyMatrix();
        for (var i = 0; i < factorizedMx.Size; i++)
        {
            var sumdi = 0.0;
            var i0 = factorizedMx.Ia[i];
            var i1 = factorizedMx.Ia[i + 1];

            for (var k = i0; k < i1; k++)
            {
                var j = factorizedMx.Ja[k];
                var j0 = factorizedMx.Ia[j];
                var j1 = factorizedMx.Ia[j + 1];

                var ik = i0;
                var kj = j0;

                var suml = 0.0;
                var sumu = 0.0;
                while (ik < k)
                {
                    if (factorizedMx.Ja[ik] == factorizedMx.Ja[kj])
                    {
                        suml += factorizedMx.LowTriangleElements[ik] * factorizedMx.UpperTriangleElements[kj];
                        sumu += factorizedMx.UpperTriangleElements[ik] * factorizedMx.LowTriangleElements[kj];
                        ik++;
                        kj++;
                    }
                    else if (factorizedMx.Ja[ik] > factorizedMx.Ja[kj]) kj++;
                    else ik++;
                }

                factorizedMx.LowTriangleElements[k] =
                    (factorizedMx.LowTriangleElements[k] - suml) / factorizedMx.DiagonalElements[j];
                factorizedMx.UpperTriangleElements[k] =
                    (factorizedMx.UpperTriangleElements[k] - sumu) / factorizedMx.DiagonalElements[j];
                sumdi += factorizedMx.LowTriangleElements[k] * factorizedMx.UpperTriangleElements[k];
            }
            factorizedMx.DiagonalElements[i] = Double.Sqrt(factorizedMx.DiagonalElements[i] - sumdi); 
        }

        return factorizedMx;
    }
}