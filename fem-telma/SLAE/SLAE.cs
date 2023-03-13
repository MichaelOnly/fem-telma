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

            for (var i = 0; i < 4; i++)
            {
                for (var j = 0; j < i; j++)
                {
                    var indexOfElement = Array.IndexOf(M.Ja, element.ElementNumbers[j],
                        M.Ia[element.ElementNumbers[i]],
                        M.Ia[element.ElementNumbers[i] + 1] -
                        M.Ia[element.ElementNumbers[i]]);
                    M.LowTriangleElements[indexOfElement] += locG[i, j];
                    M.UpTriangleElements[indexOfElement] += locG[i, j];
                }
            }
        }
        ApplyBoundary(grid);
    }

    private void ApplyBoundary(Grid.Grid grid)
    {
        var yNumberSegments = grid.YNumberOfSegments;
        var xNumberSegments = grid.XNumberOfSegments;
        var boundaryPoints = new List<int> { grid.Elements[0].ElementNumbers[0] };
        for (var i = 0; i < yNumberSegments; i++)
        {
            boundaryPoints.Add(grid.Elements[i*xNumberSegments].ElementNumbers[2]);
        }

        for (var i = 0; i < xNumberSegments; i++)
        {
            boundaryPoints.Add(grid.Elements[yNumberSegments - 1 + i].ElementNumbers[3]);
        }
        for (int i = yNumberSegments; i >=0; i--)
        {
            boundaryPoints.Add(grid.Elements[i*xNumberSegments + xNumberSegments - 1].ElementNumbers[2]);
        }

        var bigDouble = 1e30;
        foreach (var boundaryPoint in boundaryPoints)
        {
            M.DiagonalElements[boundaryPoint] = bigDouble;
            RHSVector[boundaryPoint] = 0.0;
        }
    }
}