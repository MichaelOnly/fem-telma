namespace fem_telma.SLAE;

public class Matrix
{
    public int[] Ia;
    public int[] Ja;
    public double[] LowTriangleElements;
    public double[] UpTriangleElements;
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
        UpTriangleElements = new double[Ia.Last()];
    }
}