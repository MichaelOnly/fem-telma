using System.Numerics;

namespace fem_telma;

public class ResultFunction
{
    private double[] Weights;
    private Grid.Grid Grid;

    public ResultFunction(double[] weights, Grid.Grid grid)
    {
        Weights = new double[weights.Length];
        weights.AsSpan().CopyTo(Weights);
        Grid = grid;
    }

    public double GerResultFunction(double x, double y)
    {
        var xNumberElement = 0;
        var yNumberElement = 0;
        while (!(Grid.XGrid[xNumberElement] <= x && x <= Grid.XGrid[xNumberElement+1]))
        {
            xNumberElement++;
        }

        while (!(Grid.YGrid[yNumberElement] <= y && y <= Grid.YGrid[yNumberElement+1]))
        {
            yNumberElement++;
        }

        var elementNumber = Grid.XNumberOfSegments * yNumberElement + xNumberElement;
        var x0 = (double)Grid.Points[Grid.Elements[elementNumber].ElementNumbers[0]].X;
        var x1 = (double)Grid.Points[Grid.Elements[elementNumber].ElementNumbers[1]].X;
        var y0 = (double)Grid.Points[Grid.Elements[elementNumber].ElementNumbers[0]].Y;
        var y1 = (double)Grid.Points[Grid.Elements[elementNumber].ElementNumbers[2]].Y;
        var hx = x1 - x0;
        var hy = y1 - y0;
        var X1 = (x1 - x) / hx;
        var X2 = (x - x0) / hx;
        var Y1 = (y1 - y) / hy;
        var Y2 = (y - y0) / hy;
        var firstFunc = X1 * Y1;
        var secondFunc = X2 * Y1;
        var thirdFunc = X1 * Y2;
        var fourthFunc = X2 * Y2;
        var result = Weights[Grid.Elements[elementNumber].ElementNumbers[0]] * firstFunc +
                     Weights[Grid.Elements[elementNumber].ElementNumbers[1]] * secondFunc +
                     Weights[Grid.Elements[elementNumber].ElementNumbers[2]] * thirdFunc +
                     Weights[Grid.Elements[elementNumber].ElementNumbers[3]] * fourthFunc;
        return result;
    }

    public double GetAbsB(double x, double y)
    {
        var deltax = 1e-12;
        var deltay = 1e-12;
        var deltaAx = (GerResultFunction(x + deltax, y) - GerResultFunction(x, y)) / deltax;
        var deltaAy = (GerResultFunction(x, y + deltax) - GerResultFunction(x, y)) / deltay;
        var absB = Math.Sqrt(deltaAx * deltaAx + deltaAy * deltaAy);
        return absB;
    }
}