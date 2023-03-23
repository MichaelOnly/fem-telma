using System.Numerics;

namespace fem_telma.Grid;

public class Grid
{
    public readonly Vector2[] Points;
    public readonly Element[] Elements;
    public readonly List<int>[] PointsBounds;
    public readonly int XNumberOfSegments;
    public readonly int YNumberOfSegments;
    public readonly float[] XGrid;
    public readonly float[] YGrid;

    private static List<float> Make1DGrid(float[] points, float[] steps, float[] relaxRatios,
        bool[] isReversedRelaxRatio)
    {
        var grid = new List<float> { points[0] };
        for (var i = 0; i < points.Length - 1; i++)
        {
            var step = steps[i];
            var relaxRatio = relaxRatios[i];
            if (Double.Abs(relaxRatio - 1.0) < 1e-14)
            {
                var numberSegments = (int)((points[i + 1] - points[i]) / step);
                step = (points[i + 1] - points[i]) / numberSegments;
                for (var j = 1; j < numberSegments + 1; j++)
                {
                    grid.Add(points[i] + step * j);
                }
            }
            else
            {
                var x = points[i];
                var numberSegments = 0;
                while (x < points[i + 1])
                {
                    x += step * (float)Math.Pow(relaxRatio, numberSegments);
                    numberSegments++;
                }

                if (isReversedRelaxRatio[i])
                {
                    relaxRatio = 1 / relaxRatio;
                }

                step = (points[i + 1] - points[i]) /
                       ((1 - (float)Math.Pow(relaxRatio, numberSegments)) / (1 - relaxRatio));
                for (var j = 1; j <= numberSegments; j++)
                {
                    grid.Add(points[i] + step * ((1 - (float)Math.Pow(relaxRatio, j)) / (1 - relaxRatio)));
                }
            }
        }

        return grid;
    }

    public Grid(Areas areas)
    {
        var xGrid = Make1DGrid(areas.Xs, areas.XSteps, areas.XRelaxRatios, areas.IsReversedXRelaxRatio).ToArray();
        var yGrid = Make1DGrid(areas.Ys, areas.YSteps, areas.YRelaxRatios, areas.IsReversedYRelaxRatio).ToArray();
        XGrid = xGrid;
        YGrid = yGrid;
        var xNumberPoints = xGrid.Length;
        var xNumberSegments = xNumberPoints - 1;
        var yNumberPoints = yGrid.Length;
        var yNumberSegments = yNumberPoints - 1;
        XNumberOfSegments = xNumberSegments;
        YNumberOfSegments = yNumberSegments;
        Elements = new Element[(xNumberPoints - 1) * (yNumberPoints - 1)];
        Points = new Vector2[xNumberPoints * yNumberPoints];
        for (var i = 0; i < yNumberPoints; i++)
        {
            for (var j = 0; j < xNumberPoints; j++)
            {
                Points[i * xNumberPoints + j] = new Vector2(xGrid[j], yGrid[i]);
            }
        }

        for (var i = 0; i < yNumberSegments; i++)
        {
            for (var j = 0; j < xNumberSegments; j++)
            {
                var elementPoints = new int[4]
                {
                    i * xNumberPoints + j,
                    i * xNumberPoints + j + 1,
                    (i + 1) * xNumberPoints + j,
                    (i + 1) * xNumberPoints + j + 1
                };
                Elements[i * xNumberSegments + j] = new Element(elementPoints);
            }
        }

        foreach (var element in Elements)
        {
            var xmin = Points[element.ElementNumbers[0]].X;
            var xmax = Points[element.ElementNumbers[1]].X;
            var ymin = Points[element.ElementNumbers[0]].Y;
            var ymax = Points[element.ElementNumbers[2]].Y;
            var xcenter = (xmax + xmin) / 2;
            var ycenter = (ymax + ymin) / 2;
            for (var j = 0; j < areas.Xmins.Length; j++)
            {
                if (xcenter >= areas.Xmins[j] && xcenter <= areas.Xmaxs[j] && ycenter >= areas.Ymins[j] &&
                    ycenter <= areas.Ymaxs[j])
                {
                    element.CurrentDensity = areas.CurrentDensitys[j];
                    element.MagneticPermeability = areas.MagneticPermeabilitys[j];
                }
            }
        }

        PointsBounds = new List<int>[xNumberPoints * yNumberPoints];
        for (var i = 0; i < PointsBounds.Length; i++)
        {
            PointsBounds[i] = new List<int>();
        }

        foreach (var element in Elements)
        {
            for (var i = 0; i < 4; i++)
            {
                for (var j = 0; j < 4; j++)
                {
                    if (element.ElementNumbers[j] < element.ElementNumbers[i] &&
                        !PointsBounds[element.ElementNumbers[i]].Contains(element.ElementNumbers[j]))
                    {
                        PointsBounds[element.ElementNumbers[i]].Add(element.ElementNumbers[j]);
                    }
                }
            }
        }

        foreach (var pointBounds in PointsBounds)
        {
            pointBounds.Sort();
        }
    }
}