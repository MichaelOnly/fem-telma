namespace fem_telma;

public class MuSpline
{
    private double[] _muTable;
    private double[] _BTable;

    public MuSpline()
    {
        using var reader = new StreamReader("mu.002");
        var countStrings = Int32.Parse(reader.ReadLine());
        _muTable = new double[countStrings];
        _BTable = new double[countStrings];
        for (var i = 0; i < countStrings; i++)
        {
            var tableString = reader.ReadLine().Split(' ');
            _muTable[i] = Double.Parse(tableString[1]);
            _BTable[i] = Double.Parse(tableString[2]);
        }
    }

    public double GetMu(double b)
    {
        if (b >= _BTable[^1])
        {
            return (_BTable[^1] * (1.0 / _muTable[^1] - 1.0) / b + 1.0) / (4 * Math.PI * 1e-7);
        }

        var numberSegment = 0;
        while (b >= _BTable[numberSegment] && b <= _BTable[numberSegment + 1])
        {
            numberSegment++;
        }

        var mu0 = _muTable[numberSegment];
        var mu1 = _muTable[numberSegment + 1];
        var b0 = _BTable[numberSegment];
        var b1 = _BTable[numberSegment + 1];
        return 1.0 / ((mu0 * ((b1 - b) / (b1 - b0)) + mu1 * ((b - b0) / (b1 - b0))) * 4 * Math.PI * 1e-7);
    }
}