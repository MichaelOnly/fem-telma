namespace fem_telma;

public class Areas
{
    public float[] Xmins { get; set; }
    public float[] Xmaxs { get; set; }
    public float[] Ymins { get; set; }
    public float[] Ymaxs { get; set; }
    public double[] CurrentDensitys { get; set; }
    public double[] MagneticPermeabilitys { get; set; }
    public float[] Xs { get; set; }
    public float[] XSteps { get; set; }
    public float[] XRelaxRatios { get; set; }
    public bool[] IsReversedXRelaxRatio { get; set; }
    public float[] Ys { get; set; }
    public float[] YSteps { get; set; }
    public float[] YRelaxRatios { get; set; }
    public bool[] IsReversedYRelaxRatio { get; set; }
}