namespace fem_telma;

public class Element
{
    public int[] ElementNumbers { get; init; }
    public double CurrentDensity { get; set; }
    public double MagneticPermeability { get; set; }
    public Element(int[] elementNumbers)
    {
        ElementNumbers = elementNumbers;
    }
}