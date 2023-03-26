// See https://aka.ms/new-console-template for more information
using System.Text.Json;
using fem_telma;
using fem_telma.Grid;
using fem_telma.SLAE;

var area = JsonSerializer.Deserialize<Areas>(File.ReadAllText("areas.json"));
var points = JsonSerializer.Deserialize<Points>(File.ReadAllText("points.json"));
var grid = new Grid(area);
var muSpline = new MuSpline();
var SLAE = new SLAE(grid);
SLAE.SolveWithLOSPrecond(1e-15,1000);
//SLAE.SolveWithSimpleIteration(grid, 1e-15,500, 0.2);
var resFunc = new ResultFunction(SLAE.Result,grid);
for (var i = 0; i < points.Xs.Length; i++)
{
    Console.WriteLine("x: "+ points.Xs[i] + " y: " + points.Ys[i]);
    Console.WriteLine("Az: "+resFunc.GerResultFunction(points.Xs[i], points.Ys[i]));
    Console.WriteLine("|B|: " + resFunc.GetAbsB(points.Xs[i],points.Ys[i]));
}