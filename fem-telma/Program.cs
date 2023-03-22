﻿// See https://aka.ms/new-console-template for more information
using System.Text.Json;
using fem_telma;
using fem_telma.Grid;
using fem_telma.SLAE;

var area = JsonSerializer.Deserialize<Areas>(File.ReadAllText("areas.json"));
var points = JsonSerializer.Deserialize<Points>(File.ReadAllText("points.json"));
var grid = new Grid(area);
var SLAE = new SLAE(grid);
SLAE.SolveWithLOSPrecond(1e-15,1000);
var resFunc = new ResultFunction(SLAE.Result,grid);
for (var i = 0; i < points.Xs.Length; i++)
{
    Console.WriteLine("x: "+ points.Xs[i].ToString() + " y: " + points.Ys[i].ToString());
    Console.WriteLine("Az: "+resFunc.GerResultFunction(points.Xs[i], points.Ys[i]).ToString());
}