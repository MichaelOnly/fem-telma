// See https://aka.ms/new-console-template for more information
using System.Text.Json;
using fem_telma;
using fem_telma.Grid;

var area = JsonSerializer.Deserialize<Areas>(File.ReadAllText("areas.json"));
var grid = new Grid(area);
