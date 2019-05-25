unit Unit1;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, VclTee.TeeGDIPlus,
  VCLTee.TeEngine, VCLTee.Series, Vcl.Grids, Vcl.ExtCtrls, VCLTee.TeeProcs,
  VCLTee.Chart, Vcl.Imaging.jpeg;

const
  MaxPop = 1000; { ������������ ����� ��������� }
  LenChrome = 20; { ����� ����� �� ���� ���������� �������� }
  dim = 2;   { ����������� ������������ ������ }
  PMutation = 0.01; { ����������� ������� }
  PCross = 0.9;   { ����������� ����������� }
  NN = 30; {����� ��������}

type
  Allele = boolean;  {����� - ������� � ������� ������ }
  Chromosome = array [1..LenChrome * Dim] of Allele; { ������� ������ }
  Fenotype = array [1..Dim] of double;

  Individual = record
    Chrom: Chromosome; { ������� = ������� ������ }
    x: Fenotype;
    { ������� = ������ ������������ ��������� ����� � ������������ ������ }
    Fitness: double; { �������� ������� ������� }
  end;
  Population = array [1..maxpop] of Individual;
  TForm1 = class(TForm)
    Button1: TButton;
    Chart1: TChart;
    Series15: TFastLineSeries;
    Series1: TFastLineSeries;
    Series2: TFastLineSeries;
    Series3: TFastLineSeries;
    Series4: TFastLineSeries;
    Series5: TFastLineSeries;
    Series6: TFastLineSeries;
    Series7: TFastLineSeries;
    Series8: TFastLineSeries;
    Series9: TFastLineSeries;
    Series10: TFastLineSeries;
    Series11: TFastLineSeries;
    Series12: TFastLineSeries;
    Series13: TFastLineSeries;
    Series14: TFastLineSeries;
    Series16: TFastLineSeries;
    Series17: TFastLineSeries;
    Series18: TFastLineSeries;
    Series19: TFastLineSeries;
    Series20: TFastLineSeries;
    Series21: TFastLineSeries;
    Series22: TFastLineSeries;
    Series23: TFastLineSeries;
    Series24: TFastLineSeries;
    Series25: TFastLineSeries;
    Series26: TFastLineSeries;
    Series27: TFastLineSeries;
    Series28: TFastLineSeries;
    Series29: TFastLineSeries;
    Series30: TFastLineSeries;
    Series31: TFastLineSeries;
    Series32: TFastLineSeries;
    Series33: TFastLineSeries;
    Series34: TFastLineSeries;
    Series35: TFastLineSeries;
    Series36: TFastLineSeries;
    Series37: TFastLineSeries;
    Series38: TFastLineSeries;
    Series39: TFastLineSeries;
    Series40: TFastLineSeries;
    Series41: TFastLineSeries;
    Series42: TFastLineSeries;
    Series43: TFastLineSeries;
    Series44: TFastLineSeries;
    Series45: TFastLineSeries;
    Series46: TFastLineSeries;
    Series47: TFastLineSeries;
    Series48: TFastLineSeries;
    Series49: TFastLineSeries;
    Series50: TFastLineSeries;
    Label1: TLabel;
    Label2: TLabel;
    GroupBox1: TGroupBox;
    Label3: TLabel;
    Edit1: TEdit;
    Edit2: TEdit;
    Panel1: TPanel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Panel2: TPanel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    Edit3: TEdit;
    Edit4: TEdit;
    Label10: TLabel;
    Edit5: TEdit;
    Label11: TLabel;
    Edit6: TEdit;
    procedure Button1Click(Sender: TObject);
    procedure FormActivate(Sender: TObject);



  private
    { Private declarations }
  public
    { Public declarations }
  end;


var
  Form1: TForm1;
  X: Fenotype;
  i,j,n : integer;
  xMax: Fenotype;  {������ ������������ �������� ��� ��������� ����� � ������������ ������}
  xMin: Fenotype;  {������ ����������� �������� ��� ��������� ����� � ������������ ������}
  { ��� ���������������� ��������� - ������, ����� � ������������� }
  OldPop, NewPop, IntPop: Population;
  { ���������� ����� ����������}
  PopSize, Gen, h, s, b: integer;
  { �������� �������, ����������� � ���������� ��������� }
  NMutation, NCross, NGen: integer;
  { �������������� ���������� }
  Avg, Min, Max, BestMin, BestMax, result, SumFitness: double;



implementation

{$R *.dfm}

function ObjFunc( x: Fenotype): real;
begin
  //ObjFunc := 15.5+sqr(2.3-x[1])+sqr(4.1-x[2]);
  ObjFunc := x[1]*x[1]+x[2]*x[2];
end;


{ ������������� ������� - true ���� ���� }
function Flip(Probability: double): boolean;
begin
  Flip := Random <= Probability;
end;

{ �������������  ������ � ������ ������������ ��������� � ������������ ������ }
procedure Decode(Chrom: Chromosome; var x: fenotype);
var
  i, j, f, accum: longint;
begin
  for i := 1 to Dim do
  begin
    Accum := 0;
    f := 1;
    for j := 1 + LenChrome * (i - 1) to LenChrome + LenChrome * (i - 1) do
    begin
      if Chrom[j] then
        Inc(Accum, f);
      f := f * 2;
    end;
    x[i] := xmin[i] + (xmax[i] - xmin[i]) * Accum / (f - 1);
  end;
end;

{ ������ �������������� ������� }
procedure Statistics(var Max, Avg, Min: double; Pop: Population);
var
  j: integer;
  SumFitness: double;
begin
  SumFitness := Pop[1].Fitness;
  Min := Pop[1].Fitness;
  Max := Pop[1].Fitness;
  for j := 2 to PopSize do
    with Pop[j] do
    begin
      { ���������� ����� �������� ������� ����������� }
      SumFitness := SumFitness + Fitness;
      if Fitness > Max then
        Max := Fitness; { ����� �������� Max }
      if Fitness < Min then
        Min := Fitness; { ����� �������� Min }
    end;
  Avg := SumFitness / PopSize;   { ������ �������� }
end;

{ ������������� ��������� ��������� ��������� ������� }
procedure InitPop;
var
  i, j: integer;
begin
  for i := 1 to PopSize do
    with OldPop[i] do
    begin
      for j := 1 to LenChrome * Dim do
        Chrom[j] := Flip(0.5); { ������ ������ }
      Decode(Chrom, x); { ������������� ������ }
      { ���������� ��������� �������� ������� ����������� }
      Fitness := ObjFunc(x);
    end;
end;

{ �������� ������ (��������) }
procedure Select;
var
  ipick, i: integer;

  { ��������� ������������� ��������� � �������� ������ }
  procedure Shuffle(var pop: Population);
  var
    i, j: integer;
    ind0: Individual;
  begin
    for i := 1 to PopSize do
    begin
      j := 1 + Random(i);
      { ������������ }
      ind0 := pop[i];
      pop[i] := pop[j];
      pop[j] := ind0;
    end;
  end;

  { ����� ��������� ������ ��� ��������� ��� �������� � ��������� ��������� }
  function Select1: integer;
  var
    i, j, m: integer;
  begin
    if ipick > PopSize then
    begin
      Shuffle(OldPop);
      ipick := 1;
    end;
    i := ipick;
    j := ipick + 1;
    if OldPop[j].Fitness < OldPop[i].Fitness then
      m := j
    else
      m := i;
    Inc(ipick, 2);
    Select1 := m;
  end;

begin
  ipick := 1;
  for i := 1 to PopSize do
    IntPop[i] := OldPop[Select1];
  OldPop := IntPop;
end;

{ �������� ������������ ������� }
function Mutation(alleleval: Allele; var NMutation: integer): Allele;
begin
  Mutation := alleleval;
  if Flip(PMutation) then { ������� � ������������ PMutation }
  begin
    Inc(NMutation);   { ���������� ������� ������� }
    Mutation := not alleleval; { ��������� ������� }
  end;
end;

{ �������� ������������� ����������� }
procedure Crossover(var Parent1, Parent2, Child1, Child2: Chromosome; flchrom: integer;
var NCross, NMutation: integer);
var
  i, jcross: integer;
begin
  if Flip(PCross) then { ����������� ����������� � ������������ PCross }
  begin
    { ����������� ����� ������� � ��������� ����� 1 � flchrom-1 }
    jcross := 1 + Random(flchrom);
    Inc(NCross); { ����������� �������� ����������� }
    { ������ ����� ������, 1 � 1 � 2 � 2 }
    { ���������� ����� �� ����� ������� }
    for i := 1 to jcross do
    begin
      { ������ � �������� � ������������ pmutation }
      { ������ ������� }
      Child1[i] := Mutation(Parent1[i], NMutation);
      { ������ ������� }
      Child2[i] := Mutation(Parent2[i], NMutation);
    end;
    { ������ ����� ������, 1 � 2 � 2 � 1 }
    { ���������� ����� ����� ����� ������� }
    for i := jcross + 1 to flchrom do
    begin
      { ������ � �������� � ������������ pmutation }
      { ������ ������� }
      Child1[i] := Mutation(Parent2[i], NMutation);
      { ������ ������� }
      Child2[i] := Mutation(Parent1[i], NMutation);
    end;
  end;
end;

{ ������������� ������ ��������� ��� ������ ������, ����������� � ������� }
{ ��������������, ��� ��������� ����� ������ ������ }
procedure Generation;
var
  i: integer;
begin
  Select;
  i := 1;
  repeat
  { ����������� �����, ����������� � ������� ���� ��������� ��
  ������������ ����� ��������� newpop }
    { ����������� � ������� - ������� ��������� � ��������� ����������� }
    Crossover(OldPop[i].Chrom, OldPop[i + 1].Chrom,
      NewPop[i].Chrom, NewPop[i + 1].Chrom,
      LenChrome * Dim, NCross, NMutation);
    { ������������� ������ � ���������� ����������� }
    with NewPop[i] do
    begin
      Decode(Chrom, x);
      Fitness := ObjFunc(x);
    end;
    with Newpop[i + 1] do
    begin
      Decode(Chrom, x);
      Fitness := ObjFunc(x);
    end;
    Inc(i, 2);
  until i > PopSize;
end;




//==============================================================================
procedure Pict(Chart1:TChart);
  var
    i,j,M : integer;
begin
  M:=Chart1.CountActiveSeries;         // ������, ������� � ��� �������� ������
  for i := 0 to M-1 do Chart1.Series[i].Clear;   // � ������� ��

  X[2]:=xMin[2];
  j:=0;
  while X[2]<=xMax[2] do
  begin
      X[1]:=xMin[1];
      while X[1]<=xMax[1] do              {��� ������������� X[2]
                                                ��������� Series[j]}
      begin
        Chart1.Series[j].AddXY(X[1],ObjFunc(x));     // ���������� X[1] � ObjFunc(x)
        X[1]:=X[1]+abs(xMin[1]-xMax[1])/200;
      end;                                     // while X[1]<=xMax[1]
      j:=j+1;                                  // ������� � ���������� Series[j]
      X[2]:=X[2]+abs(xMin[2]-xMax[2])/(M-1);       // � X[2]
  end;                                         // while X[2]<=xMax[2]
end;
//==============================================================================
procedure TForm1.Button1Click(Sender: TObject);
var i,j:integer;
    RezultMin, RezultMax :Double;

begin

  Randomize; { ������������� ���������� ��������� ����� }
  NGen := StrToInt(Edit5.Text); { ���������� ��������� }
  PopSize := StrToInt(Edit6.Text); {������ ��������� }
  result := 0; { ������������� ���������� ������ }


  xMax[1]:= StrToFloat(Edit2.Text);
  xMax[2]:= StrToFloat(Edit4.Text);
  xMin[1]:= StrToFloat(Edit1.Text);
  xMin[2]:= StrToFloat(Edit3.Text);
  Pict(Chart1);

  NMutation := 0;  { ������������� �������� ������� }
  NCross := 0; { ������������� �������� ����������� }
  InitPop; { �������� ��������� ��������� }
  Statistics(Max, Avg, Min, OldPop);
  BestMin := Min;
  BestMax := Max;
  Gen := 1;   { ��������� �������� ��������� � 0 }

  for b := 1 to NN do { ����������� N ��� ��� ��������� ������������� }
  begin
      repeat { ������� ������������ ���� }
        Generation;
        Statistics(Max, Avg, Min, NewPop);
        if Min < BestMin then BestMin := Min;
        if Max > BestMax then BestMax := Max;
        OldPop := NewPop;
        {������� �� ����� ��������� }
        Inc(Gen);
      until Gen > PopSize;
      RezultMin := RezultMin + BestMin;
      RezultMax := RezultMax + BestMax;
  end;
  RezultMin:= RezultMin / NN;
  RezultMax:= RezultMax / NN;
  Label1.Caption:='�������'+FloatToStrF(RezultMin,ffFixed,10,4);
  Label2.Caption:='��������'+FloatToStrF(RezultMax,ffFixed,10,4);
end;

//==============================================================================
procedure TForm1.FormActivate(Sender: TObject);
var
  i,m : integer;
begin

  M:=Chart1.CountActiveSeries; // ���������� ����� �������� ������ Series
  for i := 0 to M-1 do
    Chart1.Series[i].Clear;
end;

//=============================================================================
end.
