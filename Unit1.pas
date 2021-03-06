unit Unit1;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, VclTee.TeeGDIPlus,
  VCLTee.TeEngine, VCLTee.Series, Vcl.Grids, Vcl.ExtCtrls, VCLTee.TeeProcs,
  VCLTee.Chart, Vcl.Imaging.jpeg, StrUtils;

const
  MaxPop = 1000; { ������������ ����� ��������� }
  LenChrome = 20; { ����� ����� �� ���� ���������� �������� }
  dim = 2;   { ����������� ������������ ������ }

type
  Allele = boolean;  {����� - ������� � ������� ������ }
  Chromosome = array [1..LenChrome * Dim] of Allele; { ������� ������ }
  Fenotype = array [1..Dim] of double;

  Individual = record
    Chrom: Chromosome; { ������� = ������� ������ }
    x: Fenotype;       { ������� = ������ ������������ ��������� ����� � ������������ ������ }
    Fitness: double;   { �������� ������� ������� }
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
    Chart2: TChart;
    Series51: TLineSeries;
    RadioButton1: TRadioButton;
    RadioButton2: TRadioButton;
    Image1: TImage;
    Image2: TImage;
    Panel3: TPanel;
    Label10: TLabel;
    Edit5: TEdit;
    Edit6: TEdit;
    Label11: TLabel;
    Label12: TLabel;
    RadioButton3: TRadioButton;
    RadioButton4: TRadioButton;
    Series52: TPointSeries;
    Label2: TLabel;
    Edit7: TEdit;
    Label13: TLabel;
    Edit8: TEdit;
    Label14: TLabel;
    Label15: TLabel;
    Label16: TLabel;
    procedure Button1Click(Sender: TObject);
    procedure FormActivate(Sender: TObject);
    procedure Image1Click(Sender: TObject);
    procedure Image2Click(Sender: TObject);
    procedure Edit1KeyPress(Sender: TObject; var Key: Char);
    procedure Edit5KeyPress(Sender: TObject; var Key: Char);
  private
    { Private declarations }
  public
    { Public declarations }
  end;


var
  Form1: TForm1;
  X: Fenotype;
  i,j,n, dimDynamic : integer;
  xMax: Fenotype;  {������ ������������ �������� ��� ��������� ����� � ������������ ������}
  xMin: Fenotype;  {������ ����������� �������� ��� ��������� ����� � ������������ ������}
  { ��� ���������������� ��������� - ������, ����� � ������������� }
  OldPop, NewPop, IntPop: Population;
  { ���������� ����� ����������}
  PopSize, Gen, h, s, b: integer;
  { �������� �������, ����������� � ���������� ��������� }
  NMutation, NCross, NGen: integer;
  { �������������� ���������� }
  Avg, Min, Max, SumFitness, xMinS, xMaxS, PMutation, PCross: double;
  Search:string;

implementation

{$R *.dfm}

function ObjFunc( x: Fenotype): real;   { ������� }
begin
  if dimDynamic = 1 then
    ObjFunc := 5-24*x[1]+17*x[1]*x[1]-(11/3)*x[1]*x[1]*x[1]+(1/4)*x[1]*x[1]*x[1]*x[1]
  else
    ObjFunc := x[1]*x[1]+x[2]*x[2];
end;


{ ������������� ������� - true ���� ���� }
function Flip(Probability: double): boolean;
begin
  Flip := Random <= Probability;
end;

{ �������������  ������ � ������ ������������ ��������� }
procedure Decode(Chrom: Chromosome; var x: fenotype);
var
  i, j, f, accum: longint;
begin
  for i := 1 to dimDynamic do
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
  {������������� }
  SumFitness := Pop[1].Fitness;
  Min := Pop[1].Fitness;
  Max := Pop[1].Fitness;
  xMinS:= pop[1].x[1];
  xMaxS:= pop[1].x[1];
  {���� ��� max, min, sumfitness }
  for j := 2 to PopSize do
    with Pop[j] do
    begin
      { ���������� ����� �������� ������� ����������� }
      SumFitness := SumFitness + Fitness;
      if Fitness > Max then begin
        Max := Fitness; { ����� �������� Max }
        xMaxS := x[1];
      end;
      if Fitness < Min then begin
        Min := Fitness; { ����� �������� Min }
        xMinS := x[1];
      end;
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
      for j := 1 to LenChrome * dimDynamic do Chrom[j] := Flip(0.5); { ������ ������ }
      Decode(Chrom, x); { ������������� ������ }
      Fitness := ObjFunc(x);{ ���������� ��������� �������� ������� ����������� }
    end;
end;


{3 ������������ ���������: ������ (select), ����������� (crossover) � ������� (mutation)}
{��������� ���������� ������}
procedure Select(Search:string);
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
  function Select1(Search:string): integer;
  var
    i, j, m: integer;
  begin
    if ipick > PopSize then
    begin
      Shuffle(OldPop);
      ipick := 1;
    end;
    {������� ��� �������� �������� ������� ��������� }
    i := ipick;
    j := ipick + 1;
    {� ������ �� ��� ���������� � ������������� ���������}
    if Search='min' then
    begin
      if OldPop[j].Fitness < OldPop[i].Fitness then m := j else m := i;
    end
    else
    begin
      if OldPop[j].Fitness > OldPop[i].Fitness then m := j else m := i;
    end;

    Inc(ipick, 2);
    Select1 := m;
  end;

begin
  ipick := 1;
  for i := 1 to PopSize do
    IntPop[i] := OldPop[Select1(Search)];
  OldPop := IntPop;
end;

{ �������� ������������ ������� }
{������� ������ ���� � ������ (������) � ������������ pmutation}
{ alleleval � ��� ��� �������}
function Mutation(alleleval: Allele; var NMutation: integer): Allele;
{alleleval � ��� ��� �������;
PMutation � ����������� �������;
NMutation � ������� �������;}
begin
  if Flip(PMutation) then { ������� � ������������ PMutation }
  begin
    Inc(NMutation);   { ���������� ������� ������� }
    Mutation := not alleleval; { ��������� ������� }
  end else begin
    Mutation := alleleval;  { �� ������ ������� }
  end;
end;

{ �������� ������������� ����������� }
{ ��������� ����������� 2 ������������ �����, ��������� ���������� � 2 �������-��������}
procedure Crossover(var Parent1, Parent2, Child1, Child2: Chromosome; flchrom: integer; var NCross, NMutation: integer);
{parent1, parent2 � ��������� ���������}
{child1,child2 � ��������� ��������}
{flchrom � ����� ��������� (���������� �����)}
{ncross, nmutation � �������� ���������� ����������� � �������}
{jcross � ����� �������.}
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
procedure Generation(Search:string);
var
  i: integer;
begin
  Select(Search);
  i := 1;
  repeat
    { ����������� �����, ����������� � ������� ���� ��������� ��
    ������������ ����� ��������� newpop }
    { ����������� � ������� - ������� ��������� � ��������� ����������� }
    Crossover(OldPop[i].Chrom, OldPop[i + 1].Chrom, NewPop[i].Chrom, NewPop[i + 1].Chrom, LenChrome * dimDynamic, NCross, NMutation);
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
procedure Pict(Chart1:TChart); // ���������� ������� ��� ���� ����������
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
        Chart1.Series[j].AddXY(X[1],ObjFunc(x));// ���������� X[1] � ObjFunc(x)
        X[1]:=X[1]+abs(xMin[1]-xMax[1])/200;
      end;                                     // while X[1]<=xMax[1]
      j:=j+1;                                  // ������� � ���������� Series[j]
      X[2]:=X[2]+abs(xMin[2]-xMax[2])/(M-1);   // � X[2]
  end;                                         // while X[2]<=xMax[2]
end;

//==============================================================================


  procedure plotting(Chart1:TChart); {���������� ������� ������� ��� ����� ����������}
  var i:fenotype;
    h:Real;
  begin
    Chart1.Series[0].Clear;
    h:= (xmax[1] - xmin[1])/100;
    i[1]:=xmin[1];
    while i[1]<=xmax[1] do begin
       Chart1.Series[0].AddXY(i[1],objfunc(i));
       i[1]:=i[1]+h;
    end;
    Chart1.Update;
  end;

//==============================================================================


procedure TForm1.Button1Click(Sender: TObject);


  procedure plottingDots(popsize:integer; var pop:population);  {����� �� ������ ������� ����������� }
  var j:integer;
  begin
    Chart2.Series[1].Clear;
    for j := 1 to popsize do with pop[j] do begin
      if Search='min' then
      begin
        if Min=fitness then
        begin
          Chart2.Series[1].AddXY(xMinS,min,'',clBlue);
        end
        else
        begin
          Chart2.Series[1].AddXY(x[1],fitness);
        end;
      end
      else
      begin
        if Max=fitness then
        begin
          Chart2.Series[1].AddXY(xMaxS,Max,'',clBlue);
        end
        else
        begin
          Chart2.Series[1].AddXY(x[1],fitness);
        end;
      end;
    end;
    Chart2.Update;
    Sleep(300);
  end;


var i,j:integer;
    Result,BestResult: double;
begin

  // ��������� ���� ����� �� �������
  for i:=0 to form1.ComponentCount-1 do
  begin
    if form1.Components[i].ClassNameIs('TEdit') then
    begin
      if tedit(form1.Components[i]).Text='' then
      begin
        MessageDlg('���� "' + tedit(form1.Components[i]).Hint + '" �� ������ ���� ������!',mtWarning,[mbOK],0);
        exit;
      end;
    end;
  end;

  Randomize; { ������������� ���������� ��������� ����� }

  NGen := StrToInt(Edit6.Text); {������ ��������� }
  PopSize := StrToInt(Edit5.Text);{ ���������� ��������� } 
  PMutation := StrToFloat(Edit7.Text)/100; { ����������� ������� }
  PCross := StrToFloat(Edit8.Text)/100;;   { ����������� ����������� }

  if (PopSize <= 0) or (PopSize >= 1000) then
  begin
    MessageDlg('���������� ��������� ������ ������ 0 � ������ 1000',mtWarning,[mbOK],0);
    exit;
  end;
  if NGen<=0 then
  begin
    MessageDlg('������ ��������� ������ ������ 0',mtWarning,[mbOK],0);
    exit;
  end;

  if (PMutation < 0) or (PMutation > 100) then
  begin
    MessageDlg('����������� ������� ������ ������ 0 � ������ 100',mtWarning,[mbOK],0);
    exit;
  end;
  
  if (PCross < 0) or (PCross > 100) then
  begin
    MessageDlg('����������� ����������� ������ ������ 0 � ������ 100',mtWarning,[mbOK],0);
    exit;
  end;

  if RadioButton3.Checked = true then Search:='min' else Search:='max';

  Result := 0; { ������������� ���������� ������ }

  if dimDynamic = 2 then
  begin
    Chart1.Visible:=true;
    Chart2.Visible:=false;
    xMax[1]:= StrToFloat(Edit2.Text);
    xMax[2]:= StrToFloat(Edit4.Text);
    xMin[1]:= StrToFloat(Edit1.Text);
    xMin[2]:= StrToFloat(Edit3.Text);
    if (xMin[1] > xMax[1]) or (xMin[2] > xMax[2]) then begin
      MessageDlg('�������� ������ �����������',mtWarning,[mbOK],0);
      exit;
    end;
    Pict(Chart1);
  end
  else
  begin
    Chart1.Visible:=false;
    Chart2.Visible:=true;
    xMin[1]:= StrToFloat(Edit1.Text);
    xMax[1]:= StrToFloat(Edit2.Text);
    if (xMin[1] > xMax[1]) then begin
      MessageDlg('�������� ������ �����������',mtWarning,[mbOK],0);
        exit;
    end;
    plotting(Chart2);
  end;

  Button1.Enabled:=false;
  Button1.Update;
  NMutation := 0;  { ������������� �������� ������� }
  NCross := 0; { ������������� �������� ����������� }
  InitPop; { �������� ��������� ��������� }
  Statistics(Max, Avg, Min, OldPop);
  if dimDynamic=1 then plottingDots(PopSize,OldPop);   {����� �� ������ ������� ����������� ���������}
  if Search='min' then begin
    BestResult:= Min
  end else begin
    BestResult:= Max;
  end;
  Gen := 1;   { ��������� �������� ��������� � 0 }
  repeat { ������� ������������ ���� }
    Generation(Search);
    Statistics(Max, Avg, Min, NewPop);
    if Search='min' then  begin
      if Min < BestResult then BestResult := Min
    end
    else begin
      if Max > BestResult then BestResult := Max;
    end;
    OldPop := NewPop;
    if dimDynamic=1 then plottingDots(PopSize,OldPop);   {����� �� ������ ������� ����������� ���������}
    {������� �� ����� ��������� }
    Inc(Gen);
  until Gen > PopSize;

  if Search='min' then
  begin
    Label1.Caption:='�������: '+FloatToStrF(BestResult,ffFixed,10,4);
  end
  else
  begin
    Label1.Caption:='��������: '+FloatToStrF(BestResult,ffFixed,10,4);
  end;
  Application.ProcessMessages;
  Button1.Enabled:=true;
end;

//==============================================================================
procedure TForm1.Edit1KeyPress(Sender: TObject; var Key: Char);
begin
  if not (Key in ['0'..'9', #8,'-'])then Key:=#0;
end;

procedure TForm1.Edit5KeyPress(Sender: TObject; var Key: Char);
begin
  if not (Key in ['0'..'9', #8])then Key:=#0;
end;

procedure TForm1.FormActivate(Sender: TObject);
var
  i,m : integer;
begin
  dimDynamic:=2;
  M:=Chart1.CountActiveSeries; // ���������� ����� �������� ������ Series
  for i := 0 to M-1 do
    Chart1.Series[i].Clear;
end;

procedure TForm1.Image1Click(Sender: TObject);
begin
RadioButton1.Checked:=true;
Panel2.Visible:=true;
Panel3.Top:=79;
dimDynamic:=2;
Edit1.Text:=IntToStr(-89);
Edit2.Text:=IntToStr(89);
end;

procedure TForm1.Image2Click(Sender: TObject);
begin
RadioButton2.Checked:=true;
Panel2.Visible:=false;
Panel3.Top:=41;
dimDynamic:=1;
Edit1.Text:=IntToStr(0);
Edit2.Text:=IntToStr(8);
end;

//=============================================================================
end.
