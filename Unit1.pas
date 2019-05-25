unit Unit1;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, VclTee.TeeGDIPlus,
  VCLTee.TeEngine, VCLTee.Series, Vcl.Grids, Vcl.ExtCtrls, VCLTee.TeeProcs,
  VCLTee.Chart, Vcl.Imaging.jpeg;

const
  MaxPop = 1000; { Максимальное число поколений }
  LenChrome = 20; { Число битов на один кодируемый параметр }
  dim = 2;   { Размерность пространства поиска }
  PMutation = 0.01; { Вероятность мутации }
  PCross = 0.9;   { Вероятность скрещивания }
  NN = 10; {число прогонов}

type
  Allele = boolean;  {Алель - позиция в битовой строке }
  Chromosome = array [1..LenChrome * Dim] of Allele; { Битовая строка }
  Fenotype = array [1..Dim] of double;

  Individual = record
    Chrom: Chromosome; { Генотип = битовая строка }
    x: Fenotype;       { Фенотип = массив вещественных координат точки в пространстве поиска }
    Fitness: double;   { Значение целевой функции }
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
    procedure Button1Click(Sender: TObject);
    procedure FormActivate(Sender: TObject);
    procedure Image1Click(Sender: TObject);
    procedure Image2Click(Sender: TObject);



  private
    { Private declarations }
  public
    { Public declarations }
  end;


var
  Form1: TForm1;
  X: Fenotype;
  i,j,n, dimDynamic : integer;
  xMax: Fenotype;  {массив максимальных значений для координат точки в пространстве поиска}
  xMin: Fenotype;  {массив минимальных значений для координат точки в пространстве поиска}
  { Три непересекающихся популяции - старая, новая и промежуточная }
  OldPop, NewPop, IntPop: Population;
  { Глобальные целые переменные}
  PopSize, Gen, h, s, b: integer;
  { Счетчики мутаций, скрещиваний и количество поколений }
  NMutation, NCross, NGen: integer;
  { Статистические переменные }
  Avg, Min, Max, SumFitness: double;
  Search:string;



implementation

{$R *.dfm}

function ObjFunc( x: Fenotype): real;
begin

  if dimDynamic = 1 then
    ObjFunc := 5-24*x[1]+17*x[1]*x[1]-(11/3)*x[1]*x[1]*x[1]+(1/4)*x[1]*x[1]*x[1]*x[1]
  else
    ObjFunc := x[1]*x[1]+x[2]*x[2];
end;


{ Подбрасывание монетки - true если орел }
function Flip(Probability: double): boolean;
begin
  Flip := Random <= Probability;
end;

{ Декодирование  строки в массив вещественных координат }
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

{ Расчет статистических величин }
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
      { Накопление суммы значений функции пригодности }
      SumFitness := SumFitness + Fitness;
      if Fitness > Max then
        Max := Fitness; { Новое значение Max }
      if Fitness < Min then
        Min := Fitness; { Новое значение Min }
    end;
  Avg := SumFitness / PopSize;   { Расчет среднего }
end;

{ Инициализация начальной популяции случайным образом }
procedure InitPop;
var
  i, j: integer;
begin
  for i := 1 to PopSize do
    with OldPop[i] do
    begin
      for j := 1 to LenChrome * dimDynamic do
        Chrom[j] := Flip(0.5); { Бросок монеты }
      Decode(Chrom, x); { Декодирование строки }
      { Вычисление начальных значений функции пригодности }
      Fitness := ObjFunc(x);
    end;
end;

{ Оператор отбора (селекции) }
procedure Select(Search:string);
var
  ipick, i: integer;

  { Процедура перемешивания популяции в процессе отбора }
  procedure Shuffle(var pop: Population);
  var
    i, j: integer;
    ind0: Individual;
  begin
    for i := 1 to PopSize do
    begin
      j := 1 + Random(i);
      { Перемешиваем }
      ind0 := pop[i];
      pop[i] := pop[j];
      pop[j] := ind0;
    end;
  end;

  { Отбор наилучших особей для популяции для перехода в следующее поколение }
  function Select1(Search:string): integer;
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

{ Оператор инверсионной мутации }
function Mutation(alleleval: Allele; var NMutation: integer): Allele;
begin
  Mutation := alleleval;
  if Flip(PMutation) then { Мутация с вероятностью PMutation }
  begin
    Inc(NMutation);   { Наращиваем счетчик мутаций }
    Mutation := not alleleval; { Совершаем мутацию }
  end;
end;

{ Оператор одноточечного скрещивания }
procedure Crossover(var Parent1, Parent2, Child1, Child2: Chromosome; flchrom: integer;
var NCross, NMutation: integer);
var
  i, jcross: integer;
begin
  if Flip(PCross) then { Выполняется скрещивание с вероятностью PCross }
  begin
    { Определение точки сечения в диапазоне между 1 и flchrom-1 }
    jcross := 1 + Random(flchrom);
    Inc(NCross); { Наращивание счетчика скрещиваний }
    { Первая часть обмена, 1 в 1 и 2 в 2 }
    { Обмениваем часть до точки сечения }
    for i := 1 to jcross do
    begin
      { Заодно и мутируем с вероятностью pmutation }
      { Первый потомок }
      Child1[i] := Mutation(Parent1[i], NMutation);
      { Второй потомок }
      Child2[i] := Mutation(Parent2[i], NMutation);
    end;
    { Вторая часть обмена, 1 в 2 и 2 в 1 }
    { Обмениваем часть после точки сечения }
    for i := jcross + 1 to flchrom do
    begin
      { Заодно и мутируем с вероятностью pmutation }
      { Первый потомок }
      Child1[i] := Mutation(Parent2[i], NMutation);
      { Второй потомок }
      Child2[i] := Mutation(Parent1[i], NMutation);
    end;
  end;
end;

{ Генерирование нового поколения при помощи отбора, скрещивания и мутации }
{ Предполагается, что популяция имеет четный размер }
procedure Generation(Search:string);
var
  i: integer;
begin
  Select(Search);
  i := 1;
  repeat
  { Выполняются отбор, скрещивание и мутация пока полностью не
  сформируется новая популяция newpop }
    { Скрещивание и мутация - мутация вставлена в процедуру скрещивания }
    Crossover(OldPop[i].Chrom, OldPop[i + 1].Chrom,
      NewPop[i].Chrom, NewPop[i + 1].Chrom,
      LenChrome * dimDynamic, NCross, NMutation);
    { Декодирование строки и вычисление пригодности }
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
procedure Pict(Chart1:TChart); // построение графика для двух переменных
  var
    i,j,M : integer;
begin
  M:=Chart1.CountActiveSeries;         // Узнаем, сколько у нас хранилищ данных
  for i := 0 to M-1 do Chart1.Series[i].Clear;   // и очищаем их

  X[2]:=xMin[2];
  j:=0;
  while X[2]<=xMax[2] do
  begin
      X[1]:=xMin[1];
      while X[1]<=xMax[1] do              {при фиксированном X[2]
                                                заполняем Series[j]}
      begin
        Chart1.Series[j].AddXY(X[1],ObjFunc(x));     // значениями X[1] и ObjFunc(x)
        X[1]:=X[1]+abs(xMin[1]-xMax[1])/200;
      end;                                     // while X[1]<=xMax[1]
      j:=j+1;                                  // переход к очередному Series[j]
      X[2]:=X[2]+abs(xMin[2]-xMax[2])/(M-1);       // и X[2]
  end;                                         // while X[2]<=xMax[2]
end;

//==============================================================================


procedure plotting(Chart1:TChart); {построение графика функции для одной переменной}
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
var i,j:integer;
    Result,BestResult: double;
begin
  Application.ProcessMessages;
  Button1.Enabled:=false;
  Randomize; { Инициализация генератора случайных чисел }
  NGen := StrToInt(Edit5.Text); { Количество поколений }
  PopSize := StrToInt(Edit6.Text); {Размер популяции }
  if RadioButton3.Checked = true then Search:='min' else Search:='max';

  Result := 0; { Инициализация переменной ответа }

  if dimDynamic = 2 then
  begin
    Chart1.Visible:=true;
    Chart2.Visible:=false;
    xMax[1]:= StrToFloat(Edit2.Text);
    xMax[2]:= StrToFloat(Edit4.Text);
    xMin[1]:= StrToFloat(Edit1.Text);
    xMin[2]:= StrToFloat(Edit3.Text);
    Pict(Chart1);
  end
  else
  begin
    Chart1.Visible:=false;
    Chart2.Visible:=true;
    xMin[1]:= StrToFloat(Edit1.Text);
    xMax[1]:= StrToFloat(Edit2.Text);
    plotting(Chart2);
  end;

  for b := 1 to NN do { Прогоняется N раз для повышения достоверности }
  begin
      NMutation := 0;  { Инициализация счетчика мутация }
      NCross := 0; { Инициализация счетчика скрещиваний }
      InitPop; { Создание начальной популяции }
      Statistics(Max, Avg, Min, OldPop);

      if Search='min' then begin
        BestResult:= Min
      end else begin
        BestResult:= Max;
      end;
      Gen := 1;   { Установка счетчика поколений в 0 }
      repeat { Главный итерационный цикл }
        Generation(Search);
        Statistics(Max, Avg, Min, NewPop);
        if Search='min' then  begin
          if Min < BestResult then BestResult := Min
        end
        else begin
          if Max > BestResult then BestResult := Max;
        end;
        OldPop := NewPop;
        {переход на новое поколение }
        Inc(Gen);
      until Gen > PopSize;
      Result := Result + BestResult;
  end;
  Result:= Result / NN;
  if Search='min' then
  begin
    Label1.Caption:='Минимум: '+FloatToStrF(Result,ffFixed,10,4);
  end
  else
  begin
    Label1.Caption:='Максимум: '+FloatToStrF(Result,ffFixed,10,4);
  end;
  Button1.Enabled:=true;
end;

//==============================================================================
procedure TForm1.FormActivate(Sender: TObject);
var
  i,m : integer;
begin
  dimDynamic:=2;
  M:=Chart1.CountActiveSeries; // определяем число хранилищ данных Series
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
