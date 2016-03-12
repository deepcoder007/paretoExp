program EM_complet_generare;
const nr=50;
type vect=array[1..4000] of real;
var i,n:integer;
    s1,s2,mu1,mu2,a1,a2,a0,ab0,ab1,ab2:real; {iterative values}
    s1r,s2r,mu1r,mu2r,a1r,a2r,a0r:real; {real values}
    ai:real;ae,mse:array[1..7] of real; {average estimate,mean square error,average no iter}
    x1,x2:vect;

Procedure generare;{generates n bivariate values}
var i:integer;
    y0,y1,y2,yt:real;

 function pa(s,mu,a:real):real;
  var u:real;
  begin
   u:=exp(-ln(random)/a);
   pa:=mu+s*(u-1);
   Writeln(pa);
 end; {pa}

begin
 mu1:=1000;mu2:=1000;
 randomize;
 for i:=1 to n do begin
   y0:=pa(1,0,a0r);y1:=pa(s1r,mu1r,a1r);y2:=pa(s2r,mu2r,a2r);
   yt:=s1r*y0+mu1r;
   if yt<y1 then x1[i]:=yt else x1[i]:=y1;
   yt:=s2r*y0+mu2r;
   if yt<y2 then x2[i]:=yt else x2[i]:=y2;
   if x1[i]<mu1 then mu1:=x1[i];
   if x2[i]<mu2 then mu2:=x2[i];
 end;
end; {generare}

Procedure iterPa(var x:vect; var s,a,mu:real); {M-step; calculates sigma iteratively}
  var sb,suma:real;
      i:integer;
 begin
  sb:=s; {sb va fi valoarea noua}
  repeat
   s:=sb;
   suma:=0;
   for i:=1 to n do suma:=suma+1/(s+x[i]-mu);
   sb:=n*(a0+a)/(a0+a+1)/suma;
  until abs(s-sb)<0.0001;
  s:=sb; {from last iteration}
end; {iterPa}

Procedure EM1pas; {E-step}
 var z1,z2,y0,y1,y2,sum0,sum1,sum2:real;
     a0n,a1n,a2n:real;neg,n1,n2:real;
     i:integer;
 begin
  sum0:=0;sum1:=0;sum2:=0;
  n1:=0;n2:=0;
  for i:=1 to n do begin
   z1:=(x1[i]-mu1)/s1;z2:=(x2[i]-mu2)/s2;
   sum1:=sum1+ln(1+z1);
   sum2:=sum2+ln(1+z2);
   if z1>z2 then sum0:=sum0+ln(1+z1) else sum0:=sum0+ln(1+z2);
  end; {for}
  neg:=n*a0/(a0+a1+a2);n1:=n*a2/(a0+a1+a2);n2:=n*a1/(a0+a1+a2);
  a0n:=n/(sum0+n2*a2/(a0+a2)/a0+n1*a1/(a0+a1)/a0);
  a1n:=n/(sum1+a0*n1/(a0+a1)/a1+neg/a1);
  a2n:=n/(sum2+a0*n2/(a0+a2)/a2+neg/a2);
  a0:=a0n;a2:=a2n;a1:=a1n;
end; {EM1pas}

begin {PP}
 write('n=');readln(n);
 mu1r:=1;s1r:=1;a1r:=2; {real values}
 mu2r:=1;s2r:=1;a2r:=2;
 a0r:=1;
 s1:=1.5; a1:=1.5; {initial values}
 s2:=1.5; a2:=1.5;
 a0:=1.5;
 for i:=1 to 7 do begin ae[i]:=0;mse[i]:=0 end;
 ai:=0;
 for i:=1 to nr do begin {repeats generare+em}
  generare; {returns mu1,mu2};
  {write(mu1:10:8,' ',mu2:10:8);readln;}
  repeat
   ab0:=a0;ab1:=a1;ab2:=a2;
   EM1pas; {E-step; modifica a-urile}
   iterPa(x1,s1,a1,mu1); {M-step; new s1}
   iterPa(x2,s2,a2,mu2); {M-step; new s2}
   ai:=ai+1;
  until (abs(a0-ab0)<0.000001) and (abs(a1-ab1)<0.000001) and (abs(a2-ab2)<0.000001);
  if (i mod 5)=0 then writeln(i,' a0=',a0:5:3);
  ae[1]:=ae[1]+s1;ae[2]:=ae[2]+s2;ae[3]:=ae[3]+a1;ae[4]:=ae[4]+a2;ae[5]:=ae[5]+a0;
  ae[6]:=ae[6]+mu1;ae[7]:=ae[7]+mu2;
  mse[1]:=mse[1]+sqr(s1-s1r);mse[2]:=mse[2]+sqr(s2-s2r);
  mse[3]:=mse[3]+sqr(a1-a1r);mse[4]:=mse[4]+sqr(a2-a2r);mse[5]:=mse[5]+sqr(a0-a0r);
  mse[6]:=mse[6]+sqr(mu1-mu1r);mse[7]:=mse[7]+sqr(mu2-mu2r);
 end; {for i}
 for i:=1 to 7 do begin
  ae[i]:=ae[i]/nr;mse[i]:=mse[i]/nr;
  writeln(i,' AE=',ae[i]:10:8,' MSE=',mse[i]:10:8);
 end;
 write('AI=',ai/nr:10:4);
 readln;
end.

