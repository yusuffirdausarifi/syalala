%%Fungsi 1 Metode Bagi Dua
f1 = @(x) 5*x^5 + 4*x^4 - 7*x^3 + 10*x^2 + 1;
a = -2.1; % Setelah melihat perilaku grafik fungsi 1 menggunakan grapher, diperoleh bahwa
b = -2; % salah satu akar f1 terletak pada range [-2.1,-2]  
max_iter = 100;
eps = 1e-16;

if f1(a)*f1(b) > 0
    fprintf("Tidak punya akar di interval [a,b]")
    return
end

if f1(a) == 0
    fprintf("a merupakan salah satu solusi")
    return
elseif f1(b) == 0
    fprintf("b merupakan salah satu solusi")
    return
end

for i = 1 : max_iter
    c = (a+b)/2;
    if f1(a)*f1(c) < 0
        b = c;
    else
        a = c;
    end
if abs(f1(a)) < eps
    break
end
end
fprintf("Solusi :%f",a)
%% Fungsi 1 Metode Posisi Palsu
f1 = @(x) 5*x^5 + 4*x^4 - 7*x^3 + 10*x^2 + 1;
a = -2.1;% Setelah melihat perilaku grafik fungsi 1 menggunakan grapher, diperoleh bahwa
b = -2; % salah satu akar f1 terletak pada range [-2.1,-2]
max_iter = 100;
eps = 1e-16;

if f1(a)*f1(b) > 0
    fprintf("Tidak punya akar di interval [a,b]")
    return
end

if f1(a) == 0
    fprintf("a merupakan salah satu solusi")
    return
elseif f1(b) == 0
    fprintf("b merupakan salah satu solusi")
    return
end

clama = 2*b-a;
for i = 1 : max_iter
    c = a - f1(a)*(b-a)/(f1(b)-f1(a));
    if f1(a)*f1(c) < 0
        b = c;
    else
        a = c;
    end
if abs(c-clama) < eps
    break
end
end
fprintf("Solusi :%f",a)
%% Fungsi 1 Metode Newton Rhapson
%syms x
% f = 5*x^5 + 4*x^4 - 7*x^3 + 10*x^2 + 1;
% df = diff(f)
f1 = @(x) 5*x^5 + 4*x^4 - 7*x^3 + 10*x^2 + 1;
df1 = @(x) 25*x^4 + 16*x^3 - 21*x^2 + 20*x;
x = -2.1; %Tebakan awal dipilih -2.1 karena titik tersebut merupakan salah satu titik
%terdekat dari grafik saat memotong sumbu x
for i = 1 : max_iter
    xbaru = x - f1(x)/df1(x);
    delta = abs(xbaru-x)/abs(xbaru);
    if delta < eps
        break
    end
end
fprintf("Solusi :%f",xbaru)
%% Fungsi 2 Metode Bagi Dua
f2 = @(x) exp(2*x^2 - 2) + log(x);
a = 0.5;% Setelah melihat perilaku grafik fungsi 2 menggunakan grapher, diperoleh bahwa
b = 1; % salah satu akar f1 terletak pada range [0.5,1]
max_iter = 100;
eps = 1e-16;

if f2(a)*f2(b) > 0
    fprintf("Tidak punya akar di interval [a,b]")
    return
end

if f2(a) == 0
    fprintf("a merupakan salah satu solusi")
    return
elseif f2(b) == 0
    fprintf("b merupakan salah satu solusi")
    return
end

for i = 1 : max_iter
    c = (a+b)/2;
    if f2(a)*f2(c) < 0
        b = c;
    else
        a = c;
    end
if abs(f2(a)) < eps
    break
end
end
fprintf("Solusi :%f",a)
%% Fungsi 2 Metode Posisi Palsu
f2 = @(x) exp(2*x^2-2)+log(x);
a = 0.5;% Setelah melihat perilaku grafik fungsi 2 menggunakan grapher, diperoleh bahwa
b = 1; % salah satu akar f1 terletak pada range [0.5,1]
max_iter = 100;
eps = 1e-16;

if f2(a)*f2(b) > 0
    fprintf("Tidak punya akar di interval [a,b]")
    return
end

if f2(a) == 0
    fprintf("a merupakan salah satu solusi")
    return
elseif f2(b) == 0
    fprintf("b merupakan salah satu solusi")
    return
end

clama = 2*b-a;
for i = 1 : max_iter
    c = a - f2(a)*(b-a)/(f2(b)-f2(a));
    if f2(a)*f2(c) < 0
        b = c;
    else
        a = c;
    end
if abs(c-clama) < eps
    break
end
end
fprintf("Solusi :%f",a)
%% Fungsi 2 Metode Newton Rhapson
%syms x
% f = exp(2*x^2-2)+log(x);
% df = diff(f)
f1 = @(x) exp(2*x^2-2)+log(x);
df1 = @(x) 4*x*exp(2*x^2) + 1/x;
x = 1; %Dilihat dari perilaku grafik nya, f2 memotong sumbu x disekitar x=1

for i = 1 : max_iter
    xbaru = x - f2(x)/df2(x);
    delta = abs(xbaru-x)/abs(xbaru);
    if delta < eps
        break
    end
end
fprintf("Solusi :%f",xbaru)
%% maap be, jujur masih kagok.. gmn declare function,
%Jadinya, semua soal belum pake function..