$Title Hjemmeopgave 4.1

*suppressing output
$offsymxref offsymlist
option limcol =0, limrow=0;
option solprint=off;
*------------------------------------------------------------------------*
*Objekter og dimensioner
*------------------------------------------------------------------------*
sets t /t0*t100/;

Variables
*Endogene variable
K(t)     "Kapital"
Y(t)     "Output"
L(t)     "Arbejdsudbud"
I(t)     "Investeringer"
C(t)     "Privatforbrug"
V(t)     "Virksomhedernes v�rdi"
uc(t)    "User-cost for kapital"
W(t)     "L�nnen"
r(t)     "Renten"
*Eksogene variable
P(t)     "Prisen i indlandet"
delta(t) "Afskrivningsraten"
theta(t) "Arbejdskraftsproduktivitet"
rho      "Risikoaversion"
eta      "Tilbagediskonteringsrate"
myK(t)      "Produktionsv�gt kapital"
myL(t)      "Produktionsv�gt arbejdskraft"
E        "Elasticitet i produktionsfkt"
g        "V�kst"
* Gem basisværdier i t99
K_base    "Kapital t99, basis"
Y_base    "Output t99, basis"
C_base    "Forbrug t99, basis"
I_base    "Investering t99, basis"
L_base    "Arbejde t99, basis"
V_base    "Virksomhedsværdi t99, basis"
r_base    "Rente t99, basis"
uc_base   "User cost t99, basis"
P_base    "Pris t99, basis"
delta_base "Afskrivning t99, basis"
* Gem t99-resultater til tabel
C_t99_theta1  "Forbrug t99, θ-stød i t1"
K_t99_theta1  "Kapital t99, θ-stød i t1"
Y_t99_theta1  "Output t99, θ-stød i t1"
 C_t99_theta5  "Forbrug t99, θ-stød i t5"
    K_t99_theta5  "Kapital t99, θ-stød i t5"
    Y_t99_theta5  "Output t99, θ-stød i t5"
    C_t99_tsunami  "Forbrug t99, tsunami"
    K_t99_tsunami  "Kapital t99, tsunami"
    Y_t99_tsunami  "Output t99, tsunami"
    C_t99_delta  "Forbrug t99, lavere delta"
    K_t99_delta  "Kapital t99, lavere delta"
    Y_t99_delta  "Output t99, lavere delta"
    
;



Equations
E_K(t)
E_W(t)
E_P(t)
E_uc(t)
E_I(t)
E_C(t)
E_Y(t)
E_V(t)
E_Vterm
E_ucterm
E_Cterm
;

E_K(t)    $ (ord(t) gt 1)..                         K(t-1)/(1+g)  =e= MyK(t) * (uc(t-1)/P(t))**(-E) * Y(t);
E_W(t)    $ (ord(t) gt 1)..                         theta(t)*L(t) =e= MyL(t) * ((W(t)/theta(t))/P(t))**(-E)*Y(t);
E_P(t)    $ (ord(t) gt 1)..                         P(t)*Y(t)     =e= (uc(t-1)*K(t-1))/(1+g) + W(t)*L(t);
E_uc(t)   $ (ord(t) lt card(t))..                   uc(t)         =e= r(t+1) + delta(t);
E_I(t)    $ (ord(t) gt 1)..                         K(t)          =e= (1-delta(t))*K(t-1)/(1+g) + I(t);
E_C(t)    $ (ord(t) gt 1 and ord(t) lt card(t))..   C(t+1)*(1+g)  =e= ((1+r(t+1))/(1+eta))**(1/rho)*C(t);
E_Y(t)    $ (ord(t) gt 1)..                         Y(t)          =e= C(t) + I(t);
E_V(t)    $ (ord(t) gt 1 and ord(t) lt card(t))..   V(t+1)*(1+g)  =e= (1+r(t+1))*V(t) - (P(t+1)*Y(t+1)*(1+g) - w(t+1)*L(t+1)*(1+g) - P(t+1)*I(t+1)*(1+g));
E_Cterm..                                           C('t100')     =e= C('t99');
E_Vterm..                                           V('t100')     =e= V('t99');
E_ucterm..                                          uc('t100')    =e= uc('t99');

Model Ramsey /ALL/;
*------------------------------------------------------------------------*
* Data
*------------------------------------------------------------------------*

*$include
Sets
j "Input" /
PS   "Privat sektor"
Lon  "Lonsum"
Rest "restindkomst"
/
o "Output" /
PS   "Privat sektor"
C    "forbrug"
I    "Investeringer" /
;
Table IO(j,o) "Input-output-tabel"
        PS     C    I
PS       0   800  200
lon    700     0    0
Rest   300     0    0
;
   

* Eksogene parametre
E.fx        = 0.7;
theta.fx(t) = 1;
g.fx        = 0.02;
rho.fx      = 2 ;

*Antagelse
p.fx(t)     = 1;
w.l(t)      = 1;
r.l(t)      = 0.05;

* Initialisering
L.fx(t)     = 700;
C.l(t)      = 800;
I.l(t)      = 200;
Y.l(t)      = 1000;

* Kalibrering
eta.fx      = 0.0092272203;
delta.fx(t) = 2*r.l(t) + (-3*g.l);

* Initialisering
K.fx('t0')  =3400;
K.l(t)      =3400;
uc.l(t)     = 0.09;
V.l(t)      = 3400;

* Kalibrering
MyL.fx(t)      = L.l(t)/Y.l(t);
myK.fx(t)      = (K.l(t)/(1+g.l)/(Y.l(t)))*(r.l(t)+delta.l(t))**(E.l);

display 
K.l,    
Y.l,
L.l,
I.l,
C.l,
V.l,
uc.l,
W.l,
r.l,
P.l,
delta.l,
theta.l,
rho.l,
eta.l,
myK.l,
myL.l,
E.l,
g.l;

Solve Ramsey using CNS;

K_base.l     = K.l('t99');
Y_base.l     = Y.l('t99');
C_base.l     = C.l('t99');
I_base.l     = I.l('t99');
L_base.l     = L.l('t99');
V_base.l     = V.l('t99');
r_base.l     = r.l('t99');
uc_base.l    = uc.l('t99');
P_base.l     = P.l('t99');
delta_base.l = delta.l('t99');

* ---------- Stød 4.1a: produktivitetsstød i t1 ----------

theta.fx('t1') = 1.10;


Solve Ramsey using CNS;

C_t99_theta1.l = C.l('t99');
K_t99_theta1.l = K.l('t99');
Y_t99_theta1.l = Y.l('t99');

display 
K.l,    
Y.l,
L.l,
I.l,
C.l,
V.l,
uc.l,
W.l,
r.l,
P.l,
delta.l,
theta.l,
rho.l,
eta.l,
myK.l,
myL.l,
E.l,
g.l,
C_t99_theta1.l,
K_t99_theta1.l,
Y_t99_theta1.l
;

*-------------------------------------------------------------*
* Eksporter tidsserier til Excel (CSV-fil)
*-------------------------------------------------------------*

File out /"results_theta1.csv"/;
put out;

* Skriv kolonne-navne
put "t,Y,C,K" /;

* Loop gennem alle tidspunkter
loop(t,
    put t.tl:0 ","
        Y.l(t):12:6 ","
        C.l(t):12:6 ","
        K.l(t):12:6 /
);

putclose out;
