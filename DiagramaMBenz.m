#!/home/alexis/Descargas/Mathematica/Executables/wolframscript --file


<<FeynCalc`


GravitonVertex[m1_,n1_,p1_,m2_,n2_,p2_,m3_,n3_,p3_] = Get["./GravitonVertex_1"];


FCSetScalarProducts[{SPD[k1,k2],SPD[k1,k3],SPD[k2,k3]},{0,0,0}]


A[a1_,b1_,k1_,a2_,b2_,k2_,a3_,b3_,k3_] = FullSimplify[ScalarProductExpand[GravitonVertex[a1,b1,k1,a2,b2,k2,a3,b3,k3]]];


FCClearScalarProducts[]


(*Esta funcion tiene el parametro "a" que toma valores 1,2,3 poniendo a cero el k1,k2 y k3 respectivamente*)


GravitonVertexR[a1_,b1_,k1_,a2_,b2_,k2_,a3_,b3_,k3_] = Simplify[Simplify[
    GravitonVertex[a1,b1,k1,a2,b2,k2,a3,b3,k3]]-A[a1,b1,k1,a2,b2,k2,a3,b3,k3]];


(*Un ejemplo: poniendo a=3, mato los terminos con k3, reduciendose a 1/3 el numero de t\[EAcute]rminos. *)


GravitonPropagatorC[\[Mu]_,\[Nu]_,\[Alpha]_,\[Beta]_,k_,D_] = FullSimplify[
    1/2*I*FeynAmpDenominator[PropagatorDenominator[Momentum[k,D],0]] * 
    (2*Pair[LorentzIndex[\[Alpha],D],LorentzIndex[\[Nu],D]] *
     Pair[LorentzIndex[\[Beta],D],LorentzIndex[\[Mu],D]] - 
     (2*Pair[LorentzIndex[\[Alpha],D],LorentzIndex[\[Beta],D]] *
     Pair[LorentzIndex[\[Mu],D],LorentzIndex[\[Nu],D]])/(-2+D))];


GravitonVertexRsim[s_,a_,a1_,b1_,k1_,a2_,b2_,k2_,a3_,b3_,k3_] = If[
    s == 0,
    Simplify[GravitonVertexR[a1, b1, (1 - KroneckerDelta[a,1])*k1, a2, b2, (1-KroneckerDelta[a,2])*k2, a3, b3, (1-KroneckerDelta[a,3])*k3]],
    Simplify[A[a1, b1, (1-KroneckerDelta[a,1])*k1, a2, b2, (1-KroneckerDelta[a,2])*k2, a3, b3, (1-KroneckerDelta[a,3])*k3]]];


u = Tuples[{{0,1}, {1,1}, {0,2}, {1,2}, {0,3}, {1,3}}, 4];

(*
MBenznuevomodeloj[t_, k1_, k2_, k3_, D_] := Module[
   {fase1, fase2, fase3, DiagramaMBenz, t1, t2, t3, t4, s1, s2, s3, 
    s4, p1, p2, p3},
   {{t1, s1}, {t2, s2}, {t3, s3}, {t4, s4}} = u[[t]];
   fase1 = Simplify[Calc[GravitonVertexRsim[t1, s1, m1, n1, k1, a1, b1, p1, q1, r1, -k2]* GravitonVertexRsim[t2, s2, m2, n2, k2, a2, b2, p2, q2, r2, -k3]* GravitonVertexRsim[t3, s3, m3, n3, k3, a3, b3, p3, q3, r3, -k1]* GravitonPropagatorC[q1, r1, m2, n2, k2, D]* GravitonPropagatorC[q2, r2, m3, n3, k3, D]* GravitonPropagatorC[q3, r3, m1, n1, k1, D]]];
   fase2 = Simplify[Calc[      fase1*       GravitonVertexRsim[t4, s4, a11, b11, p1, a22, b22, p2, a33,         b33, p3]*       GravitonPropagatorC[a1, b1, a11, b11, p1, D]*       GravitonPropagatorC[a2, b2, a22, b22, p2, D]*       GravitonPropagatorC[a3, b3, a33, b33, p3, D]]];
   Clear[fase1];
   fase3 = fase2 /. {p1 -> k2 - k1, p2 -> k3 - k2, p3 -> k1 - k3};
   DiagramaMBenz = Simplify[ScalarProductExpand[fase3]]; 
   DiagramaMBenz];
*)

MBenzModelo[t_, k1_, k2_, k3_, D_] := Module[
    {i, j, suffixFase0, suffixFase1,
     fase0, fase1, fase2, fase3, DiagramaMBenz,
     t1, t2, t3, t4, s1, s2, s3, s4, p1=p1, p2=p2, p3=p3},

    tInicial = AbsoluteTime[];

    i = Quotient[t, 36, 1]*36;
    suffixFase0 = ToString[i+1]<>"-"<>ToString[i+36];
    j = Quotient[t, 6, 1]*6;
    suffixFase1 = ToString[j+1]<>"-"<>ToString[j+6];

    {{t1,s1}, {t2,s2}} = u[[t, 1;;2]];

    Print["Inicio fase0 - Item "<>ToString[t]];

    Check[
        fase0 = Get["./Fases/T10Fase0_"<>suffixFase0<>".txt"];,
        fase0 = Simplify[Calc[
            GravitonVertexRsim[t1, s1, m1, n1, k1, a1, b1, p1, q1, r1, -k2] *
            GravitonVertexRsim[t2, s2, m2, n2, k2, a2, b2, p2, q2, r2, -k3] *
            GravitonPropagatorC[q1, r1, m2, n2, k2, D]]];

            Export["./Fases/T10Fase0_"<>suffixFase0<>".txt", fase0];,
        Get::noopen];

    Print["Fin fase0 - Item: "<>ToString[t]<>":\tTiempo: "<>ToString[AbsoluteTime[]-tInicial]];

    {t3, s3} = u[[t, -2]];

    suffixSubFase1  = ToString[t3]<>"_"<>ToString[s3];
    Check[
        subFase1 = Get["./Fases/SubFases/T10SubFase1_"<>suffixSubFase1<>".txt"];,
        subFase1 = Simplify[Calc[
            GravitonVertexRsim[t3, s3, m3, n3, k3, a3, b3, p3, q3, r3, -k1] *
            GravitonPropagatorC[q2, r2, m3, n3, k3, D] *
            GravitonPropagatorC[q3, r3, m1, n1, k1, D]]];
            Export["./Fases/SubFases/T10SubFase1_"<>suffixSubFase1<>".txt", subFase1];,
        Get::noopen];

    Print["Fin subFase1 - Item: "<>ToString[t]<>":\tTiempo: "<>ToString[AbsoluteTime[]-tInicial]];

    Check[
        fase1 = Get["./Fases/T10Fase1_"<>suffixFase1<>".txt"];,
        fase1 = Simplify[Calc[fase0 * subFase1]];
            Export["./Fases/T10Fase1_"<>suffixFase1<>".txt", fase1];,
        Get::noopen];

    Print["Fin fase1 - Item: "<>ToString[t]<>":\tTiempo: "<>ToString[AbsoluteTime[]-tInicial]];

    {t4, s4} = u[[t, -1]];

    suffixSubFase2  = ToString[t4]<>"_"<>ToString[s4];
    Check[
        subFase2 = Get["./Fases/SubFases/T10SubFase2_"<>suffixSubFase2<>".txt"];,
        subFase2 = Simplify[Calc[
            GravitonVertexRsim[t4, s4, a11, b11, p1, a22, b22, p2, a33, b33, p3] *
            GravitonPropagatorC[a1, b1, a11, b11, p1, D] *
            GravitonPropagatorC[a2, b2, a22, b22, p2, D] *
            GravitonPropagatorC[a3, b3, a33, b33, p3, D]]];
            Export["./Fases/SubFases/T10SubFase2_"<>suffixSubFase2<>".txt", subFase2];,
        Get::noopen];

    Print["Fin subFase2 - Item: "<>ToString[t]<>":\tTiempo: "<>ToString[AbsoluteTime[]-tInicial]];

    fase2 = Simplify[Calc[fase1 * subFase2]];

    Print["Fin fase2 - Item: "<>ToString[t]<>":\tTiempo: "<>ToString[AbsoluteTime[]-tInicial]];

    fase3 = fase2 /. {p1 -> k2 - k1, p2 -> k3 - k2, p3 -> k1 - k3};
    DiagramaMBenz = Simplify[ScalarProductExpand[fase3]]; 

    Print["Fin fase3 - Item: "<>ToString[t]<>":\tTiempo: "<>ToString[AbsoluteTime[]-tInicial]];

    DiagramaMBenz
];


(*For[i=1,i<1297,i++,t=DiagramaMarioPazmod[i,k1,k2,k3,D];Export["MBenz"<>ToString[i]<>".txt",t]]*)


KERNELS = 2
LaunchKernels[KERNELS]

ParallelTable[<<FeynCalc`;,{i, 1, KERNELS}];

F[i_]:= Export["./TDiezTestFasesMODSubFases_"<>ToString[i]<>".txt", MBenzModelo[i, k1, k2, k3, D]];

Parallelize[Map[F,Range[1148, 1152]]];

(*
Parallelize[Map[F,{3}]];
*)
