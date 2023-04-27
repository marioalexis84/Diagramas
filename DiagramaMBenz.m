#!/home/alexis/Descargas/Mathematica/Executables/wolframscript --file


<<FeynCalc`


GravitonVertex[m1_,n1_,p1_,m2_,n2_,p2_,m3_,n3_,p3_]=Get["./GravitonVertex_1"];


FCSetScalarProducts[{SPD[k1,k2],SPD[k1,k3],SPD[k2,k3]},{0,0,0}]


A[a1_,b1_,k1_,a2_,b2_,k2_,a3_,b3_,k3_] = FullSimplify[ScalarProductExpand[GravitonVertex[a1,b1,k1,a2,b2,k2,a3,b3,k3]]];


FCClearScalarProducts[]


(*Esta funcion tiene el parametro "a" que toma valores 1,2,3 poniendo a cero el k1,k2 y k3 respectivamente*)


GravitonVertexR[a1_,b1_,k1_,a2_,b2_,k2_,a3_,b3_,k3_] = Simplify[Simplify[GravitonVertex[a1,b1,k1,a2,b2,k2,a3,b3,k3]]-A[a1,b1,k1,a2,b2,k2,a3,b3,k3];]


(*GravitonVertexR[a1,b1,k1,a2,b2,k2,a3,b3,k3]*)


(*Un ejemplo: poniendo a=3, mato los terminos con k3, reduciendose a 1/3 el numero de t\[EAcute]rminos. *)


GravitonPropagatorC[\[Mu]_,\[Nu]_,\[Alpha]_,\[Beta]_,k_,D_] = FullSimplify[1/2 I FeynAmpDenominator[PropagatorDenominator[Momentum[k,D],0]] (2*Pair[LorentzIndex[\[Alpha],D],LorentzIndex[\[Nu],D]] Pair[LorentzIndex[\[Beta],D],LorentzIndex[\[Mu],D]]-(2 Pair[LorentzIndex[\[Alpha],D],LorentzIndex[\[Beta],D]] Pair[LorentzIndex[\[Mu],D],LorentzIndex[\[Nu],D]])/(-2+D))]


GravitonVertexRsim[s_,a_,a1_,b1_,k1_,a2_,b2_,k2_,a3_,b3_,k3_]=If[s==0,Simplify[GravitonVertexR[a1,b1,(1-KroneckerDelta[a,1])*k1,a2,b2,(1-KroneckerDelta[a,2])*k2,a3,b3,(1-KroneckerDelta[a,3])*k3]],Simplify[A[a1,b1,(1-KroneckerDelta[a,1])*k1,a2,b2,(1-KroneckerDelta[a,2])*k2,a3,b3,(1-KroneckerDelta[a,3])*k3]]];


(*GravitonVertexRsim[1,1,J,R,k1,a2,b2,k2,a3,b3,k3]*)


(*Esta parte del diagrama corresponde a la union de los tres vertices en el borde, formando un circulo con tres patas salientes de momentos p1,p2,p3, con el doble indice aibi (i=1,2,3) correspondiente,*)


(*DiagramaPazfase1[t1_,t2_,t3_,s1_,s2_,s3_,a1_,b1_,p1_,a2_,b2_,p2_,a3_,b3_,p3_,k1_,k2_,k3_,D_]:=Simplify[Calc[GravitonVertexRsim[t1,s1,m1,n1,k1,a1,b1,p1,q1,r1,-k2]GravitonVertexRsim[t2,s2,m2,n2,k2,a2,b2,p2,q2,r2,-k3]GravitonVertexRsim[t3,s3,m3,n3,k3,a3,b3,p3,q3,r3,-k1]*GravitonPropagatorC[q1,r1,m2,n2,k2,D]*GravitonPropagatorC[q2,r2,m3,n3,k3,D]*GravitonPropagatorC[q3,r3,m1,n1,k1,D]]];*)


(*A continuacion uno el anterior con el vertice central. Los momentos p1 , p2 y p3 se pondran luego en terminos de k1,k2,k3.*)


(*DiagramaPazfase2[t1_,t2_,t3_,t4_,s1_,s2_,s3_,s4_,k1_,k2_,k3_,p1_,p2_,p3_,D_]:=Simplify[Calc[DiagramaPazfase1[t1,t2,t3,s1,s2,s3,a1,b1,p1,a2,b2,p2,a3,b3,p3,k1,k2,k3,D]*GravitonVertexRsim[t4,s4,a11,b11,p1,a22,b22,p2,a33,b33,p3]*GravitonPropagatorC[a1,b1,a11,b11,p1,D]GravitonPropagatorC[a2,b2,a22,b22,p2,D]GravitonPropagatorC[a3,b3,a33,b33,p3,D]]];*)


(*DiagramaMarioPaz[t1_,s1_,t2_,s2_,t3_,s3_,t4_,s4_,k1_,k2_,k3_,D_]:= Simplify[ScalarProductExpand[DiagramaPazfase2[t1,t2,t3,t4,s1,s2,s3,s4,k1,k2,k3,k2-k1,k3-k2,k1-k3,D]]]*)


u=Tuples[{{0,1},{1,1},{0,2},{1,2},{0,3},{1,3}},4];


(*DiagramaMarioPazmod[j_,k1_,k2_,k3_,D_]:=DiagramaMarioPaz[u[[j]][[1]][[1]],u[[j]][[1]][[2]],u[[j]][[2]][[1]],u[[j]][[2]][[2]],u[[j]][[3]][[1]],u[[j]][[3]][[2]],u[[j]][[4]][[1]],u[[j]][[4]][[2]],k1,k2,k3,D];*)


(*(*Una forma diferente de calcular lo mismo. Hay que meter el j cada vez. No evaluar antes de meter j*)*)


MBenzModeloLoopStep[t_, k1_, k2_, k3_, D_] := Module[
    {i, step, fase1,fase2,fase3,DiagramaMBenz,t1,t2,t3,t4,s1,s2,s3,s4,p1,p2,p3, res},

    step = 6;
    (*El item base t debe cumplir t=1 Mod(6)*)
    res = Catch[If[Mod[t-1, step] != 0,
        Throw["Numero de item base invalido: "<>ToString[t]],

        {{t1,s1},{t2,s2},{t3,s3},{t4,s4}} = u[[t]];
        Echo[ToString[t], "Inicio fase1: "];

        fase1 = Simplify[Calc[
            GravitonVertexRsim[t1,s1,m1,n1,k1,a1,b1,p1,q1,r1,-k2]
            * GravitonVertexRsim[t2,s2,m2,n2,k2,a2,b2,p2,q2,r2,-k3]
            * GravitonVertexRsim[t3,s3,m3,n3,k3,a3,b3,p3,q3,r3,-k1]
            * GravitonPropagatorC[q1,r1,m2,n2,k2,D]
            * GravitonPropagatorC[q2,r2,m3,n3,k3,D]
            * GravitonPropagatorC[q3,r3,m1,n1,k1,D]]];

        Echo[ToString[t], "Fin fase1: "];

        For[i=0, i<step, i++,
            {{t1,s1},{t2,s2},{t3,s3},{t4,s4}} = u[[t+i]];
            fase2 = Simplify[Calc[
                fase1
                * GravitonVertexRsim[t4,s4,a11,b11,p1,a22,b22,p2,a33,b33,p3]
                * GravitonPropagatorC[a1,b1,a11,b11,p1,D]
                * GravitonPropagatorC[a2,b2,a22,b22,p2,D]
                * GravitonPropagatorC[a3,b3,a33,b33,p3,D]]];

            Clear[fase1];
            Echo[ToString[t+i], "Fin fase2: "];

            fase3 = fase2/. {p1->k2-k1,p2->k3-k2,p3->k1-k3};
            DiagramaMBenz = Simplify[ScalarProductExpand[fase3]];
            Echo[ToString[t+i], "Fin fase3: "];

            Export["./TDiezPASO_"<>ToString[t+i]<>".txt", DiagramaMBenz];]
        ]];
    Echo[res];
];


(*For[i=1,i<1297,i++,t=DiagramaMarioPazmod[i,k1,k2,k3,D];Export["MBenz"<>ToString[i]<>".txt",t]]*)

(*F[i_]:=Export["./TDiez_"<>ToString[i]<>".txt", MBenznuevomodeloj[i,k1,k2,k3,D]];*)

KERNELS = 2
LaunchKernels[KERNELS]

ParallelTable[<<FeynCalc`;,{i,1,KERNELS}];

F[i_]:= MBenzModeloLoopStep[i,k1,k2,k3,D];

step = 5;

Parallelize[Map[F,Range[1, 36, step]]];

(*
Parallelize[Map[F,{1035, 1117, 1029, 1131, 1028}]];
*)



