# Wolfram-Processing-Language-
Mathematic functions in WPL - Matrices, Eigenvalues, Plotting 3D Definitive Matrices, 1st and 2nd order derivatives, profit maximization, etc. 
A positive definite matrix must have positive diagonal entires. A is positive definite if and only if all its eigenvalues are positive. A is negative definite if and only if all its eigenvalues are negative. In[36]:= mA = {{1, 1, 0}, {1, 4, 2}, {0, 2, 3}} mB = {{5,-6,-6}, {-1, 4, 2}, {3,-6,-4}} mC = {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}} Out[36]= {{1, 1, 0}, {1, 4, 2}, {0, 2, 3}} Out[37]= {{5,-6,-6}, {-1, 4, 2}, {3,-6,-4}} Out[38]= {{1, 0, 1}, {0, 1, 1}, {1, 1, 2}}
In[45]:= NegativeDefiniteMatrixQ[mA] NegativeDefiniteMatrixQ[mB] NegativeDefiniteMatrixQ[mC] 
Out[45]= False 
Out[46]= False 
Out[47]= False
In[48]:= NegativeSemidefiniteMatrixQ[mA] NegativeSemidefiniteMatrixQ[mB] NegativeSemidefiniteMatrixQ[mC]
Out[48]= False 
Out[49]= False 
Out[50]= False
PositiveSemidefiniteMatrixQ[mA] PositiveSemidefiniteMatrixQ[mB] PositiveSemidefiniteMatrixQ[mC]
Out[51]= True 
Out[52]= False 
Out[53]= True 
In[54]:= PositiveDefiniteMatrixQ[mA] PositiveDefiniteMatrixQ[mB] PositiveDefiniteMatrixQ[mC] 
Out[54]= True 
Out[55]= False 
Out[56]= False

Plotting 3D graphs 
In[59]:= Clear[mA, mB]
In[183]:= mA = {{4, 0}, {0, 3}}; mAx = mA-λ 0 0 λ y =Det[mAx] mB = {{3,-2}, {-2, 7}}; mBx = mB-(λ 0 0 λ) y1 = Det[mBx] Solve[(4-λ)(3-λ) == 0, λ] Solve[(3-λ)(7-λ) == 0, λ] mX = (x1 x2); Transpose[mX].mA.mX mX1 = (x1 x2 ); Transpose[mX1].mB.mX1 
In[62]:= PositiveDefiniteMatrixQ[mA] PositiveDefiniteMatrixQ[mB]
In[318]:= Plot3D[4 x1^2+3x2^2, {x1,-6, 6}, {x2,-6, 6}, PlotStyle -> Purple]

In[317]:= Plot3D[x1 (3 x1-2x2)+x2(-2x1+7x2), {x1,-6, 6}, {x2,-6, 6}, PlotStyle -> Purple]

Creating Hessian Matrices and Eigenvalues 
In[86]:= q1 = {x1, x2, x3}  6*x1^2+25*x2^2+9*x3^2-60*x2*x3+40*x1*x3-6*x1*x2; q2 = {x1, x2, x3} 9*x2^2+9*x3^2+6*x1*x2+10*x1*x3; 
In[88]:= (*Computing Hessian in WL*) D[q[x1, x2], {{x1, x2}, 2}]
In[91]:= Out[91]= H1 = D[q1[x1, x2], {{x1, x2}, 2}] H2 = D[q2[x1, x2], {{x1, x2}, 2}]
In[94]:= A1 = H1/2 A2 = H2/2
In[96]:= Eigenvalues[A1] Eigenvalues[A2]

First and Second Order Derivatives

In[321]:= f={x1,x2}|-> a1*x1+a2*x2 l={x1,x2,x3}|->{a1,a2,a3}.{x1,x2,x3} h={x1,x2}|->x1^3*x2^4 k={x1,x2,x3}|->x1^2*x2^4*x3^5
In[373]:= Clear[x1,x2] f={x1,x2}|->a1*x1+a2*x2; D[f[x1,x2],{{x1,x2}}]//MatrixForm D[f[x1,x2],{{x1,x2},2}]//MatrixForm

Rule of Sarrus 
I n[1]= labels = {{"a11", "a12", "a13"}, {"a21", "a22", "a23"}, {"a31", "a32", "a33"}, {"a11", "a12", "a13"}, {"a21", "a22", "a23"}}; textElements = Table[Text[Style[labels〚i, j〛, 16], {j,-i}], {i, Length[labels]}, {j, Length[labels〚i〛]}]; textElements = Flatten[textElements]; downDiagonals = {Blue, Thick, Line[{{1,-1}, {2,-2}, {3,-3}}], Line[{{2,-1}, {3,-2}, {4,-3}}], Line[{{3,-1}, {4,-2}, {5,-3}}]}; upDiagonals = {Red, Thick, Line[{{3,-1}, {2,-2}, {1,-3}}], Line[{{4,-1}, {3,-2}, {2,-3}}], Line[{{5,-1}, {4,-2}, {3,-3}}]}; Graphics[{downDiagonals, upDiagonals, Black, textElements}, Axes -> False, PlotRange -> {{0.5, 5.5}, {-3.5,-0.5}}, AspectRatio -> 1, ImageSize -> Small, PlotLabel -> "Rule of Sarrus: Visual Representation"]

In this problem we are considering the following textbook model of money supply determination.
High powered money stock: H=C+R money stock: M=C+D (a) First we are asked to express this equation system as a matrix equation in the form A.x with endogenous vecto x =<C,R,D,M> C=CD R=rD M=C+D H=C+R
Use matrix inverse tosolve for CRMD x=A-1b
ΔM /ΔH =Δ(C+D)/ (C+R) = ΔD / (Δ(C+R))
Code the Results 
I n[2]:= c =0.1; r =0.02; H =800; A ={{1, 0,-c, 0}, {0, 1,-r, 0}, {1, 1, 0, 0}, {1, 0, 1,-1}}; b ={0, 0, H, 0}; xInverse = Inverse[A].b; xSolve = LinearSolve[A, b]; vars = {"C (Currency)", "R (Reserves)", "D (Deposits)", "M (Money Supply)"}; Print["Solution using Inverse[A]:\n"]; Do[Print[vars〚i〛, " = ", NumberForm[xInverse〚i〛, {6, 2}], " billion USD"], {i, 4}]; Print["\nSolution using LinearSolve:\n"]; Do[Print[vars〚i〛, " = ", NumberForm[xSolve〚i〛, {6, 2}], " billion USD"], {i, 4}]; Print["\nMoney Multiplier (M / H) = ", NumberForm[xSolve〚4〛/H, {4, 2}]];
I n[3]= Solution to the Money Supply System:Currency held by public (C) = 727.27billion USD Reserves held by banks (R) = 145.45billion USD Deposits (D) = 7272.73 billion USD Money supply (M) = 8000.00 billion USD Money multiplier (M/H) = 10.00

Graphing Transformed Matrices
In[4]:= transformSquare[T_]:= Module[{square,transformed},square=Polygon[{{0,0},{1,0},{1,1},{0,1}}]; (*unit square*)transformed=GeometricTransformation[square,T]; Graphics[{Style[square,LightGray],Style[transformed,Blue,Thick]}, Axes->True,PlotRange->{{-2,2},{-2,2}},ImageSize->200]] 

In[5]:= matrices={{{1.5,0},{0,1}},{{1,0},{0,1.5}}, {{1.5,0},{0,1.5}},{{-1,0},{0,1}},{{-1,0},{0,-1}}}; GraphicsRow[transformSquare/@matrices,Spacings->Scaled[0.5]]
In[6 ]:= θ = 60 Degree; Cmat = RotationMatrix[θ]; MatrixQ[Cmat] (*Should return True*) transformSquare[m_] := GeometricTransformation[Polygon[{{0, 0}, {1, 0}, {1, 1}, {0, 1}}], m] Graphics[{LightGray, Polygon[{{0, 0}, {1, 0}, {1, 1}, {0, 1}}], Red, transformSquare[Cmat]}, Axes -> True, PlotRange -> {{-2, 2}, {-2, 2}}, ImageSize ->300] Det[A] (* =1.1×0.9=0.99*) Det[B] (* =1.2×0.9=1.08*) Det[C] (* =1.0—pure rotation*) Det[CA] (* =Det[C]*Det[A]=0.99*) Det[CB] (* =Det[C]*Det[B]=1.08*)
