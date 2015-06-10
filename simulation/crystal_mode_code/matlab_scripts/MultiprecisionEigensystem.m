data = Import["C:\\Users\\ACKWinDesk\\Documents\\GitHub\\ultracold-ions\\matlab_scripts\\MannyLinearization.txt","Table"]; 
{values, vectors} = Eigensystem[SetPrecision[data, 32]]; 
Export["C:\\Users\\ACKWinDesk\\Documents\\GitHub\\ultracold-ions\\matlab_scripts\\MathematicaEigenvalues_real.txt", Re[values], "List"]; 
Export["C:\\Users\\ACKWinDesk\\Documents\\GitHub\\ultracold-ions\\matlab_scripts\\MathematicaEigenvectors_real.txt",Re[vectors], "Table"]; 
Export["C:\\Users\\ACKWinDesk\\Documents\\GitHub\\ultracold-ions\\matlab_scripts\\MathematicaEigenvalues_imag.txt", Im[values], "List"]; 
Export["C:\\Users\\ACKWinDesk\\Documents\\GitHub\\ultracold-ions\\matlab_scripts\\MathematicaEigenvectors_imag.txt", Im[vectors], "Table"];
Exit[]