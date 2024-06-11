function [VelMin,VelMax] = Velublb(VarMax,VarMin)

alpha=0.5;
VelMax.r=alpha*(VarMax.r-VarMin.r);    
VelMin.r=-VelMax.r;                    
VelMax.psi=alpha*(VarMax.psi-VarMin.psi);    
VelMin.psi=-VelMax.psi;                    
VelMax.phi=alpha*(VarMax.phi-VarMin.phi);    
VelMin.phi=-VelMax.phi; 

end