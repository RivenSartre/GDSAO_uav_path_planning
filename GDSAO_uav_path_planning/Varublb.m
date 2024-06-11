function [VarMin,VarMax] = Varublb(model,nVar)

    VarMin.x=model.xmin;           
    VarMax.x=model.xmax;           
    VarMin.y=model.ymin;           
    VarMax.y=model.ymax;           
    VarMin.z=model.zmin;           
    VarMax.z=model.zmax;                 
    
    VarMax.r=2*norm(model.start-model.end)/nVar;           
    VarMin.r=0;
    
  
    AngleRange = pi/4; 
    VarMin.psi=-AngleRange;            
    VarMax.psi=AngleRange;          
    
    
    
    dirVector = model.end - model.start;
    phi0 = atan2(dirVector(2),dirVector(1));
    VarMin.phi=phi0 - AngleRange;           
    VarMax.phi=phi0 + AngleRange;           

end