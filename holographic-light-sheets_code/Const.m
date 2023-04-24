classdef Const
    properties (Constant)
      e0 = 8.8541878176e-12; %[F/m] permissividade no espaço livre 
      u0 = 4*pi*1e-7; %[H/m] permeabilidade no espaço livre
      c0 = 1/sqrt(Const.e0*Const.u0); %[m/s] velocidade da luz no espaço livre
      eta0 = sqrt(Const.u0/Const.e0); %[Ohm] impedancia caracterisitca do espaço livre      
    end
end