% Return Structural Kinematic Parameters of the parallel Robot 
% P3RPRR1G2A0
%
% Output:
% NQJ_leg [1x1]
%   Anzahl der verwendeten Gelenkkoordinaten der symmetrischen Beinsegmente.
%   Von der Basis aus gezählt

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [NQJ_leg] = P3RPRR1G2A0_structural_kinematic_parameters()

NQJ_leg = 3;
