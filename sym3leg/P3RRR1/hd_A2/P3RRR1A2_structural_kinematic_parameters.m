% Return Structural Kinematic Parameters of the parallel Robot 
% P3PPR1A1
%
% Output:
% NQJ_leg [1x1]
%   Anzahl der verwendeten Gelenkkoordinaten der symmetrischen Beinsegmente.
%   Von der Basis aus gezählt

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-03 13:01
% Revision: 73fedb8d52c8ee36c720354c027a250325027d46 (2019-01-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [NQJ_leg] = P3PPR1A1_structural_kinematic_parameters()

NQJ_leg = 2;
