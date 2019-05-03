% Return Structural Kinematic Parameters of the parallel Robot 
% P3RPP1A0
%
% Output:
% NQJ_leg [1x1]
%   Anzahl der verwendeten Gelenkkoordinaten der symmetrischen Beinsegmente.
%   Von der Basis aus gezählt

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [NQJ_leg] = P3RPP1A0_structural_kinematic_parameters()

NQJ_leg = 3;
