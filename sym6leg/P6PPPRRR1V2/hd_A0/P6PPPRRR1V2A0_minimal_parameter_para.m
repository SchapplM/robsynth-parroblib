% Return the minimum parameter vector for
% P6PPPRRR1V2A0
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d4,theta1,theta2,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% 
% Output:
% MPV [13x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-16 19:43
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = P6PPPRRR1V2A0_minimal_parameter_para(pkin, m, mrSges, Ifges, koppelP)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6),zeros(6,3)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PPPRRR1V2A0_minimal_parameter_para: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6PPPRRR1V2A0_minimal_parameter_para: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6PPPRRR1V2A0_minimal_parameter_para: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P6PPPRRR1V2A0_minimal_parameter_para: Ifges has to be [4x6] (double)'); 
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PPPRRR1V2A0_minimal_parameter_para: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From minimal_parameter_parrob_matlab.m
t1 = [m(1); m(2); m(3); Ifges(4,1); Ifges(4,4); Ifges(4,5); Ifges(4,2); Ifges(4,6); Ifges(4,3); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4);];
MPV  = t1;
